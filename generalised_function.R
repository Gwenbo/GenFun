##### Code for generalised fitness function 

######****************************************************************************************** Main generalised function code  
## Function to update the mean fitness at each timestep
meanfit_vars<-function(M,L,acq,trans,p,react,treat,nA,nL,n,type,propF){
  # Needs M: Matrix of distribution  
  #       L: Matrix of distribution of latent
  #       acq: number of acquisitions of resistant strains
  #       trans: number of transmissions of resistant strains
  #       p: proportion of transmissions that progress immediately to active TB
  #       react: number of reactivations
  #       treat: number of treatments of active resistant TB (that become latent resistant TB)
  #       nA: number with active TB
  #       nL: number with latent TB 
  #       n: number of fitness levels
  #       type: changes distribution of stochastic sampling  
  #       propF: proportion go to rel. fit = 1 (in 2 strain example)
  
  # new vector = updated distribution of fitness in active cases
  new <- matrix(0,n,1)
  rownames(new) <- seq(1/n,1,1/n)
  
  # Vector of fitness costs
  v<-seq(1/n,1,1/n)
  
  # Which column of M? last that is non-zero plus one for this timestep
  if(length(which(colSums(M)>0))){tt<-tail(which(colSums(M)>0),n=1) + 1}else{tt=2}
  
  #************************************************ At initial time point
  if(acq==0 && trans==0){meanfit = 0; return(meanfit)
  #************************************************ Otherwise look at distribution  
  }else{
    # Scale events
    a<-acq/(acq+trans+react);  b<-trans/(acq+trans+react);  c<-react/(acq+trans+react);
    
    #*** How are the acquisitions distributed? 
    # Can do a complex function here - sample from initial distribution of fitness costs
    # For now just 99% unfit, 1% fit
    if(type=="orig"){new["0.6",] <- (1-propF)*a; new["1",] <- propF*a }
    if(type=="normal"){h <- hist(rnorm(acq,0.7,0.1), breaks = c(v[1]-(v[2]-v[1]),v), plot = F,right = F); new<-a*h$counts/(acq)} # not tested yet
    if(type=="det"){new["0.8",] <- (1-propF)*a; new["0.9",] <- propF*a } # not tested yet
    
    #*** How are the transmissions distributed?
    if(trans>0){pastmean=sum(M[,tt-1]*v)}else{ pastmean = 1 }
    new = new + b * M[,tt-1] * v / pastmean
    
    #*** How are the reactivations distributed?
    if(react>0){
    new = new + c * L[,tt-1]}
    
    #*** Assign and update M
    if(nA>0){M[,tt] = M[,tt-1]*( nA/(nA+acq+trans+react) ) + new * (acq+trans+react)/(nA+acq+trans+react) 
    }else{M[,tt] =  new * (acq+trans+react)/(nA+acq+trans+react) }
    
    
    #*** Assign and update L 
    # At start have no latents.
    transL<-trans*(1/p)*(1-p)
    if((nL+transL)>0){L[,tt] = (L[,tt-1]*nL + treat*M[,tt-1] + transL*M[,tt-1]* v / pastmean )/(nL+treat+transL)}else{L[,tt]=L[,tt-1]}
    
    #*** Checks
    if(a+b+c > 1.01){print(c("Error in proportions"))}
    if(any(colSums(M)>1.01)){print(c(sum(M[,which(colSums(M)>1)]),"colsums",colSums(M),"ERROR in colsums of M"));break}
    if(any(colSums(L)>1.01)){print(c(sum(L[,which(colSums(L)>1)]),"colsums",colSums(L),"ERROR in colsums of L"));break}
    
    #*** Calculate mean fitness
    meanfit = sum(M[,tt]*v)
  }
  
  #*** Output mean fitness, new distributions of active and latent 
  return(list(meanfit=meanfit,M=M,L=L))
}

                        
######****************************************************************************************** Difference equation set up taking in above
## Function which calculates population dynamics over time
sourya_funcf_mean_vars=function(endp,home,vary,initial,M0,L0,type,propF){
  # Needs endp: length of simulation
  #       home: location
  #       vary: key parameters can change (epsilon, f and treatment differential (krdiff))
  #       initial: initial population (150/100K incidence etc) 
  #       M0,L0: initial distributions, gives information on number of fitness levels 
  #       type: type for acquisition distribution (options: "orig" "normal")
  #       propF: proportion of new acquisitons that are fit
  
  #*** Parameter assignment
  # Might not need if globally assigned
  setwd(home) # might have para in different place for different models 
  para<-read.csv("para.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE); 
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  N<-1000000; dt<-0.01
  # Assign vary after to overwrite (? not checked overruled global assignment for meanfit function)
  # Parameters assigned by length of vary so for baseline just eps = 0
  vary_n = c("eps", "f", "krdiff")
  for(i in 1:length(vary)){assign(vary_n[i],vary[i])}
  assign("kr",ks + krdiff)
  # Correct for timestep
  muA<-muA*dt;phi<-phi*dt;omega<-(1/omega)*dt;beta<-beta*dt;mu<-mu*dt;muL<-muL*dt
  
  #*** Build and intialise population
  U<-matrix(0,1,endp);Ls<-matrix(0,1,endp);Lr<-matrix(0,1,endp);
  As<-matrix(0,1,endp);Ar<-matrix(0,1,endp);
  Ls[1]<-initial[1];  Lr[1]<-initial[2];  As[1]<-initial[3];  Ar[1]<-initial[4];  U[1]<-N-Ls[1]-Lr[1]-As[1]-Ar[1]  
  
  #*** Fitness levels: Add columns to M and L if not enough time steps
  if(dim(M0)[2]<endp){M<-matrix(0,dim(M0)[1],endp+1); M[,1]<-M0[,1]; rownames(M)<-rownames(M0);
                      L<-matrix(0,dim(L0)[1],endp+1); L[,1]<-L0[,1]; rownames(L)<-rownames(L0);}else{M<-M0;L<-L0;}
  # How many fitness levels?
  nfit<-dim(M)[1]
  
  #*** Initial foi
  lambdasv<-matrix(0,1,endp);lambdarv<-matrix(0,1,endp); 
  lambdasv[1] = beta * As[1]/N;   lambdarv[1] = meanfit_vars(M,L,0,0,p,0,0,0,0,nfit,type,propF) * beta * Ar[1]/N; # function outputs just meanfit when all popns 0
  meanf<-c(0.6);propf<-c();
  
  #*** Main model dynamics
  for(i in 1:endp){
    lambdas=lambdasv[i];lambdar=lambdarv[i];
    # Dynamics
    U[i+1]  = U[i] + muL*(Ls[i]+Lr[i]) + muA*(As[i]+Ar[i]) - (lambdas+lambdar)*U[i]
    Ls[i+1] = Ls[i]   + (1-p)*lambdas*(U[i]+chi*Lr[i]) - (phi+muL+p*chi*lambdas +chi*lambdar)*Ls[i]  + omega*(1-ks)*As[i]
    Lr[i+1] = Lr[i]   + (1-p)*lambdar*(U[i]+chi*Ls[i]) - (phi+muL+p*chi*lambdar +chi*lambdas)*Lr[i]  + omega*(1-kr)*Ar[i]
    As[i+1] = As[i]   + p*lambdas*( U[i]+chi*(Ls[i]+Lr[i])) + phi*Ls[i]  - (omega*(1-ks)+muA)*As[i]  - As[i]*eps*omega
    Ar[i+1] = Ar[i]   + p*lambdar*( U[i]+chi*(Ls[i]+Lr[i])) + phi*Lr[i]  - (omega*(1-kr)+muA)*Ar[i]  + As[i]*eps*omega
    
    # Mean fitness update and foi
    X<-meanfit_vars(M,L,As[i]*eps*omega,p*lambdar*( U[i]+chi*(Ls[i]+Lr[i])),p,phi*Lr[i],(omega*(1-kr))*Ar[i],Ar[i]*(1-(omega*(1-kr)+muA)),Lr[i]*(1-(phi+muL+p*chi*lambdar +chi*lambdas)),nfit,type,propF)
    
    lambdasv[i+1] = beta * As[i+1]/N; 
    lambdarv[i+1] = X$meanfit * beta * Ar[i+1]/N; 
    M<-X$M
    L<-X$L
    meanf<-c(meanf,X$meanfit)
  }
  
  #*** Storing and output
  All<-c();D<-c()
  All<-as.data.frame(cbind(seq(1,endp+1,1),U,Ls,Lr,As,Ar,meanf))
  colnames(All)<-c("time","U","Ls","Lr","As","Ar","mf")
  # In form for plotting
  D1<-melt(All,id.vars="time")
  D<-cbind(D1,matrix(propF,endp+1,1))
  colnames(D)<-c(colnames(D1),"propF")
  # Values at 5 & 50 years
  if(endp > 100){level550<-c(100*(All[5*(1/dt),c("Ar")]) / ( All[5*(1/dt),"As"]+All[5*(1/dt),c("Ar")]), 100*(All[endp+1,c("Ar")]) / ( All[endp+1,"As"]+All[endp+1,c("Ar")]))} else {level550 = c(0,0)}
  # What to output (storage of all vector, plus individual levels, at 5 and 50, foi s, foi r, fitness distributions and mean relative fitness over time)
  return(list(D=D,U=U,Ls=Ls,Lr=Lr,As=As,Ar=Ar,level550=level550,lambdasv=lambdas,lambdarv=lambdarv,M=M,L=L,meanf=meanf))
}


######****************************************************************************************** Function for ode
## Function for ode calculation - "data" generation to check generalised function
twor_var<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dU = muL*(Ls+Lru+Lrf) + muA*(As+Aru+Arf) - (beta*As/N+f*beta*Aru/N+beta*Arf/N)*U
    dLs = (1-p)*beta*As/N*(U+chi*Lrf+chi*Lru) - (phi+muL+p*chi*beta*As/N +chi*beta*Arf/N+chi*f*beta*Aru/N)*Ls  + omega*(1-ks)*As
    dLru = (1-p)*f*beta*Aru/N*(U+chi*Ls+chi*Lrf) - (phi+muL+p*chi*f*beta*Aru/N+chi*beta*As/N+chi*beta*Arf/N)*Lru + omega*(1-kr)*Aru
    dLrf = (1-p)*beta*Arf/N*(U+chi*Ls+chi*Lru) - (phi+muL+p*chi*beta*Arf/N+chi*beta*As/N+chi*f*beta*Aru/N)*Lrf + omega*(1-kr)*Arf
    dAs  = p*beta*As/N*( U+chi*(Ls+Lru+Lrf))   + phi*Ls  - (omega*(1-ks)+muA)*As  - As*eps*omega
    dAru = p*f*beta*Aru/N*(U+chi*(Ls+Lru+Lrf)) + phi*Lru - (omega*(1-kr)+muA)*Aru + (1-propF)*As*eps*omega
    dArf = p*beta*Arf/N*(U+chi*(Ls+Lru+Lrf))   + phi*Lrf - (omega*(1-kr)+muA)*Arf + propF*As*eps*omega
    # Incident cases
    dIAru = p*f*beta*Aru/N*(U+chi*(Ls+Lru+Lrf)) +  (1-propF)*As*eps*omega + phi*Lru
    dIArf = p*beta*Arf/N*(U+chi*(Ls+Lru+Lrf)) +  propF*As*eps*omega + phi*Lrf
    
    trans = p*beta*(Arf/N)*(U+chi*(Ls+Lru+Lrf))+p*f*beta*(Aru/N)*(U+chi*(Ls+Lru+Lrf))
    acq = As*eps*omega
    react = phi*(Lru + Lrf)
    newU = p*f*beta*Aru/N*(U+chi*(Ls+Lru+Lrf)) + phi*Lru + 0.99*As*eps*omega
    newF = p*beta*Arf/N*(U+chi*(Ls+Lru+Lrf))   + phi*Lrf + 0.01*As*eps*omega
    oldR = (Aru+Arf)*(1-omega*(1-kr)-muA)
    # return the rate of change and all other indicators
    list(c(dU, dLs, dLru,dLrf,dAs,dAru,dArf,dIAru,dIArf),acq,trans,react,newU,newF,oldR)
  }) # end with(as.list ...
}


######****************************************************************************************** Original model
## Original difference equation of Sourya model 
sourya=function(endp,home,vary,initial){
  # endp<-3000
  # Parameters assigned by length of vary so for baseline just eps = 0
  vary_n = c("eps", "f", "krdiff")
  setwd(home) # might have para in different place for different models 
  N<-1000000; dt<-0.1 ### DON'T CHANGE TIMESTEP as fixed (3000) endp
  # Parameter assignment
  para<-read.csv("para.csv",header=TRUE,check.names=F,stringsAsFactors = FALSE); 
  for(i in 1:length(para[,1])){assign(para[i,1],para[i,2])}
  for(i in 1:length(vary)){assign(vary_n[i],vary[i])}
  assign("kr",ks + krdiff)
  muA<-muA*dt;phi<-phi*dt;omega<-(1/omega)*dt;beta<-beta*dt;mu<-mu*dt;muL<-muL*dt
  # Initialise
  U<-matrix(0,1,endp);Ls<-matrix(0,1,endp);Lr<-matrix(0,1,endp);
  As<-matrix(0,1,endp);Ar<-matrix(0,1,endp);
  Ls[1]<-initial[1];
  Lr[1]<-initial[2];
  As[1]<-initial[3];
  Ar[1]<-initial[4];
  U[1]<-N-Ls[1]-Lr[1]-As[1]-Ar[1]  
  # Initial foi
  lambdas = beta * As[1]/N; lambdar = f * beta * Ar[1]/N
  # Model 
  for(i in 1:endp){
    U[i+1]  = U[i] + muL*(Ls[i]+Lr[i]) + muA*(As[i]+Ar[i]) - (lambdas+lambdar)*U[i]
    Ls[i+1] = Ls[i] + (1-p)*lambdas*(U[i]+chi*Lr[i]) - (phi+muL+p*chi*lambdas+chi*lambdar)*Ls[i] + omega*(1-ks)*As[i]
    Lr[i+1] = Lr[i] + (1-p)*lambdar*(U[i]+chi*Ls[i]) - (phi+muL+p*chi*lambdar+chi*lambdas)*Lr[i] + omega*(1-kr)*Ar[i]
    As[i+1] = As[i] + p*lambdas*(U[i]+chi*(Ls[i]+Lr[i])) + phi*Ls[i] - (omega*(1-ks)+muA)*As[i] - As[i]*eps*omega#- As[i]*eps*(1/omega)
    Ar[i+1] = Ar[i] + p*lambdar*(U[i]+chi*(Ls[i]+Lr[i])) + phi*Lr[i] - (omega*(1-kr)+muA)*Ar[i] + As[i]*eps*omega #+ As[i]*eps*(1/omega)
    lambdas = beta * As[i]/N; lambdar = f * beta * Ar[i]/N
  }
  All<-c();D<-c()
  All<-as.data.frame(cbind(seq(1,endp+1,1),U,Ls,Lr,As,Ar))
  colnames(All)<-c("time","U","Ls","Lr","As","Ar")
  D<-melt(All,id.vars="time")
  # Value at 5 & 50 years
  level550<-c(100*(All[5*(1/dt),"Ar"]) / ( All[5*(1/dt),"As"]+All[5*(1/dt),"Ar"]), 100*(All[endp+1,"Ar"]) / ( All[endp+1,"As"]+All[endp+1,"Ar"]))
  return(list(D=D,U=U,Ls=Ls,Lr=Lr,As=As,Ar=Ar,level550=level550))
}


######****************************************************************************************** For multiple plots
## Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}