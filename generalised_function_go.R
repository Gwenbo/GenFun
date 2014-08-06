##### Code to run generalised function and examples

##*** Libraries needed
library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid)

##*** Locations
home<-"~/Documents/My Documents/Diagnostics/Natural history/genfun/"
plots<-paste(home,"plots",sep="")
setwd(home)

##*** Code needed
# Loads functions for generalised mean function and simulation model that uses the generalised function
# Also 2 fitness level ode function, original Sourya model as a difference model and multiplots function
source("generalised_function.R") 

##*** Setting up
# Number of discrete levels? 
nfit = 10
# Matrices of distribution of fitness
M0 <- matrix(0,nfit,10)
rownames(M0)<-seq(1/nfit,1,1/nfit) # Currently relative fitness of 1 is upper limit
L0 <- matrix(0,nfit,10)
rownames(L0)<-seq(1/nfit,1,1/nfit)

#*** Baseline model - for "fitting" initial conditions 
# Timestep (dt = 0.1) v important! 
S<-sourya(3000,home, c(0,0.6,0.035),c(0,0,100,0))
iniv<-c(S$Ls[3001],S$Lr[3001],S$As[3001],S$Ar[3001],0,0)
Ls0=S$Ls[3001]; As0=S$As[3001]; N=1000000

# Timestep
dt = 0.01

#*** For "data" generation - ode model with 2 levels 
parameters_var <- c(muL=0.02,muA=0.187,beta=7.36,phi=0.0015,p=0.14,omega=1,ks=0.05,kr=0.4,chi=0.33,eps=0.008,mu=0.02,f=0.6,propF=0.01)
times <- seq(0, 100, by = 1)
state <- c(U = N-Ls0-As0,Ls=Ls0,Lru=0,Lrf=0,As=As0,Aru=0,Arf=0,IAru=0,IArf=0) # Takes initial conditions from above baseline model
out <- data.frame(ode(y = state, times = times, func = twor_var, parms = parameters_var))
# additional output 
out$meanf<-(out[, "Arf"]/(out[, "Aru"]+out[, "Arf"])) + 0.6*(out[, "Aru"]/(out[, "Aru"]+out[, "Arf"])); out$meanf[1]<-0.6
out$propF<-0.01; out$activeR<-out[,"Arf"]+out[,"Aru"]

#############********************************************** Testing
#*** Mean fitness function is, with mock numbers
X<-meanfit_vars(M0,L0,10,0,0.14,1,0,0,0,nfit,"orig",0.01)
# Produces output for next time step's mean fitness (to generate new number of transmissions) and updates next "column" or distribution of active (M) and latent (L)
X$meanfit
X$M
X$L

#*** Run model with the above function via:
# Here inputs are timesteps, location (for parameters), varying para (the 3 key drivers), initial conditions (for pop and distributions), type of acquisition distribution and proportion acq are "fit"
Sv<-sourya_funcf_mean_vars(100*(1/dt),home, c(0.008,0.6,0.35),iniv,M0,L0,"orig",0.1)
# Quick check plot of above output (mean fitnes over time and active resistant, blue) and output from ode model with 2 hetero levels (black)
plot(Sv$meanf,ylim=c(0.5,1),type="l",col="blue")
points(seq(0,100*1/dt,1/dt),out$meanf,col="black") # Have to plot out against 10*time due to timestep of dt
legend(200,0.6,c("2 levels ode model","Generalised fit function"),c("black","blue"))
plot(Sv$"Ar",type="l",col="blue",ylim=c(0,600))
points(seq(0,100*1/dt,1/dt),out$activeR,col="black")
legend(100,150,c("2 levels ode model","Generalised fit function"),c("black","blue"))


#*** Multiple proportions of acquisitions moving to high level fitness
# Ugly plotting but just to show for now
# Store data frames
G<-data.frame(matrix(0,1,4)); colnames(G)<-c(colnames(Sv$D)) # Store generalised function output
H<-data.frame(matrix(0,1,19)); colnames(H)<-colnames(out) # Store "data" from ode model
# values for propF to try
try<-seq(0.01,0.1,0.01)
for (i in 1:length(try)){
  Smm<-sourya_funcf_mean_vars(100*(1/dt),home, c(0.008,0.6,0.35),iniv,M0,L0,"orig",try[i])
  G<-rbind(G,Smm$D)  
  parameters_var <- c(muL=0.02,muA=0.187,beta=7.36,phi=0.0015,p=0.14,omega=1,ks=0.05,kr=0.4,chi=0.33,eps=0.008,mu=0.02,f=0.6,propF=try[i])
  out <- data.frame(ode(y = state, times = times, func = twor_var, parms = parameters_var))
  out$meanf<-(out[, "Arf"]/(out[, "Aru"]+out[, "Arf"])) + 0.6*(out[, "Aru"]/(out[, "Aru"]+out[, "Arf"])); out$meanf[1]<-0.6; out$propF<-try[i];  out$activeR<-out[,"Arf"]+out[,"Aru"]
  H<-rbind(H,out)
}
setwd(plots)
# Plot time series
H[,"time"]=H[,"time"]*(1/dt) # as difference between ode solved time and time steps (0.1) of 

p1=ggplot(G[which(G[,"propF"]==0.01&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p1=p1+geom_line(data=H[which(H[,"propF"]==0.01),],aes(x=time,y=activeR))+ggtitle("propF=1%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2=ggplot(G[which(G[,"propF"]==0.02&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p2=p2+geom_line(data=H[which(H[,"propF"]==0.02),],aes(x=time,y=activeR))+ggtitle("propF=2%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3=ggplot(G[which(G[,"propF"]==0.03&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p3=p3+geom_line(data=H[which(H[,"propF"]==0.03),],aes(x=time,y=activeR))+ggtitle("propF=3%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4=ggplot(G[which(G[,"propF"]==0.04&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p4=p4+geom_line(data=H[which(H[,"propF"]==0.04),],aes(x=time,y=activeR))+ggtitle("propF=4%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5=ggplot(G[which(G[,"propF"]==0.05&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p5=p5+geom_line(data=H[which(H[,"propF"]==0.05),],aes(x=time,y=activeR))+ggtitle("propF=5%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p6=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[7]&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p6=p6+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[7]),],aes(x=time,y=activeR))+ggtitle("propF=6%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p7=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[8]&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p7=p7+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[8]),],aes(x=time,y=activeR))+ggtitle("propF=7%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p8=ggplot(G[which(G[,"propF"]==0.08&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p8=p8+geom_line(data=H[which(H[,"propF"]==0.08),],aes(x=time,y=activeR))+ggtitle("propF=8%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p9=ggplot(G[which(G[,"propF"]==0.09&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p9=p9+geom_line(data=H[which(H[,"propF"]==0.09),],aes(x=time,y=activeR))+ggtitle("propF=9%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p10=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[11]&G[,"variable"]=="Ar"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Number with Ar")
p10=p10+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[11]),],aes(x=time,y=activeR))+ggtitle("propF=10%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("overtimefits_Ar (dt=0.01).pdf")
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,cols=5)
dev.off()

p1=ggplot(G[which(G[,"propF"]==0.01&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p1=p1+geom_line(data=H[which(H[,"propF"]==0.01),],aes(x=time,y=meanf))+ggtitle("propF=1%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2=ggplot(G[which(G[,"propF"]==0.02&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p2=p2+geom_line(data=H[which(H[,"propF"]==0.02),],aes(x=time,y=meanf))+ggtitle("propF=2%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3=ggplot(G[which(G[,"propF"]==0.03&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p3=p3+geom_line(data=H[which(H[,"propF"]==0.03),],aes(x=time,y=meanf))+ggtitle("propF=3%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4=ggplot(G[which(G[,"propF"]==0.04&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p4=p4+geom_line(data=H[which(H[,"propF"]==0.04),],aes(x=time,y=meanf))+ggtitle("propF=4%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5=ggplot(G[which(G[,"propF"]==0.05&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p5=p5+geom_line(data=H[which(H[,"propF"]==0.05),],aes(x=time,y=meanf))+ggtitle("propF=5%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p6=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[7]&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p6=p6+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[7]),],aes(x=time,y=meanf))+ggtitle("propF=6%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p7=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[8]&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p7=p7+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[8]),],aes(x=time,y=meanf))+ggtitle("propF=7%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p8=ggplot(G[which(G[,"propF"]==0.08&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p8=p8+geom_line(data=H[which(H[,"propF"]==0.08),],aes(x=time,y=meanf))+ggtitle("propF=8%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p9=ggplot(G[which(G[,"propF"]==0.09&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p9=p9+geom_line(data=H[which(H[,"propF"]==0.09),],aes(x=time,y=meanf))+ggtitle("propF=9%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p10=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[11]&G[,"variable"]=="mf"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p10=p10+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[11]),],aes(x=time,y=meanf))+ggtitle("propF=10%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("overtimefits_mf (dt=0.01).pdf")
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,cols=5)
dev.off()
# Plot latent resistant population over time 
H[,"latentR"]=H[,"Lru"]+H[,"Lrf"]
p1=ggplot(G[which(G[,"propF"]==0.01&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p1=p1+geom_line(data=H[which(H[,"propF"]==0.01),],aes(x=time,y=latentR))+ggtitle("propF=1%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2=ggplot(G[which(G[,"propF"]==0.02&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p2=p2+geom_line(data=H[which(H[,"propF"]==0.02),],aes(x=time,y=latentR))+ggtitle("propF=2%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3=ggplot(G[which(G[,"propF"]==0.03&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p3=p3+geom_line(data=H[which(H[,"propF"]==0.03),],aes(x=time,y=latentR))+ggtitle("propF=3%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4=ggplot(G[which(G[,"propF"]==0.04&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p4=p4+geom_line(data=H[which(H[,"propF"]==0.04),],aes(x=time,y=latentR))+ggtitle("propF=4%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5=ggplot(G[which(G[,"propF"]==0.05&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p5=p5+geom_line(data=H[which(H[,"propF"]==0.05),],aes(x=time,y=latentR))+ggtitle("propF=5%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p6=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[7]&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p6=p6+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[7]),],aes(x=time,y=latentR))+ggtitle("propF=6%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p7=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[8]&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p7=p7+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[8]),],aes(x=time,y=latentR))+ggtitle("propF=7%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p8=ggplot(G[which(G[,"propF"]==0.08&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p8=p8+geom_line(data=H[which(H[,"propF"]==0.08),],aes(x=time,y=latentR))+ggtitle("propF=8%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p9=ggplot(G[which(G[,"propF"]==0.09&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p9=p9+geom_line(data=H[which(H[,"propF"]==0.09),],aes(x=time,y=latentR))+ggtitle("propF=9%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p10=ggplot(G[which(G[,"propF"]==unique(G[,"propF"])[11]&G[,"variable"]=="Lr"),],aes(x=time,y=value))+geom_line(col="red")+scale_y_continuous("Mean fitness")
p10=p10+geom_line(data=H[which(H[,"propF"]==unique(G[,"propF"])[11]),],aes(x=time,y=latentR))+ggtitle("propF=10%")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

setwd(plots)
pdf("overtimefits_Lr (dt=0.01).pdf")
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,cols=5)
dev.off()
#*** Nicer combination fits
# Plot fit to active resistant over time and propF
Gty<-G[which(G[,"variable"]=="Ar"),]
ggplot(Gty,aes(x=time,y=value,colour=factor(propF)))+geom_line()+geom_point(data=H,aes(x=time,y=activeR,colour=factor(propF)))+scale_color_discrete("propF")+scale_y_continuous("Number with active cases with resistant strain")
ggsave("fit_Ar (dt=0.01).pdf")
# Plot fit to mean relative fitness over time and propF
Gty<-G[which(G[,"variable"]=="mf"),]
ggplot(Gty,aes(x=time,y=value,colour=factor(propF)))+geom_line()+geom_point(data=H,aes(x=time,y=meanf,colour=factor(propF)))+scale_color_discrete("propF")+scale_y_continuous(lim=c(0.5,1),"Mean relative fitness of resistant strains")
ggsave("fit_mf (dt=0.01).pdf")
