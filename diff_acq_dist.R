Sv<-sourya_funcf_mean_vars(100*(1/dt),home, c(0.008,0.6,0.35),iniv,M0,L0,"orig",try[i])
plot(Sv$meanf,ylim=c(0,1),type="l",col="blue")
points(seq(0,100*1/0.1,1/0.1),out$meanf,col="black")
legend(200,0.6,c("2 levels ode model","Generalised fit function"),c("blue","black"))
plot(Sv$"Ar",type="l",col="blue")
points(seq(0,100*1/0.1,1/0.1),out$activeR,col="black")
legend(100,150,c("2 levels ode model","Generalised fit function"),c("blue","black"))

mm<-as.data.frame(cbind(seq(1,dim(Sv$M)[2],1),t(Sv$M)))
colnames(mm)<-c("t",seq(1/n,1,1/n))

mmt<-melt(mm,id.vars="t")
ggplot(mmt,aes(x=t,y=value,colour=variable))+geom_line()
