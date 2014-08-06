### Distributions


##*** Libraries needed
library(plyr); library(ggplot2);library(reshape2);library(deSolve);library(grid)

##*** Locations
home<-"~/Documents/My Documents/Diagnostics/Natural history/genfun/GenFun/data"
plots<-paste(home,"/plots",sep="")
setwd(home)

##*** Gagneaux 2006
g<-read.csv("dist_gagn2006.csv"); d<-c()
for (i in 1:dim(g)[1]){d<-c(d,matrix(1,1,g[i,2])*g[i,1])}

setwd(plots)
pdf("data_gagn2006.pdf")
h <- hist(d, breaks = seq(0,1,0.01),main="rpoB mutants, CDC1551 (Gagneux, 2006)",xlab = "Relative fitness (competitive growth)")
dev.off()

setwd(home)
g1<-read.csv("dist_gagn2006t85.csv")[1:4,]; d1<-c()
for (i in 1:dim(g1)[1]){d1<-c(d1,matrix(1,1,g1[i,2])*g1[i,1])}

setwd(plots)
pdf("data_gagn2006t85.pdf")
h1 <- hist(d1, breaks = seq(0,1,0.01),ylim=c(0,25),main="rpoB mutants, T85 (Beijing) \n (Gagneux, 2006)",xlab = "Relative fitness (competitive growth)")
dev.off()

pdf("gagn2006_diff.pdf")
plot( h, col=rgb(0,0,1,1/4), xlim=c(0,1),ylim=c(0,25),main="rpoB mutations (red = Beijing) \n (Gagneux 2006)",xlab = "Relative fitness (competitive growth)")  # first histogram
plot( h1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T)  # second
dev.off()

## Pool
gg<-cbind(c(g[,1],g1[,1]),c(g[,2],g1[,2])); d2<-c()
for (i in 1:dim(gg)[1]){d2<-c(d2,matrix(1,1,gg[i,2])*gg[i,1])}
pdf("gagn2006_pool.pdf")
h2 <- hist(d2, breaks = seq(0,1,0.01),ylim=c(0,25),main="rpoB mutants \n (Gagneux, 2006)",xlab = "Relative fitness (competitive growth)")
dev.off()


##************************************************************************ Mariam 2004
setwd(home)
g<-read.csv("dist_mariam2004.csv"); d<-c();d1<-c(); g[,3]<-g[,3]
for (i in 1:dim(g)[1]){d<-c(d,matrix(1,1,g[i,3])*g[i,2]);d1<-c(d1,matrix(g[i,1],1,g[i,3]))}
bb<-as.data.frame(d); bb$exp<-d1; colnames(bb)<-c("count","Experiment")
p<-ggplot(bb, aes(count, fill = Experiment)) + geom_histogram(alpha = 0.8, aes(y = ..density..), position = 'identity',binwidth=0.01)+ geom_density(alpha=0.5)
p<-p+scale_x_continuous(lim=c(0,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("rpoB mutants (Harlingen) \n (Mariam 2004)")
p

setwd(plots)
ggsave("mariam_pool.pdf")


##************************************************************************ Pool all rpoB data
setwd(home)
g<-read.csv("rpob_all.csv")[1:22,]; d1<-c();d2<-c();d3<-c();d4<-c();d<-c()
for (i in 1:dim(g)[1]){d<-c(d,matrix(1,1,g[i,3])*g[i,2]);d1<-c(d1,matrix(g[i,1],1,g[i,3]))
                       d2<-c(d2,matrix(g[i,4],1,g[i,3]))
                       d3<-c(d3,matrix(g[i,5],1,g[i,3]))
                       d4<-c(d4,matrix(g[i,6],1,g[i,3]))}
bb<-as.data.frame(d); bb$exp<-d1;bb$paper<-d2;bb$rpoB<-d3;bb$strain<-d4; colnames(bb)<-c("count","Experiment","Ref","rpoBmut","Strain")
p1<-ggplot(bb,aes(count,fill=Experiment))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)+geom_density(alpha=0.5)
p1<-p1+scale_x_continuous(lim=c(0,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All rpoB mutants")
setwd(plots)
ggsave("all_rpoB_exp.pdf")
w<-which(bb[,"Experiment"]=="Independent")
p2<-ggplot(bb[w,],aes(count,fill=Ref))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)+geom_density(alpha=0.5)
p2<-p2+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All rpoB mutants, \n independent growth experiments")
p2
ggsave("all_rpoB_ref.pdf")

w<-which(bb[,"Experiment"]=="Independent")
p3<-ggplot(bb[w,],aes(count,fill=rpoBmut))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)
p3<-p3+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All rpoB mutants, \n independent growth experiments")
p3
ggsave("all_rpoB_mut.pdf")

##************************************************************************ Sander 2002 M. smeg.
setwd(home)
g<-read.csv("dist_sander2002.csv"); d<-c()
for (i in 1:dim(g)[1]){d<-c(d,matrix(1,1,g[i,2])*g[i,1])}

setwd(plots)
pdf("data_sander2002.pdf")
h <- hist(d, breaks = seq(0,1,0.01),main="Ribosomal nucleic acids mutation, \n M. smegmatis (Sander, 2002)",xlab = "Relative fitness (competitive growth)")
dev.off()

##************************************************************************ Pool all smeg data
setwd(home)
g<-read.csv("all_smeg.csv")[1:16,]; d1<-c();d2<-c();d3<-c();d4<-c();d<-c()
for (i in 1:dim(g)[1]){d<-c(d,matrix(1,1,g[i,3])*g[i,2]);d1<-c(d1,matrix(g[i,1],1,g[i,3]))
                       d2<-c(d2,matrix(g[i,4],1,g[i,3]))
                       d3<-c(d3,matrix(g[i,5],1,g[i,3]))
                       d4<-c(d4,matrix(g[i,6],1,g[i,3]))}
bb<-as.data.frame(d); bb$exp<-d1;bb$paper<-d2;bb$res<-d3;bb$mutation<-d4; colnames(bb)<-c("count","Experiment","Ref","Resistance","Mutation")
p1<-ggplot(bb,aes(count,fill=Experiment))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)+geom_density(alpha=0.5)
p1<-p1+scale_x_continuous(lim=c(0,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All smegmatis data")
setwd(plots)
ggsave("all_smeg_exp.pdf")
w<-which(bb[,"Experiment"]=="Independent")
p2<-ggplot(bb[w,],aes(count,fill=Ref))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)
p2<-p2+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All smegmatis data, \n independent growth experiments")
p2
ggsave("all_smeg_ref.pdf")

w<-which(bb[,"Experiment"]=="Independent")
p3<-ggplot(bb[w,],aes(count,fill=Mutation))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)
p3<-p3+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All smegmatis data, \n independent growth experiments")
p3
ggsave("all_smeg_mut.pdf")

w<-which(bb[,"Experiment"]=="Independent")
p3<-ggplot(bb[w,],aes(count,fill=Resistance))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)
p3<-p3+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All smegmatis data, \n independent growth experiments")
p3
ggsave("all_smeg_res_ind.pdf")

w<-which(bb[,"Experiment"]=="Competition")
p3<-ggplot(bb[w,],aes(count,fill=Resistance))+geom_histogram(alpha=0.8,position = 'identity',binwidth=0.01)
p3<-p3+scale_x_continuous(lim=c(0.25,1),"Relative fitness")+scale_y_continuous("Frequency")+ggtitle("All smegmatis data, \n competition experiments")
p3
ggsave("all_smeg_res_comp.pdf")
