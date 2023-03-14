## ----label=ddataprep----------------------------------------------------------
set.seed(125)
# retrieve the data
library(MASS)
data(birthwt)
birthwt$low=NULL 
# remove column for the low indicator 
# which is 1 i a child has low birtwt
rnd=round(rnorm(nrow(birthwt),mean=10,sd=2),2)
# rnd just contains random data
birthwt=cbind(birthwt,rnd=rnd) # adding it
head(birthwt)

## ----label=plotsnha, fig.width=10,fig.height=5,fig.cap="Birth weight data variable interactions"----
opar=par(mfrow=c(1,2),mai=c(0.8,0.8,0.1,0.2))
library(snha)
# retrieve some data
pca=prcomp(t(scale(birthwt)))
summary(pca)
plot(pca$x[,1:2],xlab='PC1', ylab='PC2',pch=19,cex=5,col='salmon')
text(pca$x[,1:2],colnames(birthwt))
text(-5,-10,"PCA",cex=2) 
as=snha(birthwt,method="spearman",alpha=0.1)
par(mai=c(0.8,0.2,0.1,0.2))
plot(as,layout="sam",vertex.size=7,lwd=3,edge.width=3)
text(-1.5,-1.8,"SNHA",cex=2)
box()
par(opar)

## -----------------------------------------------------------------------------
round(snha_rsquare(as),2)
as$theta

## ----label=dec,fig.width=9,fig.height=6,fig.cap="SNHA - Decathlon Data 1988"----
### data loading
data(decathlon88)
head(decathlon88)
A=snha(decathlon88,method="spearman",alpha=0.1)
cols=rep("salmon",10)
cols[names(A$data) %in% c("jave","shot","disc","pole")]="skyblue"
plot(A,layout="sam",vertex.color=cols,vertex.size=8,cex=1.1,edge.width=5)
snha_rsquare(A)
mn=mean(snha_rsquare(A))
title(paste("R-square = ",round(mn,2)))

## -----------------------------------------------------------------------------
round(A$sigma,2)
round(A$p.value,3)

## ----label=dec2,fig.width=9,fig.height=6,fig.cap="Decathlon Data 1988 (p-value Graph)"----
B = A$theta
B[]=0
B[A$p.value<0.05]=1
diag(B)=0
plot.snha(B,layout='sam',vertex.color=cols,vertex.size=8,cex=1.1,edge.width=5)

## ----label=plot,fig.width=10,fig.height=5,fig.cap="Swiss data variable associations"----
library(snha)
data(swiss)
head(swiss,4)
### shorter names useful for display later in the graph
colnames(swiss)=abbreviate(colnames(swiss))
head(swiss,4)
opar=par(mfrow=c(1,2))
### options(warn=-1)
as=snha(swiss,method="pearson")
### store layout for reuse in two graphs
lay = snha_layout(as,mode="sam")
plot(as,layout=lay,vertex.size=8,main="Pearson")
as=snha(swiss,method="spearman")
plot(as,layout=lay,vertex.size=8,main="Spearman")
par(opar)

## ----label=theta,results='asis'-----------------------------------------------
knitr::kable(as$theta)

## ---- result="asis"-----------------------------------------------------------
### prepare a test returning only p-values
mtest = function (x) { return(shapiro.test(x)$p.value)  }
df=data.frame(orig=round(apply(swiss,2,mtest),3))
df=cbind(df,log2=round(apply(log2(swiss),2,mtest),3))
knitr::kable(df)

## ----label=corrplot,fig.width=10,fig.height=5,fig.cap="Correlation and Network plot with correlation values on the edges"----
opar=par(mfrow=c(1,2),mai=rep(0.2,4))
sw=snha(swiss,method="spearman",alpha=0.1)
plot(sw,type="corrplot")
plot(as,edge.text=round(as$sigma,2),edge.pch=15,layout='sam')
par(opar)

## ----label=chains-------------------------------------------------------------
snha_get_chains(as)

## ----label=ll-----------------------------------------------------------------
snha_ll(as)

## ----label=boot,fig.width=9,fig.height=4,fig.cap="Boostrap Example"-----------
opar=par(mfrow=c(1,2),mai=c(0.1,0.1,0.7,0.1))
as.boot=snha(swiss,method="spearman",prob=TRUE)
lay=snha_layout(as.boot,method="sam")
plot(as,layout=lay,vertex.size=6,main="Single Run")
plot(as.boot,layout=lay,vertex.size=6,main="Bootstrap Run")
par(opar)

## ----werner-------------------------------------------------------------------
W=matrix(0,nrow=6,ncol=6,dimnames=list(LETTERS[1:6],LETTERS[1:6]))
W[1:2,3]=1
W[3,4]=1
W[4,5:6]=1
W[5,6]=1
W

## -----------------------------------------------------------------------------
data=snha_graph2data(W)
dim(data)
round(cor(t(data)),2)

## ----label=wplot,fig.width=8,fig.height=3,out.width=900,fig.cap="True graph, correlations and predicted graph (left to right)"----
opar=par(mfrow=c(1,3),mai=rep(0.2,4))
plot.snha(W)
plot.snha(cor(t(data)),type="cor")
plot.snha(snha(t(data)))
par(opar)

## ----label=start,fig.width=6,fig.height=1.7,fig.cap="An association chain"----
opar=par(mai=c(0.1,0.1,0.1,0.0))
plot(1,xlab="",ylab="",axes=FALSE,type="n",xlim=c(0.5,4.5),ylim=c(0.8,1.2))
arrows(1:3,rep(1,3),1:3+0.8,rep(1,3),lwd=3,length=0.1)
points(1:4,rep(1,4),pch=19,col="salmon",cex=6)
text(1:4,1,LETTERS[1:4],cex=2)
par(opar)

## ----result='asis'------------------------------------------------------------
C=matrix(c(1,0.7,0.5,0.3,
           0.7,1,0.7,0.5,
           0.5,0.7,1,0.7,
           0.3,0.5,0.7,1),
           nrow=4,byrow=TRUE)
rownames(C)=colnames(C)=LETTERS[1:4]           
knitr::kable(C)

## ----label=corplot,fig.width=6,fig.height=3,fig.cap="Visualization of correlation matrix sigma and the adjacency matrix theta"----
set.seed(123)
opar=par(mfrow=c(1,2),mai=c(0.1,0.1,0.1,0.1))
C=C+rnorm(length(C),mean=0,sd=0.1)
C[lower.tri(C)]=t(C)[lower.tri(C)]
diag(C)=1
as=snha(C)
round(as$sigma,3)
plot(as,type="corplot")
as$theta
plot(as)
par(opar)

## -----------------------------------------------------------------------------
print(sessionInfo())

