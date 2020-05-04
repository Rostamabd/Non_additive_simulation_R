rm(list=ls())
## we set the seed here to tell the program that always start from the same random number. 
set.seed(12345)
#==========# Loading files ####################==================================
library(BGLR)

### replace your genotype data with mice genotype data

data(mice)

X <- mice.X[1:1000,1:5000]

### here no. of total QTL is 100, you can change them to match with your simulation

nqtl <-100 ; m <-ncol(X); n<-nrow(X) ; h2<-0.30

### randomly selected nqtl SNPs as QTL
whichQTL<- sort(sample(1:m,nqtl))

### create Genotype matrix
qtlmatrix <-X[,whichQTL]
dim(qtlmatrix)

##########  Calculate the allelic frequency ####
p <- (colMeans(qtlmatrix))/2
q <- 1-p
 
K<-sum(apply(FUN=var,X=X[,whichQTL],MARGIN=2))

#### generating the additive QTL effects from a normal distribution##=============

b0 <-rep(0,m)

### additive effects: approximately 50% negative and 50% positive

a <- rnorm(n=nqtl, mean=0,sd=sqrt(h2/K))

### dominance is set to zero here since only additive architecture is considered here

d=0

b0[whichQTL]<- a+d*(q-p)

signal_A<-X%*%b0

### achieve desirable h2
vg <- var(signal_A)
vp <- vg/h2
var.env <- vp-vg; error<-rnorm(n,sd=sqrt(var.env))

###===================== Simulate additive phenotypes ###========================
y <- signal_A+ error

var(signal_A)/var(y)

###========== generating the domonance effects##==================================

###scaled measure of the dominance

k <- rnorm(n=nqtl, 0.5,1)

### creat the design matrix
D <- X
D[D==2]<- 0

d0 <- rep(0,m)
d=a*abs(k)

d0[whichQTL] <- d


signal_D<-D%*%d0

### achieve desirable h2
h2= h2+0.10
vg <- var(signal_A+ signal_D)
vp <- vg/h2
var.env <- vp-vg; error <- rnorm(n,sd=sqrt(var.env))

###=============# generate additive+dominance phenotypes ###=====================
y <- signal_A+ signal_D+error

var(signal_A)/var(y)

var(signal_D)/var(y)

###====================## Epistatic effects: each SNP have intreaction with three adjucent SNPs
v1 <- cbind(rep(1:99,1),rep(2:100,1))
v2 <- cbind(rep(1:98,1),rep(3:100,1))
v3 <- cbind(rep(1:97,1),rep(4:100,1))

epis <- rbind(v1,v2,v3)
epis <- epis[order(epis[,1]),]
n.epis <- dim(epis)[1]
### define the epistatic effects

### aditive by additive
add_add <-rgamma(n.epis,0.1,10)
signaa <-rep(c(-1,1),n.epis/2)
add_add <- signaa*add_add

### addititive by dominance
add_dom <-rgamma(n.epis,0.1,10)
signaa <-rep(c(-1,1),n.epis/2)

add_dom <- signaa*add_dom

### dominance by additive
dom_add <-rgamma(n.epis,0.1,10)
signaa <-rep(c(-1,1),n.epis/2)

dom_add <- signaa*dom_add

### Dominance by dominance
dom_dom <-rgamma(n.epis,0.1,10)
signaa <-rep(c(-1,1),n.epis/2)

dom_dom  <- signaa*dom_dom

#### Creat the final matrix for the epistatic effects
qtl.epis <- cbind(add_add,add_dom,dom_add,dom_dom)

qtl.epis <- cbind(epis,qtl.epis)			   


### new matrices for constracting add*add and add*dom epistatics effects
 
 xa1 <- qtlmatrix[ ,qtl.epis[,1]]
 xd1 <- D[ ,qtl.epis[,1]]
 xa2 <- qtlmatrix[ ,qtl.epis[,2]]
 xd2 <- D[ ,qtl.epis[,2]]

###================ generate additive+dominance+epistatistic architecture ##==================

signal_Ipis<- (xa1*xa2)%*%qtl.epis[,3]+(xa1*xd2)%*%qtl.epis[,4]+(xd1*xa2)%*%qtl.epis[,5]+(xd1*xd2)%*%qtl.epis[,6]

#######  achieve desirable h2 
h2=h2+0.30
vg <- var(signal_A+ signal_D+signal_Ipis)
vp <- vg/h2
var.env <- vp-vg; error<-rnorm(n,sd=sqrt(var.env))


y <- signal_A+ signal_D+signal_Ipis+error

########################################
### Genomic evaluation by BayesB     ###
########################################


#### Evaluation using all markers (Qausal markers ignored)
#X <- X[,-whichQTL]

### Evaluation only using cusal variants
X <- X[,whichQTL]

#####================ #run the statistical model #================================

#### testing set

### this evaluation is based on 10 replications and 5-fold cross-evalution.
reps <- 10
nFolds <- 5
ntest <- n/nFolds


### some jars to save the results
mse_A <- numeric()
predcor_Ap <- numeric()
predcor_As <- numeric()

### beggining of evaluations

for (i in 1:reps){
cat("Replication:", i,"\n")

### Splitting the data into differents folds

tst <-sample(1:n,size=ntest,replace=FALSE)

yNA <- y

### modify y to yNA
yNA[tst] <- NA

#### Additive model

ETA<-list(list(X=X, model='BayesB'))

fmA<-BGLR(y=yNA,ETA=ETA,thin=5, nIter=20000,burnIn=5000,verbose = FALSE,saveAt="Pred_")


mse_A[i] <- mean((signal_A[tst]-fmA$yHat[tst])^2)
predcor_Ap [i] <- cor(signal_A[tst],fmA$yHat[tst])
predcor_As[i] <- cor(signal_A[tst],fmA$yHat[tst],method ="spearman")

}

Results <-  cbind(predcor_Ap,predcor_As, mse_A) 

### Total criteria for prediction
trainfile <- paste(i,"output.txt",sep="")

write.table(Results,file=trainfile, quote = FALSE, col.names = c("Pred_Corr_Pearson", "Pred_Corr_Spearman","MSE"))



