setwd("~/OneDrive - The Ohio State University/[2019-2023 OSU]/2-2021-Spring/STAT7730 Advanced Comoputational Statistics (Dr. Chang)/Final project")

#############
## Example ##
#############

# q : number of nodes (including the response Y)
# n : number of observations  


## Generate a DAG
library(pcalg)
library(extraDistr)
library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(gridExtra)

set.seed(123)
source("obayes_inference_causal_effects/OBES_mcmc.r")
source("obayes_inference_causal_effects/bayes_causal.r")

load("R codes and results/datan100q10p01distinct.RData")

isugMAT_ <- function(A_) {
  .Call('_gRbase_isugMAT_', PACKAGE = 'gRbase', A_)
}

shd.new <- function(m1, m2){
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}


q = 5
n = 50
M    = 20
T    = 600
burn = 100
S = 50
p.star = 0.1
graphnum<-40

OBMA_Int_eff.list<-rep(0,graphnum)
OBMED_Int_eff.list<-rep(0,graphnum)
Int_eff.list<-rep(0,graphnum)
Int_eff_true.list<-rep(0,graphnum)
SHD.list<-rep(0,graphnum)

for(d in 1:graphnum){
  
  true.dag = randomDAG(n = q, prob = p.star)
  # plot(true.dag)
  true.dag.mat = t(as(true.dag, "matrix"))
  
  true.beta = true.dag.mat*matrix(runif(q*(q-1), -2, -1), q, q)
  diag(true.beta) = 1
  
  sigma.cond = diag(rep(1, q))
  Sigma = t(solve(true.beta))%*%sigma.cond%*%solve(true.beta)
  
  Omega = solve(Sigma)
  
  mu_0 = c(rep(0, q))
  X = rmvnorm(n, mu_0, Sigma)
  
  m = colMeans(X)
  s = apply(X = X, FUN = sd, MARGIN = 2)
  
  X = t((t(X) - m)/s)
  
  out = OBES_mcmc_chain(Y = X, M, T, verbose = FALSE)
  MCMC_chain = out[[1]]
  
  post_probs = matrix((rowMeans(MCMC_chain[,(burn + 1):(T)])),q,q)
  
  rownames(post_probs) = colnames(post_probs) = 1:q
  median_graph = (post_probs > 0.5)*1
  
  #proj_median_graph = dag2essgraph(median_graph)
  proj_median_graph = as(median_graph, "matrix")*1
  
  x.pos = rdunif(1,2,q)
  
  OBMA_output = OBMA_function(X, MCMC_chain, M, x.pos, y.pos = 1, S)
  
  OBMA_Int_eff.list[d]<-OBMA_output$OBMA_average
  
  EG_estimate = proj_median_graph
  
  OBMED_output = OBMED_function(X, EG_estimate, x.pos, y.pos = 1, S = S)
  
  OBMED_Int_eff.list[d]<-OBMED_output$OBMED_average
  
  fa_i = fa(x.pos, true.dag)
  Int_eff_true.list[d] = (Sigma[1, fa_i]%*%solve(Sigma[fa_i, fa_i]))[1,1]
  
  Distinc_graph<-distinct_causal(x.pos, y.pos = 1,  true.dag)
  Gammas = lapply(FUN = posterior_gamma, X = Distinc_graph, Y = X, x.pos = x.pos, y.pos = 1, S = S) 
  average_Gammas = lapply(FUN = mean, X = Gammas)
  Int_eff.list[d] = mean(unlist(average_Gammas))
  
  SHD.list[d]<-shd.new(true.dag.mat, EG_estimate)
  
  print(d)
  
}

dist_OBMA<-abs(OBMA_Int_eff.list-Int_eff_true.list)
dist_OBMED<-abs(OBMED_Int_eff.list-Int_eff_true.list)

save.image("/Users/xinyuxuan/Desktop/7730 projet papers/obayes_inference_causal_effects/obayes_inference_causal_effects/datan50q10p01distinct.RData")




#### Try implementing other methods


PC.graph <- function(niter = 40, a=0.1, q=10, n=50){
  PC <- NULL; PC.graph <- list()
  set.seed(123)
  for(i in 1:niter){
    true.dag = randomDAG(n = q, prob = p.star)
    # plot(true.dag)
    true.dag.mat = t(as(true.dag, "matrix"))
    
    true.beta = true.dag.mat*matrix(runif(q*(q-1), -2, -1), q, q)
    diag(true.beta) = 1
    
    sigma.cond = diag(rep(1, q))
    Sigma = t(solve(true.beta))%*%sigma.cond%*%solve(true.beta)
    
    Omega = solve(Sigma)
    
    mu_0 = c(rep(0, q))
    X = rmvnorm(n, mu_0, Sigma)
    
    m = colMeans(X)
    s = apply(X = X, FUN = sd, MARGIN = 2)
    
    X = t((t(X) - m)/s)
    
    n <- nrow(X)
    V <- as.character(1:ncol(X))
    
    res <- pcalg::skeleton(list(C = cor(X), n = n), alpha=a, labels = V,
                           indepTest = gaussCItest)
    PC[i] <- shd(true.dag, res@graph)
    PC.graph[[i]] <- res@graph
  }
  return(list(PC.shd=PC, graph=PC.graph))
}

PC.0.1.n50 <- PC.graph(a = 0.1, n=50)
PC.0.1.n100 <- PC.graph(a = 0.1, n=100)
PC.0.1.n200 <- PC.graph(a = 0.1, n=200)

PC.0.05.n50 <- PC.graph(a = 0.05, n=50)
PC.0.05.n100 <- PC.graph(a = 0.05, n=100)
PC.0.05.n200 <- PC.graph(a = 0.05, n=200)

PC.0.01.n50 <- PC.graph(a = 0.01, n=50)
PC.0.01.n100 <- PC.graph(a = 0.01, n=100)
PC.0.01.n200 <- PC.graph(a = 0.01, n=200)


mean(PC.0.1.n50$PC.shd); sd(PC.0.1.n50$PC.shd)
mean(PC.0.1.n100$PC.shd); sd(PC.0.1.n100$PC.shd)
mean(PC.0.1.n200$PC.shd); sd(PC.0.1.n200$PC.shd)

mean(PC.0.05.n50$PC.shd); sd(PC.0.05.n50$PC.shd)
mean(PC.0.05.n100$PC.shd); sd(PC.0.05.n100$PC.shd)
mean(PC.0.05.n200$PC.shd); sd(PC.0.05.n200$PC.shd)

mean(PC.0.01.n50$PC.shd); sd(PC.0.01.n50$PC.shd)
mean(PC.0.01.n100$PC.shd); sd(PC.0.01.n100$PC.shd)
mean(PC.0.01.n200$PC.shd); sd(PC.0.01.n200$PC.shd)



#######################################################################################

plot(Int_eff.list, type="l", col=1,main="Causal Effects Comparison",xlab="node",ylab="causal effect",ylim=c(-2,3))
lines(OBMED_Int_eff.list,col=2)
lines(OBMA_Int_eff.list,col=4)
lines(Int_eff_true.list,col=6)
legend("topright",c("True_distinct","OBMED","OBMA","True"),col=c(1,2,4,6),lty=1)

plot(dist_OBMA,type="l",col=2,main="Absolute Distance of Causal Effects",xlab="node",ylab="distance",ylim=c(-2,5))
lines(dist_OBMED,col=4)
legend("topright",c("OBMED","OBMA"),col=c(2,4),lty=1)


#######################################################################################



#######################################################################################
rm(list=ls())
setwd("~/OneDrive - The Ohio State University/[2019-2023 OSU]/2-2021-Spring/STAT7730 Advanced Comoputational Statistics (Dr. Chang)/Final project/R codes and results")





library(tidyverse)
library(ggplot2)

### Create a boxplot
load("datan50q10p01distinct.RData")
dat <- as.tibble(cbind(OB_MED = dist_OBMED, OB_MA = dist_OBMA))
data<-cbind(c(mean(dist_OBMA,na.rm=T),mean(dist_OBMED,na.rm=T),mean(SHD.list,na.rm=T)),c(sd(dist_OBMA,na.rm=T),sd(dist_OBMED,na.rm=T),sd(SHD.list,na.rm=T)))
colnames(data)<-c("mean","Standard Deviation")
rownames(data)<-c("Distance for OBMA","Distance for OBMED","SHD")
ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_boxplot(aes(distance, fill=Method,linetype=Method), alpha=0.5) + 
  scale_fill_grey() + xlim(0, 4)+
  theme_classic()


load("datan100q10p01distinct.RData")
dat <- as.tibble(cbind(OB_MED = dist_OBMED, OB_MA = dist_OBMA))
data<-cbind(c(mean(dist_OBMA,na.rm=T),mean(dist_OBMED,na.rm=T),mean(SHD.list,na.rm=T)),c(sd(dist_OBMA,na.rm=T),sd(dist_OBMED,na.rm=T),sd(SHD.list,na.rm=T)))
colnames(data)<-c("mean","Standard Deviation")
rownames(data)<-c("Distance for OBMA","Distance for OBMED","SHD")
ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_boxplot(aes(distance, fill=Method,linetype=Method), alpha=0.5) + 
  scale_fill_grey() + xlim(0, 4)+
  theme_classic()

load("datan200q10p01distinct.RData")
dat <- as.tibble(cbind(OB_MED = dist_OBMED, OB_MA = dist_OBMA))
data<-cbind(c(mean(dist_OBMA,na.rm=T),mean(dist_OBMED,na.rm=T),mean(SHD.list,na.rm=T)),c(sd(dist_OBMA,na.rm=T),sd(dist_OBMED,na.rm=T),sd(SHD.list,na.rm=T)))
colnames(data)<-c("mean","Standard Deviation")
rownames(data)<-c("Distance for OBMA","Distance for OBMED","SHD")
ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_boxplot(aes(distance, fill=Method,linetype=Method), alpha=0.5) + 
  scale_fill_grey() + xlim(0, 4)+
  theme_classic()
###



dat <- as.tibble(cbind(OB_MED = dist_OBMED, OB_MA = dist_OBMA))
data<-cbind(c(mean(dist_OBMA,na.rm=T),mean(dist_OBMED,na.rm=T),mean(SHD.list,na.rm=T)),c(sd(dist_OBMA,na.rm=T),sd(dist_OBMED,na.rm=T),sd(SHD.list,na.rm=T)))
colnames(data)<-c("mean","Standard Deviation")
rownames(data)<-c("Distance for OBMA","Distance for OBMED","SHD")


ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_boxplot(aes(distance, fill=Method,linetype=Method), alpha=0.5) + 
  scale_fill_grey() + xlim(0, 4)+
  theme_classic()


# pdf("q10n200.pdf",onefile = T)
#layout(nf)
#par(mfrow = c(2, 2)) 

ggplot() +
  geom_boxplot( aes(SHD.list),alpha=0.5) + 
  scale_fill_grey() +
  theme_classic()+
  ggtitle("SHD for OBMED")

ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_density(aes(distance, linetype=Method, fill=Method), alpha=0.5) + 
  scale_fill_grey() +
  theme_classic()+
  ylim(0, 5)+
  xlim(0, 5)+
  ggtitle("Comparison of Abs Distance of Effects by density plot")

ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_boxplot(aes(distance, fill=Method,linetype=Method), alpha=0.5) + 
  scale_fill_grey() + xlim(0, 4)+
  theme_classic()
  # ggtitle("Comparison of Abs Distance of Effects by boxplot plot")


ggplot(data=dat %>% gather(key = "Method", value = "distance")) +
  geom_histogram(aes(distance, fill=Method), color=1, alpha=0.5,position = "stack") + 
  scale_fill_grey() +
  theme_classic()+
  ggtitle("Comparison of Abs Distance of Effects by Histogram")

# ggplot()
# grid.table(data)

# dev.off()




