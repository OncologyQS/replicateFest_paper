### interp package to interpolate variance for expanded clones 
require(interp)

########### load matrix
countDat=readRDS("countMatrixForSim.rds") 
countDat=as.matrix(countDat) 
countDat <- countDat[rowSums(countDat) > 0, ]### this clone is all zeros
dim(countDat)

###### distinguish peptides from replicates
trimRep=function(x) return(substr(x,1,nchar(x)-2))
peps=sapply(colnames(countDat),trimRep)
peps=sub("Es8","ES8",peps)


### get means and variances for each clone across all conditions
ctMns=apply(countDat,1,mean)
ctVars=apply(countDat,1,var)

### make this a distribution of means and variances, then we can sample from this distribution to get realistic parameters for our simulation
ctMnsS=sort(unique(ctMns))
ctVarsS=sort(tapply(ctVars,ctMns,mean))

### loess version of distribution
set.seed(124323)
loessFit=loess(log(ctVars)~log(ctMns),span=0.5)
ctMnsL=exp(loessFit$x[order(loessFit$x)])
ctVarsL=exp(loessFit$fitted[order(loessFit$x)])

pdf("meanVarPlot.pdf")
plot(ctMns,ctVars,pch=16,xlab="mean",ylab="variance",log="xy")
lines(ctMnsL,ctVarsL,col="blue",lwd=2)
dev.off()

# Create data frames for plotting
data_points <- data.frame(mean = ctMns, variance = ctVars)
fitted_line <- data.frame(mean = ctMnsL, variance = ctVarsL)

# Create plot
p <- ggplot(data_points, aes(x = mean, y = variance)) +
  geom_point(aes(color = "Data points"), size = 2) +
  geom_line(data = fitted_line, aes(color = "LOESS fit"), linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Data points" = "black", "LOESS fit" = "blue")) +
  labs(x = "mean", y = "variance", color = NULL) +
  theme_minimal() +
  theme(legend.position = "top")


pdf("meanVarPlot_ggplot.pdf", width = 5, height = 5)
print(p)
dev.off()
ggsave("meanVarPlot_ggplot.pdf", plot = p, width = 5, height = 4)
ggsave("meanVarPlot_ggplot.png", plot = p, width = 5, height = 4, dpi = 300)


pdf("meanVarPlot2.pdf")
plot(ctMns,ctVars,pch=16,xlab="mean",ylab="variance",log="xy")
lines(ctMnsS,ctVarsS,col="red",lwd=2)
lines(ctMnsL,ctVarsL,col="blue",lwd=2)
dev.off()


#### chose starting means/variances and expansion factors for simulation
meanSim=quantile(ctMns,probs=c(0.5,0.75,0.95))
varSim=quantile(ctVars,probs=c(0.5,0.75,0.95))
expSim=c(2,4,8)

#### parameterization of negative binomial distribution, 
#### size is the dispersion parameter, mu is the mean
#### the size function reparameterizes from var to size. 

size=function(mean,var) return(mean^2/(var-mean) )
sizeSim=size(meanSim,varSim)

##### sequence of functions to simulate in steps


#### first simulate the size to go with the mean and expansion factor

sizeExpSim=function(mean,exp){
  return(size(mean*exp,exp(aspline(
    x=log(ctMnsL),
    y=log(ctVarsL),
    xo=log(mean*exp))$y)))
}

## then simulate data with these parameters
sim=function(mean,size,exp,n=10){
  normDat=matrix(rnbinom(30*n,size=size,mu=mean),nrow=n)
  expDat=matrix(rnbinom(3*n,size=sizeExpSim(mean,exp),mu=mean*exp),nrow=n)
  return(cbind(normDat,expDat))
}

#### a little wrapper to sim from a vector of parameters
simP=function(params,n=10){
  return(sim(mean=params[1],size=params[2],exp=params[3],n=n))
}

#### create the vectors
params=cbind(rep(meanSim,3),rep(sizeSim,3),rep(expSim,each=3))

### simulate data for each set of parameters
getSimData = function(params, countData)
{
  simSet=simP(params[1,])
  for(i in 2:nrow(params)) simSet=rbind(simSet,simP(params[i,]))
  
  levels=apply(cbind(rep(c("L1","L2","L3"),3),rep(c("E1","E2","E3"),each=3)),1,paste,collapse="_")
  clones=apply(cbind(rep(LETTERS[1:10],9),rep(levels,each=10)),1,paste,collapse="_")
  rownames(simSet)=clones
  
  dat=rbind(countDat,simSet)
  dat
}  

# generate 10 datasets 
set.seed(12345)
for (i in 1:10)
{
  dat = getSimData(params, countData)
    #### save result
  saveRDS(dat,paste0("Simulated_matrices/countMatrixFromSim_",i,".rds"))
}

