rm(list=ls())

devtools::install_github("merrillrudd/LIME", dependencies=TRUE)
library(LIME)

library(FishLife)


library(mvQuad)
library(mvtnorm)

library(dplyr)
library(ggplot2)
library(lsd)
library(foreach)
library(doParallel)

## directories
wd <- file.path("C:\\merrill\\bioensembles")

sim <- file.path(wd, "sim")
dir.create(sim, showWarnings=FALSE)

figs <- file.path(wd, "figures")
dir.create(figs, showWarnings=FALSE)

funs <- list.files(file.path(wd, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(wd,"R", funs[x])))


####################################
###  call FishLife               ###
####################################
## Puerto Rico hogfish Lachnolaimu maximus

png(file.path(figs, "FishLife_sim.png"), width=10, height=8, units="in", res=200)
sp <- Plot_taxa( Search_species(Genus="Lachnolaimus", Species="maximus")$match_taxonomy, mfrow=c(2,2) )
dev.off()


## species-level
mp <- sp[[1]]$Mean_pred
cov <- sp[[1]]$Cov_pred

## family-level
mpf <- sp[[3]]$Mean_pred
covf <- sp[[3]]$Cov_pred

### from life history studies, Adams and Rios 2015
study_vec <- c(NA, "Puerto Rico", "Cuba", "GOM", "GOM", "GOM", "GOM", NA, "Dry Tortugas", "FL Keys", "S. Florida", "GOM", "Marquesas Key FL", "S. Florida", "GOM", "GOM", "Florida", "W. Florida", "SE US")
loo_vec <- c(566, 913, 850, 966, NA, 381, 896, 566, 651, 336, 428, 917, 397, 448, 939, NA, NA, 849, NA)
vbk_vec <- c(0.19, 0.08, 0.10, 0.08, NA, 0.56, 0.09, 0.19, 0.11, 0.56, 0.26, 0.08, 0.17, 0.23, 0.07, NA, NA, 0.11, NA)
t0_vec <- c(-0.78, -1.78, -1.38, -1.77, NA, -0.16, -1.98, -2.33, -3.34, -0.19, -0.93, -1.84, -3.74, -1.02, -2.01, NA, NA, -1.33, NA)
M_vec <- c(0.25, 0.13, rep(NA, 16), mean(c(0.16,0.29)))
Lm_vec <- c(NA, 249, NA, 152, 152, NA, NA, 198, NA, NA, NA, NA, NA, 163, 166, 177, NA, NA, mean(c(152, 193)))

lhstudies <- data.frame("Location"=study_vec, "Loo"=loo_vec, "K"=vbk_vec, "t0"=t0_vec, "M"=M_vec, "Lm"=Lm_vec)

a_vec <- c(2.55e-5, 1.52e-2, 2.37e-2, 9.5e-5)
b_vec <- c(2.97, 3.11, 2.95, 2.75)
###############################
## 2 parameters
###############################

# choose parameters
paramMK <- c("K","M")

MeanMK <- mp[which(names(mp) %in% paramMK)]
CovMK <- cov[which(rownames(cov) %in% paramMK), which(colnames(cov) %in% paramMK)]
gridMK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridMK, m=MeanMK, C=(CovMK+t(CovMK))/2, dec.type=1)

## local studies

png(file.path(figs, "MK_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(gridMK, xlab="log(k)", ylab="log(M)", pch=19, cex=1.5)
abline(h=log(lhstudies$M)[-which(is.na(lhstudies$M))][1], col="blue", lty=3)
dev.off()

nodesMK <- getNodes(gridMK)
colnames(nodesMK) <- names(MeanMK)

weightsMK <- getWeights(gridMK)

saveRDS(nodesMK, file.path(sim, "nodes_MK_hogfish.rds"))
saveRDS(weightsMK, file.path(sim, "weights_MK_hogfish.rds"))


# choose parameters
paramML <- c("Loo","M")

MeanML <- mp[which(names(mp) %in% paramML)]
CovML <- cov[which(rownames(cov) %in% paramML), which(colnames(cov) %in% paramML)]
gridML <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridML, m=MeanML, C=(CovML+t(CovML))/2, dec.type=1)

png(file.path(figs, "ML_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(gridML, xlab="log(Loo)", ylab="log(M)", pch=19, cex=1.5)
dev.off()

nodesML <- getNodes(gridML)
colnames(nodesML) <- names(MeanML)

weightsML <- getWeights(gridML)

saveRDS(nodesML, file.path(sim, "nodes_ML_hogfish.rds"))
saveRDS(weightsML, file.path(sim, "weights_ML_hogfish.rds"))


# choose parameters
paramLK <- c("Loo","K")

MeanLK <- mp[which(names(mp) %in% paramLK)]
CovLK <- cov[which(rownames(cov) %in% paramLK), which(colnames(cov) %in% paramLK)]
gridLK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridLK, m=MeanLK, C=(CovLK+t(CovLK))/2, dec.type=1)

png(file.path(figs, "LK_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(gridLK, xlab="log(Loo)", ylab="log(k)", pch=19, cex=1.5)
dev.off()

nodesLK <- getNodes(gridLK)
colnames(nodesLK) <- names(MeanLK)

weightsLK <- getWeights(gridLK)

saveRDS(nodesLK, file.path(sim, "nodes_LK_hogfish.rds"))
saveRDS(weightsLK, file.path(sim, "weights_LK_hogfish.rds"))


###############################
## 3 parameters
###############################

## choose parameters
param3 <- c("K","M","Loo")

Mean3 <- mp[which(names(mp) %in% param3)]
Cov3 <- cov[which(rownames(cov) %in% param3), which(colnames(cov) %in% param3)]

myGrid3 <- createNIGrid(dim=3, type="GHe", level=4,  ndConstruction="sparse")
# !(all(Cov == t(Cov)) & all(eigen(Cov)$values > 0))
rescale(myGrid3, m=Mean3, C=(Cov3+t(Cov3))/2, dec.type=1)

nodes3 <- getNodes(myGrid3)
colnames(nodes3) <- names(Mean3)
weights3 <- getWeights(myGrid3)

saveRDS(nodes3, file.path(sim, "nodes_3param_hogfish.rds"))
saveRDS(weights3, file.path(sim, "weights_3param_hogfish.rds"))



png(file.path(figs, "3param_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(3,1))
plot(nodes3[,"Loo"], nodes3[,"K"], xlab="log(Loo)", ylab="log(k)", pch=19, cex=1.5)
plot(nodes3[,"Loo"], nodes3[,"M"], xlab="log(Loo)", ylab="log(k)", pch=19, cex=1.5)
plot(nodes3[,"M"], nodes3[,"K"], xlab="log(M)", ylab="log(k)", pch=19, cex=1.5)
dev.off()





##############################
## Run model - 2 parameters
##############################


itervec <- 1:100

## predictive stacking
stack_dir <- file.path(sim, "stack")
dir.create(stack_dir, showWarnings=FALSE)

ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)	


start <- Sys.time()
foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=stack_dir, iter=loop, nodes=nodesMK, param=paramMK, mean=Mean3, cov=Cov3, modname="MK", rewrite=FALSE)
end <- Sys.time() - start

stopCluster(cl)


### means only
means_dir <- file.path(sim, "means")
dir.create(means_dir, showWarnings=FALSE)

ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=means_dir, iter=itervec[loop], nodes=t(as.matrix(MeanMK)), param=paramMK, mean=Mean3, cov=Cov3, modname="MK", rewrite=FALSE)

stopCluster(cl)


##############################
## Run model - 3 parameters
##############################


itervec <- 1:100

## predictive stacking
stack_dir <- file.path(sim, "stack")
dir.create(stack_dir, showWarnings=FALSE)

ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)	


start <- Sys.time()
foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=stack_dir, iter=loop, nodes=nodes3, param=param3, mean=Mean3, cov=Cov3, modname="MKL", rewrite=FALSE)
end <- Sys.time() - start

stopCluster(cl)




#########################
## results
#########################


sprMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$SPR

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	if(true != gen2$SPR) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$SPR_t[length(res2[[1]]$Report$SPR_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="SPR", "Npar"=2)
	return(out)
})
sprMK <- do.call(rbind, sprMK) 

fMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$F_t[length(gen$F_t)]

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	if(true != gen2$F_t[length(gen2$F_t)]) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$F_t[length(res2[[1]]$Report$F_t)]
	}
	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="F", "Npar"=2)

	return(out)
})
fMK <- do.call(rbind, fMK)

dMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$D_t[length(res[[x]]$Report$D_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$D_t[length(gen$D_t)]

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	if(true != gen2$D_t[length(gen2$D_t)]) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$D_t[length(res2[[1]]$Report$D_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="Depl", "Npar"=2)
	return(out)
})
dMK <- do.call(rbind, dMK)

ffMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			F30 <- uniroot(calc_ref, ages=res[[x]]$Inputs$Data$ages, Mat_a=res[[x]]$Inputs$Data$Mat_a, W_a=res[[x]]$Inputs$Data$W_a, M=res[[x]]$Inputs$Data$M, S_a=res[[x]]$Report$S_a, ref=0.3, interval=c(0,100))$root
			FF30 <- res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)]/F30
			return(FF30)
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)

	trueF30 <- uniroot(calc_ref, ages=gen$ages, Mat_a=gen$Mat_a, W_a=gen$W_a, M=gen$M, S_a=gen$S_a, ref=0.3, interval=c(0,100))$root
	true <- gen$F_t[length(gen$F_t)]/trueF30

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	trueF302 <-  uniroot(calc_ref, ages=gen2$ages, Mat_a=gen2$Mat_a, W_a=gen2$W_a, M=gen2$M, S_a=gen2$S_a, ref=0.3, interval=c(0,100))$root
	if(true != gen2$F_t[length(gen2$F_t)]/trueF302) stop("Stack and mean true values don't match")
	F302 <- uniroot(calc_ref, ages=res2[[1]]$Inputs$Data$ages, Mat_a=res2[[1]]$Inputs$Data$Mat_a, W_a=res2[[1]]$Inputs$Data$W_a, M=res2[[1]]$Inputs$Data$M, S_a=res2[[1]]$Report$S_a, ref=0.3, interval=c(0,100))$root
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$F_t[length(res2[[1]]$Report$F_t)]/F302
	}


	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="FF30", "Npar"=2)
	return(out)
})
ffMK <- do.call(rbind, ffMK)

frefMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			F30 <- uniroot(calc_ref, ages=res[[x]]$Inputs$Data$ages, Mat_a=res[[x]]$Inputs$Data$Mat_a, W_a=res[[x]]$Inputs$Data$W_a, M=res[[x]]$Inputs$Data$M, S_a=res[[x]]$Report$S_a, ref=0.3, interval=c(0,100))$root
			return(F30)
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)

	trueF30 <- uniroot(calc_ref, ages=gen$ages, Mat_a=gen$Mat_a, W_a=gen$W_a, M=gen$M, S_a=gen$S_a, ref=0.3, interval=c(0,100))$root

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	trueF302 <-  uniroot(calc_ref, ages=gen2$ages, Mat_a=gen2$Mat_a, W_a=gen2$W_a, M=gen2$M, S_a=gen2$S_a, ref=0.3, interval=c(0,100))$root
	if(trueF30 != trueF302) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
	est2 <- uniroot(calc_ref, ages=res2[[1]]$Inputs$Data$ages, Mat_a=res2[[1]]$Inputs$Data$Mat_a, W_a=res2[[1]]$Inputs$Data$W_a, M=res2[[1]]$Inputs$Data$M, S_a=res2[[1]]$Report$S_a, ref=0.3, interval=c(0,100))$root
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="F30", "Npar"=2)

	return(out)
})
frefMK <- do.call(rbind, frefMK)


par2 <- rbind(sprMK, dMK, fMK) %>%
		dplyr::mutate("RE"= (Estimate - True)/True)

tpar2 <- par2 %>%
		dplyr::group_by(Estimate_type, Value) %>%
		dplyr::select(Estimate_type, Value, RE) %>%
		dplyr::summarise_all(funs(mre=median(., na.rm=TRUE), mare=median(abs(.), na.rm=TRUE)))
write.csv(tpar2, file.path(figs, "RE_table_MK.csv"))


p <- ggplot(par2) +
	geom_violin(aes(x=Value, y=RE), colour="steelblue", fill="steelblue", scale="width") + 
	facet_grid(Estimate_type ~ .) +
	geom_hline(yintercept = 0, lwd=1.5, alpha=0.7) +
	ylab("Relative error") + 
	theme_lsd() +
	coord_cartesian(ylim=c(-1,10))
ggsave(file.path(figs, "RE_est_vs_stacked.png"), p)