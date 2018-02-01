rm(list=ls())

###############################
### packages 				###
###############################
devtools::install_github("merrillrudd/LIME", dependencies=TRUE, ref="master")
library(LIME)

library(FishLife)


library(mvQuad)
library(mvtnorm)

library(dplyr)
library(ggplot2)
library(lsd)
library(foreach)
library(doParallel)

###############################
### directories 			###
###############################
wd <- file.path("C:\\merrill\\bioensembles")

sim <- file.path(wd, "sim")
dir.create(sim, showWarnings=FALSE)

figs <- file.path(wd, "figures")
dir.create(figs, showWarnings=FALSE)

funs <- list.files(file.path(wd, "R"))
ignore <- sapply(1:length(funs), function(x) source(file.path(wd,"R", funs[x])))

res_dir <- file.path(sim, "results")
dir.create(res_dir, showWarnings=FALSE)


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

## genus-level
mpg <- sp[[2]]$Mean_pred
covg <- sp[[2]]$Cov_pred

## family-level
mpf <- sp[[3]]$Mean_pred
covf <- sp[[3]]$Cov_pred

## calculate maximum age
tmax_all <- ceiling(quantile(rlnorm(1000, mean=mp["tmax"], sd=sqrt(cov["tmax","tmax"])), prob=0.95))

###############################
### M and K
###############################

# choose parameters
paramMK <- c("K","M")

## means and grid
## species level
MeanMK <- mp[which(names(mp) %in% paramMK)]
CovMK <- cov[which(rownames(cov) %in% paramMK), which(colnames(cov) %in% paramMK)]
gridMK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridMK, m=MeanMK, C=(CovMK+t(CovMK))/2, dec.type=1)

## family level
MeanfMK <- mpf[which(names(mpf) %in% paramMK)]
CovfMK <- covf[which(rownames(covf) %in% paramMK), which(colnames(covf) %in% paramMK)]
gridfMK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridfMK, m=MeanfMK, C=(CovfMK+t(CovfMK))/2, dec.type=1)


## get nodes and weights
nodesMK <- getNodes(gridMK)
colnames(nodesMK) <- names(MeanMK)
weightsMK <- getWeights(gridMK)
saveRDS(nodesMK, file.path(sim, "nodes_MK_hogfish.rds"))
saveRDS(weightsMK, file.path(sim, "weights_MK_hogfish.rds"))

## family level
nodesfMK <- getNodes(gridfMK)
colnames(nodesfMK) <- names(MeanfMK)
weightsfMK <- getWeights(gridfMK)
saveRDS(nodesfMK, file.path(sim, "nodes_MK_labridae.rds"))
saveRDS(weightsfMK, file.path(sim, "weights_MK_labridae.rds"))

png(file.path(figs, "MK_dist.png"), width=10, height=8, units="in", res=200)
col <- rev(brewer.pal(3,"Set1"))
par(mfrow=c(1,1), mar=c(4,5,2,2))
plot(exp(nodesfMK[,"K"]), exp(nodesfMK[,"M"]), xlab="von Bertalanffy k", ylab="Natural mortality", pch=19, cex=2, col=paste0(col[2],"90"), cex.lab=2, cex.axis=1.5)
# points(exp(nodesgMK[,"K"]), exp(nodesgMK[,"M"]), pch=19, cex=2, col=paste0(col[2],"90"))
points(exp(nodesMK[,"K"]), exp(nodesMK[,"M"]), pch=19, cex=2, col=paste0(col[3],"90"))
legend("topleft", legend=c("Family", "Species"), col=paste0(col[c(2,3)],"90"), pch=19, cex=2)
# legend("topleft", legend=c("Family", "Genus", "Species"), col=paste0(col,"90"), pch=19, cex=2)
dev.off()

###############################
## 3 parameters
###############################

## choose parameters
param3 <- c("K","M","Loo")

## means
Mean3 <- mp[which(names(mp) %in% param3)]
Cov3 <- cov[which(rownames(cov) %in% param3), which(colnames(cov) %in% param3)]

Meanf3 <- mpf[which(names(mpf) %in% param3)]
Covf3 <- covf[which(rownames(covf) %in% param3), which(colnames(covf) %in% param3)]


##############################
## Run model - 2 parameters
##############################

itervec <- 1:100


### effort dynamics -- endogenous
dyn_dir <- file.path(res_dir, "harvestdyn")
dir.create(dyn_dir, showWarnings=FALSE)

##############
## species
#############
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)	

start <- Sys.time()
foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=dyn_dir, iter=loop, tmax=tmax_all, nodes=nodesMK, param=paramMK, mean=Mean3, cov=Cov3, modname="MK", level="species", input_data=FALSE, Fscenario="harvestdyn", rewrite=FALSE)
end <- Sys.time() - start

stopCluster(cl)

##############
## family
#############
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)	

start <- Sys.time()
foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=dyn_dir, iter=loop, tmax=tmax_all, nodes=nodesfMK, param=paramMK, mean=Mean3, cov=Cov3, modname="MK_family", input_data=FALSE, Fscenario="harvestdyn", rewrite=FALSE)
end <- Sys.time() - start

stopCluster(cl)


#########################
## Results - 2parameters
#########################
	### visual summary of predictive stacking method

## good seeds for visualizing: 1212, 111, 4444
png(file.path(figs, "Predictive_stacking_MK_effortdyn.png"), height=10, width=15, res=200, units="in")
stacking_figure(res_dir=file.path(res_dir,"harvestdyn"), itervec=1:100, seed=555)
dev.off()


sprMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_MK.rds"))
	mres <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_IterTrue.rds"))
	gen <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$SPR

	if(is.null(mres$df)){
		est2 <- NA
	} else {
		est2 <- mres$Report$SPR_t[length(mres$Report$SPR_t)]
	}

	if(is.null(ires$df)){
		est3 <- NA
	} else {
		est3 <- ires$Report$SPR_t[length(ires$Report$SPR_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(stack, est2, est3), "True"=true, "Estimate_type"=c(rep("Stacked", length(stack)), rep("Means", length(est2)), rep("Truth", length(est3))), "Value"="SPR", "Npar"=2)
	return(out)
})
sprMK <- do.call(rbind, sprMK) 

fMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_MK.rds"))
	mres <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_IterTrue.rds"))
	gen <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$F_t[length(gen$F_t)]

	if(is.null(mres$df)){
		est2 <- NA
	} else {
		est2 <- mres$Report$F_t[length(mres$Report$F_t)]
	}

	if(is.null(ires$df)){
		est3 <- NA
	} else {
		est3 <- ires$Report$F_t[length(ires$Report$F_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(stack, est2, est3), "True"=true, "Estimate_type"=c(rep("Stacked", length(stack)), rep("Means", length(est2)), rep("Truth", length(est3))), "Value"="F", "Npar"=2)
	return(out)
})
fMK <- do.call(rbind, fMK) 

ffrefMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_MK.rds"))
	mres <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_IterTrue.rds"))
	gen <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			F <- res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)]
			F30 <- uniroot(calc_ref, 
				ages=0:(length(res[[x]]$Inputs$Data$L_a)-1),
				Mat_a=res[[x]]$Inputs$Data$Mat_a,
				W_a=res[[x]]$Inputs$Data$W_a,
				M=res[[x]]$Inputs$Data$M,
				S_a=res[[x]]$Report$S_a,
				lower=0,
				upper=200,
				ref=0.3)$root
			return(F/F30)
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	trueF30 <- uniroot(calc_ref, 
				ages=gen$ages,
				Mat_a=gen$Mat_a,
				W_a=gen$W_a,
				M=gen$M,
				S_a=gen$S_a,
				lower=0,
				upper=200,
				ref=0.3)$root
	true <- gen$F_t[length(gen$F_t)]/trueF30

	if(is.null(mres$df)){
		est2 <- NA
	} else {
		mresF30 <- uniroot(calc_ref, 
				ages=0:(length(mres$Inputs$Data$L_a)-1),
				Mat_a=mres$Inputs$Data$Mat_a,
				W_a=mres$Inputs$Data$W_a,
				M=mres$Inputs$Data$M,
				S_a=mres$Report$S_a,
				lower=0,
				upper=200,
				ref=0.3)$root
		est2 <- mres$Report$F_t[length(mres$Report$F_t)]/mresF30
	}

	if(is.null(ires$df)){
		est3 <- NA
	} else {
		iresF30 <- uniroot(calc_ref, 
				ages=0:(length(ires$Inputs$Data$L_a)-1),
				Mat_a=ires$Inputs$Data$Mat_a,
				W_a=ires$Inputs$Data$W_a,
				M=ires$Inputs$Data$M,
				S_a=ires$Report$S_a,
				lower=0,
				upper=200,
				ref=0.3)$root
		est3 <- ires$Report$F_t[length(ires$Report$F_t)]/iresF30
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(stack, est2, est3), "True"=true, "Estimate_type"=c(rep("Stacked", length(stack)), rep("Means", length(est2)), rep("Truth", length(est3))), "Value"="FF30", "Npar"=2)
	return(out)
})
ffrefMK <- do.call(rbind, ffrefMK) 


dMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_MK.rds"))
	mres <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_IterTrue.rds"))
	gen <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$D_t[length(res[[x]]$Report$D_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$D_t[length(gen$D_t)]

	if(is.null(mres$df)){
		est2 <- NA
	} else {
		est2 <- mres$Report$D_t[length(mres$Report$D_t)]
	}

	if(is.null(ires$df)){
		est3 <- NA
	} else {
		est3 <- ires$Report$D_t[length(ires$Report$D_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(stack, est2, est3), "True"=true, "Estimate_type"=c(rep("Stacked", length(stack)), rep("Means", length(est2)), rep("Truth", length(est3))), "Value"="RelSB", "Npar"=2)
	return(out)
})
dMK <- do.call(rbind, dMK) 

rMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_MK.rds"))
	mres <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "res_IterTrue.rds"))
	gen <- readRDS(file.path(res_dir, "harvestdyn", itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$R_t[length(res[[x]]$Report$R_t)])
	}})
	stack <- sum(est * weightsMK, na.rm=TRUE)
	true <- gen$R_t[length(gen$R_t)]

	if(is.null(mres$df)){
		est2 <- NA
	} else {
		est2 <- mres$Report$R_t[length(mres$Report$R_t)]
	}

	if(is.null(ires$df)){
		est3 <- NA
	} else {
		est3 <- ires$Report$R_t[length(ires$Report$R_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(stack, est2, est3), "True"=true, "Estimate_type"=c(rep("Stacked", length(stack)), rep("Means", length(est2)), rep("Truth", length(est3))), "Value"="Rec", "Npar"=2)
	return(out)
})
rMK <- do.call(rbind, rMK) 

par2 <- rbind(dMK, sprMK, fMK, ffrefMK, rMK) %>%
		dplyr::mutate("RE"= (Estimate - True)/True)

tpar2 <- par2 %>%
		dplyr::group_by(Estimate_type, Value) %>%
		dplyr::select(Estimate_type, Value, RE) %>%
		dplyr::summarise_all(funs(mre=median(., na.rm=TRUE), mare=median(abs(.), na.rm=TRUE)))
write.csv(tpar2, file.path(figs, "RE_table_MK.csv"))


col <- brewer.pal(4, "Set1")[c(2,1,3)]
p <- ggplot(par2) +
	geom_violin(aes(x=Value, y=RE, colour=Estimate_type, fill=Estimate_type), scale="width") + 
	scale_colour_manual(values=col) +
	scale_fill_manual(values=col) +
	facet_grid(Estimate_type ~ .) +
	geom_hline(yintercept = 0, lwd=1.5, alpha=0.7) +
	ylab("Relative error") + 
	theme_lsd() +
	coord_cartesian(ylim=c(-1,5))
ggsave(file.path(figs, "RE_est_vs_stacked_MK.png"), p)




######################################
## case study - Puerto Rican hogfish
######################################

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
## 3 parameter grid
###############################

## grid
myGrid3 <- createNIGrid(dim=3, type="GHe", level=4,  ndConstruction="sparse")
rescale(myGrid3, m=Mean3, C=(Cov3+t(Cov3))/2, dec.type=1)

## get  nodes and weights
nodes3 <- getNodes(myGrid3)
colnames(nodes3) <- names(Mean3)
weights3 <- getWeights(myGrid3)

saveRDS(nodes3, file.path(sim, "nodes_3param_hogfish.rds"))
saveRDS(weights3, file.path(sim, "weights_3param_hogfish.rds"))

png(file.path(figs, "3param_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(3,1))
plot(exp(nodes3[,"Loo"]), exp(nodes3[,"K"]), xlab="von Bertalanffy asymptotic length", ylab="von Bertalanffy k", pch=19, cex=1.5)
plot(exp(nodes3[,"Loo"]), exp(nodes3[,"M"]), xlab="von Bertalanffy asymptotic length", ylab="von Bertalanffy k", pch=19, cex=1.5)
plot(exp(nodes3[,"M"]), exp(nodes3[,"K"]), xlab="Natural mortality", ylab="von Bertalanffy k", pch=19, cex=1.5)
dev.off()




##############################
## Run model - 3 parameters
##############################


# itervec <- 1:5

# ## predictive stacking
# res_dir <- file.path(sim, "stack")
# dir.create(res_dir, showWarnings=FALSE)

# ncores <- 5
# cl <- makeCluster(ncores)
# registerDoParallel(cl)	


# start <- Sys.time()
# foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
# 		runstack(savedir=res_dir, iter=loop, nodes=nodes3, param=param3, mean=Mean3, cov=Cov3, modname="MKL", rewrite=FALSE)
# end <- Sys.time() - start

# stopCluster(cl)






sprMKL <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})
	stack <- sum(est * weights3, na.rm=TRUE)
	true <- gen$SPR

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	# if(true != gen2$SPR) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$SPR_t[length(res2[[1]]$Report$SPR_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="SPR", "Npar"=2)
	return(out)
})
sprMKL <- do.call(rbind, sprMKL) 

fMKL <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})
	stack <- sum(est * weights3, na.rm=TRUE)
	true <- gen$F_t[length(gen$F_t)]

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	# if(true != gen2$F_t[length(gen2$F_t)]) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$F_t[length(res2[[1]]$Report$F_t)]
	}
	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="F", "Npar"=2)

	return(out)
})
fMKL <- do.call(rbind, fMKL)

dMKL <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$D_t[length(res[[x]]$Report$D_t)])
	}})
	stack <- sum(est * weights3, na.rm=TRUE)
	true <- gen$D_t[length(gen$D_t)]

	res2 <- readRDS(file.path(means_dir, itervec[y], "res_MK.rds"))
	gen2 <- readRDS(file.path(means_dir, itervec[y], "True.rds"))
	# if(true != gen2$D_t[length(gen2$D_t)]) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$D_t[length(res2[[1]]$Report$D_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="Depl", "Npar"=2)
	return(out)
})
dMKL <- do.call(rbind, dMKL)


par3 <- rbind(sprMKL, dMKL, fMKL) %>%
		dplyr::mutate("RE"= (Estimate - True)/True)

tpar3 <- par3 %>%
		dplyr::group_by(Estimate_type, Value) %>%
		dplyr::select(Estimate_type, Value, RE) %>%
		dplyr::summarise_all(funs(mre=median(., na.rm=TRUE), mare=median(abs(.), na.rm=TRUE)))
write.csv(tpar3, file.path(figs, "RE_table_MKL.csv"))


p <- ggplot(par3) +
	geom_violin(aes(x=Value, y=RE), colour="steelblue", fill="steelblue", scale="width") + 
	facet_grid(Estimate_type ~ .) +
	geom_hline(yintercept = 0, lwd=1.5, alpha=0.7) +
	ylab("Relative error") + 
	theme_lsd() +
	coord_cartesian(ylim=c(-1,10))
ggsave(file.path(figs, "RE_est_vs_stacked.png"), p)
