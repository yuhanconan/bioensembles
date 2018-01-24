rm(list=ls())

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

###############################
## 2 parameters
###############################

# choose parameters
paramMK <- c("K","M")

MeanMK <- mp[which(names(mp) %in% paramMK)]
CovMK <- cov[which(rownames(cov) %in% paramMK), which(colnames(cov) %in% paramMK)]
gridMK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridMK, m=MeanMK, C=(CovMK+t(CovMK))/2, dec.type=1)


png(file.path(figs, "MK_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(gridMK, xlab="log(k)", ylab="log(M)", pch=19, cex=1.5)
# abline(h=log(lhstudies$M)[-which(is.na(lhstudies$M))][1], col="blue", lty=3)
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
	### figure 3 - visual summary of predictive stacking method
png(file.path(figs, "Predictive_stacking.png"), height=8, width=10, res=200, units="in")
par(mfrow=c(3,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
col <- brewer.pal(3, "Set1")
plot_iter <- 1:3
for(i in 1:length(plot_iter)){
	choose_iter=plot_iter[i]
	res <- readRDS(file.path(stack_dir, choose_iter, "res_MK.rds"))
	find <- lapply(1:length(res), function(x){	

		if(is.null(res[[x]]$df)==FALSE){
			sub <- res[[x]]$Report
			mle <- sub$D_t		

			sdrep <- summary(res[[x]]$Sdreport)
			sd <- sdrep[which(rownames(sdrep)=="lD_t"),2]
		}
		if(is.null(res[[x]]$df)){
			mle <- rep(NA, 20)
			sd <- rep(NA, 20)
		}		
			df <- data.frame("Node"=x, "Year"=1:length(mle), "MLE"=mle, "SE"=sd)
			return(df)
	})
	find <- do.call(rbind, find)

	yrs <- unique(find$Year)
	stack <- sapply(yrs, function(x){
		yrsub <- find %>% filter(Year==yrs[x])
		stack <- sum(yrsub[,"MLE"] * weightsMK, na.rm=TRUE)
	})	

	# quant <- sapply(1:nrow(find), function(x){
	# 	return(quantile(find[x,], prob=c(0.05,0.5,0.95)))
	# })
	plot(x=1, y=1, type="n", xlim=c(1,20), ylim=c(0,2), ylab="", xlab="", xaxt="n", las=2, xaxs="i", yaxs="i")
	if(i==length(plot_iter)){
		axis(1)
		mtext(side=1, "Year", line=3)
	}
	## estimates for each node
	for(x in 1:length(res)){
		if(is.null(res[[x]]$df) == FALSE) sub <- res[[x]]$Report
		lines(x=1:length(sub$D_t), y=sub$D_t, lwd=3, col="#AAAAAA80")
	}	

	# polygon(x=c(1:ncol(quant), ncol(quant):1), y=c(quant[1,], rev(quant[3,])), col="#55555550", border=NA)
	# lines(x=1:ncol(quant), y=quant[2,], lwd=3, col="black")
	true <- readRDS(file.path(stack_dir, choose_iter, "True.rds"))
	mres <- readRDS(file.path(means_dir, choose_iter, "res_MK.rds"))
	lines(x=1:length(true$D_t), y=true$D_t, col=col[1], lwd=3)
	lines(x=1:length(mres[[1]]$Report$D_t), y=mres[[1]]$Report$D_t, col=col[2], lwd=3)
	## predictive stacking
	lines(x=1:length(stack), y=stack, col=col[3], lwd=3)
	if(i==1) legend("topleft", legend=c("True", "FishLife Means", "Quadrature Nodes", "Predictive stacking"), col=c(col[1], col[2], "#AAAAAA80", col[3]), lwd=3)

	plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(-1.5,2.5), ylim=c(0,2), ylab="", xaxt="n", yaxt="n", las=2)
	axis(4, las=2)
	if(i==length(plot_iter)){
		axis(1)
		mtext(side=1, "Relative spawning biomass\n in terminal year", line=4)
	}
	for(x in 1:length(res)){
		sub <- find %>% filter(Node==x)
		term <- sub %>% filter(Year==max(Year))
		if(is.na(term[,"MLE"])==FALSE & is.na(term[,"SE"])==FALSE){
			xx <- density(rnorm(10000, mean=(term[,"MLE"]), sd=term[,"SE"]))
			polygon(x=c(xx$x, rev(xx$x)), y=c(xx$y, rep(0,length(xx$x))), col="#AAAAAA80", border="#222222")
		}
	}
	abline(v=stack[length(stack)], col=col[3], lwd=4)
	abline(v=true$D_t[length(true$D_t)], col=col[1], lwd=3)
	abline(v=mres[[1]]$Report$D_t[length(stack)], col=col[2], lwd=3)
}
mtext(side=2, "Relative spawning biomass", outer=TRUE, line=4)
mtext(side=4, "Density", outer=TRUE, line=4)
dev.off()


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
	# if(true != gen2$SPR) stop("Stack and mean true values don't match")
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
	# if(true != gen2$F_t[length(gen2$F_t)]) stop("Stack and mean true values don't match")
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
	# if(true != gen2$D_t[length(gen2$D_t)]) stop("Stack and mean true values don't match")
	if(is.null(res2[[1]]$df)){
		est2 <- NA
	} else {
		est2 <- res2[[1]]$Report$D_t[length(res2[[1]]$Report$D_t)]
	}

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack, est2), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack)), rep("Estimate_at_mean", length(est2))), "Value"="Depl", "Npar"=2)
	return(out)
})
dMK <- do.call(rbind, dMK)



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
ggsave(file.path(figs, "RE_est_vs_stacked_MK.png"), p)


sprMKL <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
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
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
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
	res <- readRDS(file.path(stack_dir, itervec[y], "res_MKL.rds"))
	gen <- readRDS(file.path(stack_dir, itervec[y], "True.rds"))
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
