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

##############################
## maximum age from FishLife
##############################

tmax_all <- ceiling(quantile(rlnorm(1000, mean=mp["tmax"], sd=sqrt(cov["tmax","tmax"])), prob=0.95))

###############################
## 2 parameters
###############################

# choose parameters
paramMK <- c("K","M")

MeanMK <- mp[which(names(mp) %in% paramMK)]
CovMK <- cov[which(rownames(cov) %in% paramMK), which(colnames(cov) %in% paramMK)]
gridMK <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(gridMK, m=MeanMK, C=(CovMK+t(CovMK))/2, dec.type=1)

nodesMK <- getNodes(gridMK)
colnames(nodesMK) <- names(MeanMK)

weightsMK <- getWeights(gridMK)
# weightsMK_rescale <- weightsMK + 0.8
# weightsMK_rescale2 <- weightsMK_rescale
# weightsMK_rescale2[which(weightsMK < 0)] <- weightsMK_rescale[which(weightsMK < 0)] * 0.5
# weightsMK_rescale2 <- floor(100*weightsMK_rescale2)

saveRDS(nodesMK, file.path(sim, "nodes_MK_hogfish.rds"))
saveRDS(weightsMK, file.path(sim, "weights_MK_hogfish.rds"))

# plot(nodesMK[,1], nodesMK[,2], pch=16, cex=2, col=paste0("#111111", weightsMK_rescale2))


png(file.path(figs, "MK_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(gridMK, xlab="log(k)", ylab="log(M)", pch=19, cex=1.5)
# rug(log(vbk_vec), side=1, col="steelblue", lwd=3)
# rug(log(M_vec), side=2, col="tomato", lwd=3)
# abline(h=log(lhstudies$M)[-which(is.na(lhstudies$M))][1], col="blue", lty=3)
dev.off()


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
plot(nodes3[,"Loo"], nodes3[,"M"], xlab="log(Loo)", ylab="log(M)", pch=19, cex=1.5)
plot(nodes3[,"M"], nodes3[,"K"], xlab="log(M)", ylab="log(k)", pch=19, cex=1.5)
dev.off()





##############################
## Run model - 2 parameters
##############################


itervec <- 1:5

## predictive stacking
res_dir <- file.path(sim, "results")
dir.create(res_dir, showWarnings=FALSE)

equil_dir <- file.path(res_dir, "equil")
dir.create(equil_dir, showWarnings=FALSE)


ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)	


start <- Sys.time()
foreach(loop=1:length(itervec), .packages=c('TMB','LIME')) %dopar% 
		runstack(savedir=equil_dir, iter=loop, tmax=tmax_all, nodes=nodesMK, param=paramMK, mean=Mean3, cov=Cov3, modname="MK", input_data=FALSE, Fscenario="equil", rewrite=TRUE)
end <- Sys.time() - start

stopCluster(cl)



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




#########################
## results
#########################
	### figure 3 - visual summary of predictive stacking method
png(file.path(figs, "Predictive_stacking.png"), height=8, width=10, res=200, units="in")
par(mfrow=c(5,3), mar=c(0,3,0,0), omi=c(0.5,0.5,0.5,0.5))
col <- brewer.pal(4, "Set1")
plot_iter <- 1:5
for(i in 1:length(plot_iter)){
	choose_iter=plot_iter[i]
	choose_dir <- file.path(res_dir, "F_Constant", choose_iter)
	res <- readRDS(file.path(choose_dir, "res_MK.rds"))
	true <- readRDS(file.path(choose_dir, "True.rds"))
	mres <- readRDS(file.path(choose_dir, "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(choose_dir, "res_IterTrue.rds"))


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


	plot(x=exp(nodesMK[,"K"]), y=exp(nodesMK[,"M"]), pch=16, xlim=c(0.9*min(exp(nodesMK[,"K"])), 1.1*max(exp(nodesMK[,"K"]))), ylim=c(0.9*min(exp(nodesMK[,"M"])), 1.1*max(exp(nodesMK[,"M"]))), ylab="", xlab="", xaxt="n", las=2, xaxs="i", yaxs="i")
	points(true$vbk, true$M, col=col[1], pch=19, cex=2)
	points(exp(MeanMK["K"]), exp(MeanMK["M"]), col=col[3], pch=19, cex=2)
	if(i==length(plot_iter)){
		axis(1)
		mtext(side=1, "von Bertalanffy k", line=3)
	}
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
	lines(x=1:length(true$D_t), y=true$D_t, col=col[1], lwd=3)
	lines(x=1:length(ires$Report$D_t), y=ires$Report$D_t, col=col[2], lwd=3)
	lines(x=1:length(mres$Report$D_t), y=mres$Report$D_t, col=col[3], lwd=3)
	## predictive stacking
	lines(x=1:length(stack), y=stack, col=col[4], lwd=3)


	plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0,2), ylim=c(0,8), ylab="", xaxt="n", yaxt="n", las=2)
	xx <- seq(0,2,by=0.001)
	axis(2, las=2)
	stack_dens <- sapply(1:length(res), function(x){
		sub <- find %>% filter(Node==x)
		term <- sub %>% filter(Year==max(Year))
		if(is.na(term[,"MLE"])==FALSE & is.na(term[,"SE"])==FALSE){
			dxx <- dlnorm(xx, mean=log(term[,"MLE"]), sd=term[,"SE"])
			polygon(x=c(xx, rev(xx)), y=c(dxx, rep(0, length(dxx))), col="#AAAAAA80", border="#222222")
			# xx <- density(rlnorm(10000, mean=log(term[,"MLE"]), sd=term[,"SE"]))
			# polygon(x=c(xx$x, rev(xx$x)), y=c(xx$y, rep(0,length(xx$x))), col="#AAAAAA80", border="#222222")
		} else { dxx <- rep(NA, length(xx))}
		return(dxx)
	})
	wdens <- sapply(1:length(res), function(x){
		return(weightsMK[x] * stack_dens[,x])
	})
	ysum <- sapply(1:nrow(wdens), function(x){
		sub <- wdens[x,]
		return(sum(sub, na.rm=TRUE))
	})
	polygon(x=c(xx, rev(xx)), y=c(ysum, rep(0, length(xx))), col=paste0(col[4],"80"), border=col[4])
	abline(v=stack[length(stack)], col=col[4], lwd=4)
	abline(v=true$D_t[length(true$D_t)], col=col[1], lwd=3)
	abline(v=ires$Report$D_t[length(stack)], col=col[2], lwd=3)
	abline(v=mres$Report$D_t[length(stack)], col=col[3], lwd=3)

	if(i==1) legend("topright", legend=c("True", "Run at Truth", "FishLife Means", "Quadrature Nodes", "Predictive stacking"), col=c(col[1], col[2], col[3], "#AAAAAA80", col[4]), lwd=3)
	if(i==length(plot_iter)){
		axis(1)
		mtext(side=1, "Terminal year estimate", line=3)
	}
}
dev.off()


sprMK <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(res_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
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
	res <- readRDS(file.path(res_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
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
	res <- readRDS(file.path(res_dir, itervec[y], "res_MK.rds"))
	gen <- readRDS(file.path(res_dir, itervec[y], "True.rds"))
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
