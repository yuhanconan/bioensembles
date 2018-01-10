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


###############################
## 2 parameters
###############################

# choose parameters
param2 <- c("K","M")#,"Loo")

Mean2 <- mp[which(names(mp) %in% param2)]
Cov2 <- cov[which(rownames(cov) %in% param2), which(colnames(cov) %in% param2)]
myGrid2 <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
rescale(myGrid2, m=Mean2, C=(Cov2+t(Cov2))/2, dec.type=1)

png(file.path(figs, "2param_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(myGrid2, lwd=3, xlab="log(k)", ylab="log(M)")
dev.off()

# param2alt1 <- c("K","Loo")

# Mean2alt1 <- mp[which(names(mp) %in% param2alt1)]
# Cov2alt1 <- cov[which(rownames(cov) %in% param2alt1), which(colnames(cov) %in% param2alt1)]
# myGrid2alt1 <- createNIGrid(dim=2, type="GHe", level=4,  ndConstruction="sparse")
# rescale(myGrid2alt1, m=Mean2alt1, C=(Cov2alt1+t(Cov2alt1))/2, dec.type=1)

# png(file.path(figs, "2alt1param_dist.png"), width=10, height=8, units="in", res=200)
# par(mfrow=c(1,1))
# plot(myGrid2alt1, lwd=3, xlab="log(Linf)", ylab="log(k)")
# dev.off()

nodes2 <- getNodes(myGrid2)
colnames(nodes2) <- names(Mean2)

weights2 <- getWeights(myGrid2)

saveRDS(nodes2, file.path(sim, "nodes_2param_hogfish.rds"))
saveRDS(weights2, file.path(sim, "weights_2param_hogfish.rds"))


###############################
## 3 parameters
###############################

## choose parameters
param3 <- c("K","M","Loo")

Mean3 <- mp[which(names(mp) %in% param3)]
Cov3 <- cov[which(rownames(cov) %in% param3), which(colnames(cov) %in% param3)]

############## NEW CODE
# Change #1 :  change dim=2 to dim=3 so that it matches exected syntax
# Change #2 : chanve Cov to (Cov+t(Cov))/2 so that its guarunteed to be symmetric

myGrid3 <- createNIGrid(dim=3, type="GHe", level=4,  ndConstruction="sparse")
# !(all(Cov == t(Cov)) & all(eigen(Cov)$values > 0))
rescale(myGrid3, m=Mean3, C=(Cov3+t(Cov3))/2, dec.type=1)

############## END NEW CODE

png(file.path(figs, "3param_dist.png"), width=10, height=8, units="in", res=200)
par(mfrow=c(1,1))
plot(myGrid3, lwd=3, xlab="log(Linf)", ylab="log(k)", zlab="log(M)")
dev.off()

nodes3 <- getNodes(myGrid3)
weights3 <- getWeights(myGrid3)

saveRDS(nodes3, file.path(sim, "nodes_3param_hogfish.rds"))
saveRDS(weights3, file.path(sim, "weights_3param_hogfish.rds"))

####################################
###  simulate truth              ###
####################################
## based on Lutjanus guttatus
plist <- create_lh_list(linf=64.6, vbk=0.21, t0=-0.01, 
						lwa=0.0245, lwb=2.79, 
						M=0.43, 
						M50=34, maturity_input="length",
						S50=c(30), S95=c(35), selex_input="length",
						SigmaF=0.2, SigmaR=0.737)

## generate data
itervec <- 1:50
data <- generate_data(modpath=sim, itervec=itervec, 
						Fdynamics="Endogenous", Rdynamics="AR", 
						lh=plist, 
						Nyears=20, Nyears_comp=20, comp_sample=200,
						init_depl=c(0.10,0.90), 
						seed=itervec,
						nburn=50,
						rewrite=FALSE)

##############################
## Run model - 2 parameters
##############################

itervec <- 1:10

ncores <- 5
registerDoParallel(cores=ncores)	

runstack <- function(savedir, iter, nodes){
	iterpath <- file.path(savedir, itervec[iter])	
	
	## true values
	set.seed(iter)
	vbk_choose <- exp(rnorm(1, mean=Mean2["K"], sd=sqrt(Cov2["K","K"])))
	M_choose <- exp(rnorm(1, mean=Mean2["M"], sd=sqrt(Cov2["M","M"])))
	plist <- create_lh_list(linf=exp(Mean3["Loo"]), vbk=vbk_choose, t0=-0.01,
							lwa=0.053, lwb=2.706,
							M=M_choose,
							M50=16.9, maturity_input="length",
							S50=16.9, S95=42.6, selex_input="length",
							SigmaF=0.2, SigmaR=0.737)

		## generate data
		data <- generate_data(modpath=savedir, itervec=iter, 
						Fdynamics="Endogenous", Rdynamics="AR", 
						lh=plist, 
						Nyears=20, Nyears_comp=20, comp_sample=200,
						init_depl=c(0.10,0.90), 
						seed=iter+1000,
						nburn=50,
						rewrite=FALSE)

		data <- readRDS(file.path(iterpath, "True.rds"))
		data_list <- list("years"=data$years, "LF"=data$LF)	
		res <- lapply(1:nrow(nodes), function(x){
		 		lhinp <- with(plist, 
		 				 create_lh_list(linf=linf, vbk=exp(nodes[x,"K"]), t0=t0,
									lwa=lwa, lwb=lwb,
									M=exp(nodes[x,"M"]),
									M50=M50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=SigmaF, SigmaR=SigmaR))		
		 		out <- run_LIME(modpath=NULL, lh=lhinp, input_data=data_list, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE)	
			return(out)
		})
		saveRDS(res, file.path(iterpath, "res_2param.rds"))
		return(paste0("Ran iter ", iter, " in ", savedir))
}

## predictive stacking
stack_dir <- file.path(sim, "stack")
dir.create(stack_dir, showWarnings=FALSE)
foreach(iter=1:length(itervec), .packages=c('TMB','LIME')) %dopar% runstack(savedir=stack_dir, iter=itervec[iter], nodes=nodes2)

### means only
means_dir <- file.path(sim, "means")
dir.create(means_dir, showWarnings=FALSE)
foreach(iter=1:length(itervec), .packages=c('TMB','LIME') %dopar% runstack(savedir=means_dir, iter=itervec[iter], nodes=Mean2)


spr2 <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})
	stack <- sum(est * weights2)
	true <- gen$SPR

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack))), "Value"="SPR", "Npar"=2)

	return(out)
})
spr2 <- do.call(rbind, spr2)

f2 <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})
	stack <- sum(est * weights2)
	true <- gen$F_t[length(gen$F_t)]

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack))), "Value"="F", "Npar"=2)

	return(out)
})
f2 <- do.call(rbind, f2)

d2 <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))
	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			return(res[[x]]$Report$D_t[length(res[[x]]$Report$D_t)])
	}})
	stack <- sum(est * weights2)
	true <- gen$D_t[length(gen$D_t)]

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack))), "Value"="Relative SB", "Npar"=2)

	return(out)
})
d2 <- do.call(rbind, d2)

ff2 <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			F30 <- uniroot(calc_ref, ages=res[[x]]$Inputs$Data$ages, Mat_a=res[[x]]$Inputs$Data$Mat_a, W_a=res[[x]]$Inputs$Data$W_a, M=res[[x]]$Inputs$Data$M, S_a=res[[x]]$Report$S_a, ref=0.3, interval=c(0,100))$root
			FF30 <- res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)]/F30
			return(FF30)
	}})
	stack <- sum(est * weights2)

	trueF30 <- uniroot(calc_ref, ages=gen$ages, Mat_a=gen$Mat_a, W_a=gen$W_a, M=gen$M, S_a=gen$S_a, ref=0.3, interval=c(0,100))$root
	true <- gen$F_t[length(gen$F_t)]/trueF30

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack))), "Value"="F/F30", "Npar"=2)

	return(out)
})
ff2 <- do.call(rbind, ff2)

fref2 <- lapply(1:length(itervec), function(y){
	res <- readRDS(file.path(sim, itervec[y], "res_2param.rds"))
	gen <- readRDS(file.path(sim, itervec[y], "True.rds"))

	est <- sapply(1:length(res), function(x){
		if(is.null(res[[x]]$df)){
			return(NA)
		} else{
			F30 <- uniroot(calc_ref, ages=res[[x]]$Inputs$Data$ages, Mat_a=res[[x]]$Inputs$Data$Mat_a, W_a=res[[x]]$Inputs$Data$W_a, M=res[[x]]$Inputs$Data$M, S_a=res[[x]]$Report$S_a, ref=0.3, interval=c(0,100))$root
			return(F30)
	}})
	stack <- sum(est * weights2)

	true <- uniroot(calc_ref, ages=gen$ages, Mat_a=gen$Mat_a, W_a=gen$W_a, M=gen$M, S_a=gen$S_a, ref=0.3, interval=c(0,100))$root

	out <- data.frame("Iteration" = y, "Estimate"=c(est, stack), "True"=true, "Estimate_type"=c(rep("Estimate", length(est)), rep("Stacked", length(stack))), "Value"="F30", "Npar"=2)

	return(out)
})
fref2 <- do.call(rbind, fref2)


par2 <- rbind(spr2, d2, f2, fref2, ff2) %>%
		dplyr::mutate("RE"= (Estimate - True)/True)


p <- ggplot(par2) +
	geom_violin(aes(x=Value, y=RE), colour="steelblue", fill="steelblue", scale="width") + 
	facet_grid(Estimate_type ~ .) +
	geom_hline(yintercept = 0, lwd=1.5, alpha=0.7) +
	ylab("Relative error") + 
	theme_lsd()
ggsave(file.path(figs, "RE_est_vs_stacked.png"), p)