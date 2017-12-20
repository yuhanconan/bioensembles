rm(list=ls())

devtools::install_github("merrillrudd/LIME", dependencies=TRUE)#, ref="multifleet")
library(LIME)

library(FishLife)


library(mvQuad)
library(mvtnorm)

wd <- file.path("C:\\merrill\\bioensembles")

Quadrature = function( myGrid, f, ... ){
  Nodes = getNodes(myGrid)
  Weights = getWeights(myGrid)
  # Compute integral
  Integral = rep(NA, nrow(Nodes))
  for( i in seq_along(Integral) ){
    Integral[i] = f( Nodes[i,], ... ) * Weights[i]
  }
  return( sum(Integral) )
}

Sampling = function( dim, min=-Inf, max=Inf, SampVar, f, n=1000, ... ){
  Nodes = matrix( runif(n=n*dim, min=min, max=max), ncol=dim )
  Weights = rep( (max-min)^dim/n, length=n )
  # Compute integral
  Integral = rep(NA, nrow(Nodes))
  for( i in seq_along(Integral) ){
    Integral[i] = f( Nodes[i,], ... ) * Weights[i]
  }
  return( sum(Integral) )
}

## based on Lutjanus guttatus
plist <- create_lh_list(linf=64.6, vbk=0.21, t0=-0.01, 
						lwa=0.0245, lwb=2.79, 
						M=0.43, 
						M50=34, maturity_input="length",
						S50=c(30), S95=c(35), selex_input="length",
						SigmaF=0.2, SigmaR=0.737)

data <- generate_data(modpath=NULL, itervec=1, 
						Fdynamics="Endogenous", Rdynamics="AR", 
						lh=plist, 
						Nyears=10, Nyears_comp=10, comp_sample=200,
						init_depl=0.5,
						seed=4, 
						nburn=50)

## use FishLife
sp <- Plot_taxa( Search_species(Genus="Lutjanus", Species="guttatus")$match_taxonomy, mfrow=c(2,2) )

mp <- sp[[1]]$Mean_pred
cov <- sp[[1]]$Cov_pred

## choose parameters
param <- c("K","M")

Mean <- mp[which(names(mp) %in% param)]
Cov <- cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)]
f = dmvnorm

myGrid <- createNIGrid(dim=2, type="GHe", level=5,  ndConstruction="sparse")
rescale(myGrid, m=Mean, C=Cov, dec.type=1)

plot(myGrid)
quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Quadrature(f=f, mean=Mean, sigma=Cov, myGrid)
Sampling(f=f, mean=Mean, sigma=Cov, dim=2, min=-2, max=4, n=10000)

nodes <- getNodes(myGrid)
weights <- getWeights(myGrid)

## multifleet
# LFarray <- array(data$LF, dim=c(nrow(data$LF), ncol(data$LF), plist$nfleets))
# rownames(LFarray) <- rownames(data$LF)
# colnames(LFarray) <- colnames(data$LF)
# data_list <- list("years"=data$years, "LF"=LFarray)

## no multifleet
data_list <- list("years"=data$years, "LF"=data$LF)

res <- list()
for(i in 1:nrow(nodes)){
	lhinp <- with(plist, create_lh_list(linf=linf, vbk=exp(nodes[i,1]), t0=t0,
					lwa=lwa, lwb=lwb,
					M=exp(nodes[i,2]),
					M50=M50, maturity_input="length",
					S50=SL50, S95=SL95, selex_input="length",
					SigmaF=SigmaF, SigmaR=SigmaR))

	# res[[i]] <- run_LIME(modpath=NULL, lh=lhinp, input_data=data_list,
	# 				est_param="log_sigma_R", data_avail="LC")
	res[[i]] <- run_LIME(modpath=NULL, lh=lhinp, input_data=data_list,
					est_sigma="log_sigma_R", data_avail="LC")


	inp <- res[[i]]$Inputs
	inp$Data$n_f <- 1
	rep <- res[[i]]$Report
	sdrep <- res[[i]]$Sdreport
	df <- res[[i]]$df

	plot_output(Inputs=inp, Report=rep, Sdreport=sdrep, lh=lhinp, set_ylim=list("Fish"=c(0,3)), True=data, true_years=data$years, lc_years=data$years, all_years=data$years)
}

saveRDS(res, file.path(wd, "res_master.rds"))
saveRDS(nodes, file.path(wd, "nodes_master.rds"))
saveRDS(weights, file.path(wd, "weights_master.rds"))

spr_mle <- sapply(1:length(res), function(x){
	if(is.null(res[[x]]$df)){
		return(NA)
	} else{
		return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)])
	}})

f_mle <- sapply(1:length(res), function(x){
	if(is.null(res[[x]]$df)){
		return(NA)
	} else{
		return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)])
	}})

spr_mle_wt <- sapply(1:length(res), function(x){
	if(is.null(res[[x]]$df)){
		return(NA)
	} else{
		return(res[[x]]$Report$SPR_t[length(res[[x]]$Report$SPR_t)] * exp(weights[x]))
	}})

f_mle_wt <- sapply(1:length(res), function(x){
	if(is.null(res[[x]]$df)){
		return(NA)
	} else{
		return(res[[x]]$Report$F_t[length(res[[x]]$Report$F_t)] * exp(weights[x]))
	}})



par(mfrow=c(2,2))
boxplot(spr_mle, col="gray", main="SPR", ylim=c(0,1))
abline(h=data$SPR, col="red", lwd=3)

boxplot(f_mle, col="gray", main="F", ylim=c(0, 2))
abline(h=data$F_t[length(data$F_t)], col="red", lwd=3)

boxplot(spr_mle_wt, col="steelblue", ylim=c(0,1))
abline(h=data$SPR, col="red", lwd=3)

boxplot(f_mle_wt, col="steelblue", ylim=c(0,2))
abline(h=data$F_t[length(data$F_t)], col="red", lwd=3)