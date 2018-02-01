stacking_figure <- function(res_dir, itervec, seed){
	set.seed(seed)
par(mfrow=c(3,3), mar=c(0,5,0,0), omi=c(0.75,0.5,0.5,0.5))
col <- rev(brewer.pal(4, "Set1"))
plot_iter <- sample(itervec,3)
for(i in 1:length(plot_iter)){
	choose_iter=plot_iter[i]
	choose_dir <- file.path(res_dir, choose_iter)
	res <- readRDS(file.path(choose_dir, "res_MK.rds"))
	true <- readRDS(file.path(choose_dir, "True.rds"))
	mres <- readRDS(file.path(choose_dir, "res_FishLifeMeans.rds"))
	ires <- readRDS(file.path(choose_dir, "res_IterTrue.rds"))


	plot(x=exp(nodesMK[,"K"]), y=exp(nodesMK[,"M"]), pch=16, xlim=c(0.9*min(exp(nodesMK[,"K"])), 1.1*max(exp(nodesMK[,"K"]))), ylim=c(0.9*min(exp(nodesMK[,"M"])), 1.1*max(exp(nodesMK[,"M"]))), ylab="", cex.lab=2, xlab="", xaxt="n", las=2, xaxs="i", yaxs="i", cex=1.5, cex.axis=1.5, col="gray")
	points(true$vbk, true$M, col=col[1], pch=19, cex=2.5)
	points(exp(MeanMK["K"]), exp(MeanMK["M"]), col=col[3], pch=19, cex=2.5)
	if(i==length(plot_iter)){
		axis(1, cex.axis=1.5)
		mtext(side=1, "von Bertalanffy k", line=3, cex=1.3)
	}
	if(i==2){
		mtext(side=2, "Natural mortality", line=3, cex=1.3)
	}

	## RELATIVE BIOMASS
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
	plot(x=1, y=1, type="n", xlim=c(1,20), ylim=c(0,1.5), ylab="", cex.lab=2, xlab="", xaxt="n", las=2, xaxs="i", yaxs="i", cex.axis=1.5)
	# if(i==1){
	# 	mtext(side=3, "Relative spawning biomass", line=1, cex=1.5)
	# }
	if(i==length(plot_iter)){
		axis(1, cex.axis=1.5)
		mtext(side=1, "Year", line=3, cex=1.3)
	}
	if(i==2){
		mtext(side=2, "Relative spawning biomass", line=3, cex=1.3)
	}
	## estimates for each node
	for(x in 1:length(res)){
		if(is.null(res[[x]]$df) == FALSE) sub <- res[[x]]$Report
		lines(x=1:length(sub$D_t), y=sub$D_t, lwd=3, col="#AAAAAA80")
	}	

	# polygon(x=c(1:ncol(quant), ncol(quant):1), y=c(quant[1,], rev(quant[3,])), col="#55555550", border=NA)
	# lines(x=1:ncol(quant), y=quant[2,], lwd=3, col="black")
	lines(x=1:length(ires$Report$D_t), y=ires$Report$D_t, col=col[2], lwd=3)
	lines(x=1:length(mres$Report$D_t), y=mres$Report$D_t, col=col[3], lwd=3)
	lines(x=1:length(true$D_t), y=true$D_t, col=col[1], lwd=3)

	## predictive stacking
	lines(x=1:length(stack), y=stack, col=col[4], lwd=3)


	plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0,2), ylim=c(0,8.5), ylab="", cex.lab=2, xaxt="n", las=2, cex.axis=1.5)
	if(i==length(plot_iter)){
		axis(1, cex.axis=1.5)
		mtext(side=1, "Terminal year estimate", line=3, cex=1.3)
	}
	if(i==2){
		mtext(side=2, "Density", line=3, cex=1.3)
	}

	xx <- seq(0,2,by=0.001)
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
	polygon(x=c(xx, rev(xx)), y=c(ysum, rep(0, length(xx))), col=paste0(col[4],"99"), border=col[4])
	points(x=true$D_t[length(true$D_t)], y=0, xpd=NA, col=paste0(col[1],"99"), pch=16, cex=3)
	points(x=ires$Report$D_t[length(stack)], y=0, xpd=NA, col=paste0(col[2],"99"), pch=16, cex=3)
	points(x=mres$Report$D_t[length(stack)], y=0, xpd=NA, col=paste0(col[3],"99"), pch=16, cex=3)


	# ### RECRUITMENT
	# find <- lapply(1:length(res), function(x){	
	# 	if(is.null(res[[x]]$df)==FALSE){
	# 		sub <- res[[x]]$Report
	# 		mle <- sub$R_t		

	# 		sdrep <- summary(res[[x]]$Sdreport)
	# 		sd <- sdrep[which(rownames(sdrep)=="lR_t"),2]
	# 	}
	# 	if(is.null(res[[x]]$df)){
	# 		mle <- rep(NA, 20)
	# 		sd <- rep(NA, 20)
	# 	}		
	# 		df <- data.frame("Node"=x, "Year"=1:length(mle), "MLE"=mle, "SE"=sd)
	# 		return(df)
	# })
	# find <- do.call(rbind, find)
	# yrs <- unique(find$Year)
	# stack <- sapply(yrs, function(x){
	# 	yrsub <- find %>% filter(Year==yrs[x])
	# 	stack <- sum(yrsub[,"MLE"] * weightsMK, na.rm=TRUE)
	# })

	# plot(x=1, y=1, type="n", xlim=c(1,20), ylim=c(0,2.75), ylab="Estimate", cex.lab=2, xlab="", xaxt="n", las=2, xaxs="i", yaxs="i", cex.axis=1.5)
	# if(i==1){
	# 	mtext(side=3, "Recruitment", line=1, cex=1.5)
	# }
	# if(i==length(plot_iter)){
	# 	axis(1, cex.axis=1.5)
	# 	mtext(side=1, "Year", line=3, cex=2)
	# }
	# ## estimates for each node
	# for(x in 1:length(res)){
	# 	if(is.null(res[[x]]$df) == FALSE) sub <- res[[x]]$Report
	# 	lines(x=1:length(sub$R_t), y=sub$R_t, lwd=3, col="#AAAAAA80")
	# }	

	# # polygon(x=c(1:ncol(quant), ncol(quant):1), y=c(quant[1,], rev(quant[3,])), col="#55555550", border=NA)
	# # lines(x=1:ncol(quant), y=quant[2,], lwd=3, col="black")
	# lines(x=1:length(ires$Report$R_t), y=ires$Report$R_t, col=col[2], lwd=3)
	# lines(x=1:length(mres$Report$R_t), y=mres$Report$R_t, col=col[3], lwd=3)
	# lines(x=1:length(true$R_t), y=true$R_t, col=col[1], lwd=3)

	# ## predictive stacking
	# lines(x=1:length(stack), y=stack, col=col[4], lwd=3)


	# plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=c(0,2.5), ylim=c(0,3.5), ylab="Density", cex.lab=2, xaxt="n", yaxt="n", las=2, cex.axis=1.5)
	# if(i==length(plot_iter)){
	# 	axis(1, cex.axis=1.5)
	# 	mtext(side=1, "Terminal year estimate", line=3, cex=2)
	# }
	# xx <- seq(0,2.5,by=0.001)
	# axis(2, las=2, cex.axis=1.5)
	# stack_dens <- sapply(1:length(res), function(x){
	# 	sub <- find %>% filter(Node==x)
	# 	term <- sub %>% filter(Year==max(Year))
	# 	if(is.na(term[,"MLE"])==FALSE & is.na(term[,"SE"])==FALSE){
	# 		dxx <- dlnorm(xx, mean=log(term[,"MLE"]), sd=term[,"SE"])
	# 		polygon(x=c(xx, rev(xx)), y=c(dxx, rep(0, length(dxx))), col="#AAAAAA80", border="#222222")
	# 		# xx <- density(rlnorm(10000, mean=log(term[,"MLE"]), sd=term[,"SE"]))
	# 		# polygon(x=c(xx$x, rev(xx$x)), y=c(xx$y, rep(0,length(xx$x))), col="#AAAAAA80", border="#222222")
	# 	} else { dxx <- rep(NA, length(xx))}
	# 	return(dxx)
	# })
	# wdens <- sapply(1:length(res), function(x){
	# 	return(weightsMK[x] * stack_dens[,x])
	# })
	# ysum <- sapply(1:nrow(wdens), function(x){
	# 	sub <- wdens[x,]
	# 	return(sum(sub, na.rm=TRUE))
	# })
	# polygon(x=c(xx, rev(xx)), y=c(ysum, rep(0, length(xx))), col=paste0(col[4],"99"), border=col[4])
	# points(x=true$R_t[length(true$R_t)], y=0, xpd=NA, col=paste0(col[1],"99"), pch=16, cex=3)
	# points(x=ires$Report$R_t[length(stack)], y=0, xpd=NA, col=paste0(col[2],"99"), pch=16, cex=3)
	# points(x=mres$Report$R_t[length(stack)], y=0, xpd=NA, col=paste0(col[3],"99"), pch=16, cex=3)
	if(i==1) legend("topright", legend=c("True", "Run at truth", "FishLife means", "Quadrature nodes", "Predictive stacking"), col=c(col[1], col[2], col[3], "#AAAAAA80", col[4]), lwd=4, cex=2)

}
}