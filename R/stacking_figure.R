stacking_figure <- function(savedir, modname, model, nodes, mean, iters, value="Depletion"){

ex <- readRDS(file.path(savedir, 1, paste0("True.rds")))
Nyears <- ex$Nyear

if(Nyears>1) par(mfrow=c(length(iters),3), mar=c(0,5,0,0), omi=c(0.75,0.5,0.5,0.5))
if(Nyears==1) par(mfrow=c(length(iters),2), mar=c(0,5,0,0), omi=c(0.75,0.5,0.5,0.5))
col <- rev(brewer.pal(4, "Set1"))

for(i in 1:length(iters)){
	summary <- readRDS(file.path(savedir, iters[i], paste0("results_summary_", modname,"_",model,".rds")))
	summary <- do.call(rbind, summary)
	noderes <- readRDS(file.path(savedir, iters[i], paste0("results_nodes_", modname, "_", model, ".rds")))
	density <- readRDS(file.path(savedir, iters[i], paste0("results_density_", modname, "_", model, ".rds")))
	stackinterp <- readRDS(file.path(savedir, iters[i], paste0("results_stackinterpolation_dens_", modname, "_", model, ".rds")))
	true <- readRDS(file.path(savedir, iters[i], paste0("True.rds")))
	ires <- readRDS(file.path(savedir, iters[i], paste0("res_IterTrue_", model, ".rds")))
	mres <- readRDS(file.path(savedir, iters[i], paste0(modname, "_res_Means_", model,".rds")))

	## true value vs. FishLife means on nodes
	plot(x=exp(nodes[,"K"]), y=exp(nodes[,"M"]), pch=16, xlim=c(0.9*min(exp(nodes[,"K"])), 1.1*max(exp(nodes[,"K"]))), ylim=c(0.9*min(exp(nodes[,"M"])), 1.1*max(exp(nodes[,"M"]))), ylab="", cex.lab=2, xlab="", xaxt="n", las=2, xaxs="i", yaxs="i", cex=1.5, cex.axis=1.5, col="gray")
	points(true$vbk, true$M, col=col[1], pch=19, cex=2.5)
	points(exp(mean["K"]), exp(mean["M"]), col=col[3], pch=19, cex=2.5)
	if(i==length(iters)){
		axis(1, cex.axis=1.5)
		mtext(side=1, "von Bertalanffy k", line=3, cex=1.3)
	}
	if(i==which(iters==median(iters))) mtext(side=2, "Natural mortality", line=4, cex=1.3)

	if(value=="Depletion") ylim <- c(0,2)
	if(value=="SPR") ylim <- c(0,1)
	if(value=="F") ylim <- c(0,1)

  if(Nyears>1){
	plot(x=1, y=1, type="n", xlim=c(1,length(true$SPR_t)), ylim=ylim, ylab="", cex.lab=2, xlab="", xaxt="n", las=2, xaxs="i", yaxs="i", cex.axis=1.5)
	# if(i==1){
	# 	mtext(side=3, "Relative spawning biomass", line=1, cex=1.5)
	# }
	if(i==length(iters)){
		axis(1, cex.axis=1.5)
		mtext(side=1, "Year", line=3, cex=1.3)
	}
		if(value=="Depletion")	if(i==which(iters==median(iters))) mtext(side=2, "Relative spawning biomass", line=3, cex=1.3)
		if(value=="F") 	if(i==which(iters==median(iters))) mtext(side=2, "Fishing mortality", line=3, cex=1.3)
		if(value=="SPR") 	if(i==which(iters==median(iters))) mtext(side=2, "SPR", line=3, cex=1.3)
	## estimates for each node
	nodevec <- 1:nrow(nodes)
	for(n in 1:length(nodevec)){
		sub <- noderes[[value]] %>% filter(Node == nodevec[n])
		lines(x=sub$Year, y=sub$MLE, col="#AAAAAA80", lwd=3)
	}
	if(value=="Depletion"){
		lines(x=1:length(ires$Report$D_t), y=ires$Report$D_t, col=col[2], lwd=3)
		lines(x=1:length(mres$Report$D_t), y=mres$Report$D_t, col=col[3], lwd=3)
		lines(x=1:length(true$D_t), y=true$D_t, col=col[1], lwd=3)
	}
	if(value=="SPR"){
		if(model=="LBSPR"){
			lines(x=1:length(ires@SPR), y=ires@SPR, col=col[2], lwd=3)
			lines(x=1:length(mres@SPR), y=mres@SPR, col=col[3], lwd=3)
			lines(x=1:length(true$SPR_t), y=true$SPR_t, col=col[1], lwd=3)
		}
		if(model=="LIME"){
			lines(x=1:length(ires$Report$SPR_t), y=ires$Report$SPR_t, col=col[2], lwd=3)
			lines(x=1:length(mres$Report$SPR_t), y=mres$Report$SPR_t, col=col[3], lwd=3)
			lines(x=1:length(true$SPR_t), y=true$SPR_t, col=col[1], lwd=3)
		}
	}
	if(value=="F"){
		lines(x=1:length(ires$Report$F_t), y=ires$Report$F_t, col=col[2], lwd=3)
		lines(x=1:length(mres$Report$F_t), y=mres$Report$F_t, col=col[3], lwd=3)
		lines(x=1:length(true$F_t), y=true$F_t, col=col[1], lwd=3)		
	}

	## predictive stacking
	stack <- t(summary %>% filter(Variable==value) %>% filter(Type=="StackByInterpolation") %>% select(Value))
	lines(x=1:length(stack), y=stack, col=col[4], lwd=3)
  }

	## density
	dsub <- density[[value]][[length(density[[value]])]]
	dseq <- as.numeric(rownames(dsub))

	plot(x=1, y=1, type="n", xaxs="i", yaxs="i", xlim=ylim, ylim=c(0,min(c(20,max(dsub, na.rm=TRUE)*1.1))), xaxt="n", ylab="", xlab="", cex.lab=2, las=2, cex.axis=1.5)
	if(i==length(iters)){
		axis(1, cex.axis=1.5)
		if(Nyears>1) mtext(side=1, "Terminal year estimate", line=3, cex=1.3)
		if(Nyears==1) mtext(side=1, "Estimate", line=3, cex=1.3)
	}
	if(i==which(iters==median(iters))) mtext(side=2, "Density", line=3, cex=1.3)
	for(d in 1:ncol(dsub)){
			polygon(x=c(dseq, rev(dseq)), y=c(dsub[1:length(dseq),d], rep(0, length(dsub[1:length(dseq),d]))), col="#AAAAAA80", border="#222222")
	}


	if(value=="Depletion"){
		yrs <- 1:length(ires$Report$D_t)
		idens <- lapply(1:length(yrs), function(y){
			mle <- ires$Report$D_t[y]
			sd <- summary(ires$Sdreport)[which(rownames(summary(ires$Sdreport))=="lD_t"),2][y]
				## density at each node
				dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
				return(dxx)
		})
		mdens <- lapply(1:length(yrs), function(y){
			mle <- mres$Report$D_t[y]
			sd <- summary(mres$Sdreport)[which(rownames(summary(mres$Sdreport))=="lD_t"),2][y]
				## density at each node
				dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
				return(dxx)
		})
		polygon(x=c(dseq,rev(dseq)), y=c(idens[[length(idens)]], rep(0,length(dseq))), col=paste0(col[2], "99"), border=col[2])
		polygon(x=c(dseq,rev(dseq)), y=c(mdens[[length(mdens)]], rep(0,length(dseq))), col=paste0(col[3], "99"), border=col[3])

		stackdens <- stackinterp[[value]][[length(stackinterp[[value]])]]
		polygon(x=c(dseq, rev(dseq)), y=c(stackdens[1:length(dseq)], rep(0, length(dseq))), col=paste0(col[4],"99"), border=col[4])


		abline(v=true$D_t[length(true$D_t)], col=col[1], lwd=4)
		abline(v=mres$Report$D_t[length(mres$Report$D_t)], col=col[3], lwd=4)
		abline(v=ires$Report$D_t[length(ires$Report$D_t)], col=col[2], lwd=4)
		abline(v=stack[length(stack)], col=col[4], lwd=4)
	}
	if(value=="SPR"){
		if(model=="LIME"){
			yrs <- 1:length(ires$Report$SPR_t)
			idens <- lapply(1:length(yrs), function(y){
				mle <- ires$Report$SPR_t[y]
				sd <- summary(ires$Sdreport)[which(rownames(summary(ires$Sdreport))=="SPR_t"),2][y]
					## density at each node
					dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
					return(dxx)
			})
			mdens <- lapply(1:length(yrs), function(y){
				mle <- mres$Report$SPR_t[y]
				sd <- summary(mres$Sdreport)[which(rownames(summary(mres$Sdreport))=="SPR_t"),2][y]
					## density at each node
					dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
					return(dxx)
			})
			polygon(x=c(dseq,rev(dseq)), y=c(idens[[length(idens)]], rep(0,length(dseq))), col=paste0(col[2], "99"), border=col[2])
			polygon(x=c(dseq,rev(dseq)), y=c(mdens[[length(mdens)]], rep(0,length(dseq))), col=paste0(col[3], "99"), border=col[3])


			abline(v=true$SPR_t[length(true$SPR_t)], col=col[1], lwd=4)
			abline(v=mres$Report$SPR_t[length(mres$Report$SPR_t)], col=col[3], lwd=4)
			abline(v=ires$Report$SPR_t[length(ires$Report$SPR_t)], col=col[2], lwd=4)
		}
		if(model=="LBSPR"){
			idens <- lapply(1:Nyears, function(y){
				mle <- ires@SPR[y]
				sd <- sqrt(ires@Vars[,"SPR"])
					dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
					return(dxx)
			})
			mdens <- lapply(1:Nyears, function(y){
				mle <- mres@SPR[y]
				sd <- sqrt(mres@Vars[,"SPR"])
					dxx <- dlnorm(dseq, mean=log(mle)-sd^2/2, sd=sd)
					return(dxx)
			})			
			polygon(x=c(dseq,rev(dseq)), y=c(idens[[length(idens)]], rep(0,length(dseq))), col=paste0(col[2], "99"), border=col[2])
			polygon(x=c(dseq,rev(dseq)), y=c(mdens[[length(mdens)]], rep(0,length(dseq))), col=paste0(col[3], "99"), border=col[3])

			abline(v=true$SPR_t[length(true$SPR_t)], col=col[1], lwd=4)
			abline(v=mres@SPR[length(mres@SPR)], col=col[3], lwd=4)
			abline(v=ires@SPR[length(ires@SPR)], col=col[2], lwd=4)

		}

	}	
	## predictive stacking
	stack <- t(summary %>% filter(Variable==value) %>% filter(Type=="Stacking") %>% select(Value))
	stackdens <- stackinterp[[value]][[length(stackinterp[[value]])]]
	polygon(x=c(dseq, rev(dseq)), y=c(stackdens[1:length(dseq)], rep(0, length(dseq))), col=paste0(col[4],"99"), border=col[4])
	abline(v=stack[length(stack)], col=col[4], lwd=4)



}


}