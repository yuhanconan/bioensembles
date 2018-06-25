stackdensity <- function(savedir, modname, model, iter, vals, nodes, param, mean, cov, weights=NULL, rewrite=FALSE, rewrite_summary=TRUE){

	Integrate = function( Nodes_i, h_i, Mean, Cov, n_interp=100 ){

		if(ncol(Nodes_i)>1){
	  		# Interpolate PDF and get integral
			logdens_i = dmvnorm( Nodes_i, mean=Mean, sigma=Cov, log=TRUE )

			interp_logdens = akima::interp( x=Nodes_i[,1], y=Nodes_i[,2], z=logdens_i, nx=n_interp, ny=n_interp, duplicate= "median")
	    	cell_area = mean(diff(interp_logdens$x)) * mean(diff(interp_logdens$y))
	     	sumdens = sum( exp(interp_logdens$z), na.rm=TRUE ) * cell_area	    
 
    	  	# Interpolate product of function and PDF
    	  	interp_h = akima::interp( x=Nodes_i[,1], y=Nodes_i[,2], z=h_i, nx=n_interp, ny=n_interp, duplicate="median" )
    	  	interp_z = exp(interp_logdens$z) * interp_h$z / sumdens	
		}

		if(ncol(Nodes_i)==1){
			logdens_i = dmvnorm( Nodes_i, mean=Mean, sigma=Cov, log=TRUE )
			interp_logdens <- approx(x=Nodes_i[,1], y=logdens_i )
	    	cell_area = mean(diff(interp_logdens$x))
	    	sumdens = sum( exp(interp_logdens$y), na.rm=TRUE) * cell_area

	    	interp_h <- approx( x=Nodes_i[,1], y=h_i )
	    	interp_z <- exp(interp_logdens$y) * interp_h$y / sumdens
		}

	      # Calculate interal and return
	      Integral = sum( interp_z, na.rm=TRUE ) * cell_area
	  return( Integral )
	}

	# est_summary <- list()
	# for(y in 1:length(iters)){
		# print(iters[y])

		if(is.null(iter)==FALSE) iterpath <- file.path(savedir, iter)
		if(is.null(iter)) iterpath <- savedir

		pres <- readRDS(file.path(iterpath, paste0(modname, "_res_stacking_", model, ".rds")))
		if(file.exists(file.path(iterpath, paste0(modname, "_res_Means_", model,".rds")))){
			mres <- readRDS(file.path(iterpath, paste0(modname, "_res_Means_", model, ".rds")))
		} else { mres <- NULL}
		if(all(is.null(iter))==FALSE){
			if(file.exists(file.path(iterpath, paste0("res_IterTrue_", model, ".rds")))){
				ires <- readRDS(file.path(iterpath, paste0("res_IterTrue_", model, ".rds")))
			} else { ires <- NULL }
			true <- readRDS(file.path(iterpath, "True.rds"))
		}

		fnodes <- sapply(1:length(pres), function(x) pres[[x]]@SPR)
		sum( fnodes * weights )






		if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("results_summary_", modname,"_", model, ".rds")))==FALSE){

		node_results <- list()
		for(i in 1:length(vals)){
			if(vals[i]=="Depletion"){
				if(model=="LIME"){
					SB_nodes <- lapply(1:length(pres), function(x){
						if(is.null(pres[[x]]$df)==FALSE){
							sub <- pres[[x]]$Report
							mle <- sub$D_t
							years <- 1:length(mle)
							sdrep <- summary(pres[[x]]$Sdreport)
							sd <- sdrep[which(rownames(sdrep)=="lD_t"),2]

							up <- mle + 0.674*sd
							low <- mle - 0.674*sd
							df <- data.frame("Node"=x, "Year"=years, "MLE"=mle, "SE"=sd, "Value"="Depletion")

						}
						if(is.null(pres[[x]]$df)){
							df <- data.frame("Node"=NULL, "Year"=NULL, "MLE"=NULL, "SE"=NULL, "Value"=NULL)
						}
						return(df)
					})
					SB_nodes <- do.call(rbind, SB_nodes)
					node_results[[i]] <- SB_nodes
				}
			}
			if(vals[i]=="F"){

				if(model=="LIME"){
					F_nodes <- lapply(1:length(pres), function(x){
						if(is.null(pres[[x]]$df)==FALSE){
							sub <- pres[[x]]$Report
							mle <- sub$F_ft[1,]
							years <- 1:length(mle)
							sdrep <- summary(pres[[x]]$Sdreport)
							sd <- sdrep[which(rownames(sdrep)=="lF_ft"),2]

							up <- mle + 0.674*sd
							low <- mle - 0.674*sd
							df <- data.frame("Node"=x, "Year"=years, "MLE"=mle, "SE"=sd, "Value"="F")

						}
						if(is.null(pres[[x]]$df)){
							df <- data.frame("Node"=NULL, "Year"=NULL, "MLE"=NULL, "SE"=NULL, "Value"=NULL)
						}
						return(df)
					})
					F_nodes <- do.call(rbind, F_nodes)
					node_results[[i]] <- F_nodes
				}			
			}
			if(vals[i]=="SPR"){
				SPR_nodes <- lapply(1:length(pres), function(x){
					if(model=="LIME"){
						if(is.null(pres[[x]]$df)==FALSE){
							sub <- pres[[x]]$Report
							mle <- sub$SPR_t
							years <- 1:length(mle)
							sdrep <- summary(pres[[x]]$Sdreport)
							sd <- sdrep[which(rownames(sdrep)=="SPR_t"),2]

							up <- mle + 0.674*sd
							low <- mle - 0.674*sd
							df <- data.frame("Node"=x, "Year"=years, "MLE"=mle, "SE"=sd, "Value"="SPR")
						}
						if(is.null(pres[[x]]$df)){
							df <- data.frame("Node"=NULL, "Year"=NULL, "MLE"=NULL, "SE"=NULL, "Value"=NULL)
						}
					}

					if(model=="LBSPR"){
						mle <- pres[[x]]@SPR
						sd <- sqrt(pres[[x]]@Vars[,"SPR"])
						years <- pres[[x]]@Years

						up <- mle + 0.674*sd
						low <- mle - 0.674*sd
						df <- data.frame("Node"=x, "Year"=years, "MLE"=mle, "SE"=sd, "Value"="SPR")
					}
					return(df)
				})
				SPR_nodes <- do.call(rbind, SPR_nodes)		
				node_results[[i]] <- SPR_nodes	
			}
		}
		names(node_results) <- vals
		saveRDS(node_results, file.path(iterpath, paste0("results_nodes_", modname, "_", model, ".rds")))

		## density at each node
		density_results <- list()
		for(i in 1:length(vals)){
			if(vals[i]=="Depletion"){
				yrs <- unique(node_results[["Depletion"]]$Year)
				SB_seq <- seq(0,4,by=0.001)
				if(max(node_results[["Depletion"]][,"MLE"]) > max(SB_seq)){
					SB_seq <- seq(0, max(SB_seq)+1, by=0.001)
				}

				SB_dens <- lapply(1:length(yrs), function(y){
					findYear <- SB_nodes %>% filter(Year==yrs[y])
					byYear <- sapply(1:length(pres), function(x){
						sub <- findYear %>% filter(Node==x)
						if(nrow(sub)==0) return(rep(NA, length(SB_seq)))
						if(nrow(sub)>0){
							## density at each node
							dxx <- dlnorm(SB_seq, mean=log(sub[,"MLE"])-sub[,"SE"]^2/2, sd=sub[,"SE"])
							return(dxx)							
						}
					})
					return(byYear)
					rownames(byYear) <- SB_seq
				})
				density_results[[i]] <- SB_dens
			}
			if(vals[i]=="F"){
				yrs <- unique(node_results[["F"]]$Year)
				F_seq <- seq(0.001,2,by=0.001)
				if(max(node_results[["F"]][,"MLE"]) > max(F_seq)){
					F_seq <- seq(0.001,max(F_seq)+1,by=0.001)
				}
				
				F_dens <- lapply(1:length(yrs), function(y){
					findYear <- F_nodes %>% filter(Year==yrs[y])
					byYear <- sapply(1:length(pres), function(x){
						sub <- findYear %>% filter(Node==x)
						if(nrow(sub)==0) return(rep(NA, length(F_seq)))
						if(nrow(sub)>0){
							## density at each node
							dxx <- dlnorm(F_seq, mean=log(sub[,"MLE"])-sub[,"SE"]^2/2, sd=sub[,"SE"])
							return(dxx)
						}
					})
					return(byYear)
					rownames(byYear) <- F_seq
				})
				density_results[[i]] <- F_dens
			}
			if(vals[i]=="SPR"){
				yrs <- unique(node_results[["SPR"]]$Year)
				SPR_seq <- seq(0,0.999,by=0.001)
				# if(max(node_results[["SPR"]][,"MLE"],  rm.na=TRUE) > max(SPR_seq)){
				# 	SPR_seq <- seq(0, max(SPR_seq)+1, by=0.001)
				# }

				SPR_dens <- lapply(1:length(yrs), function(y){
					findYear <- SPR_nodes %>% filter(Year==yrs[y])
					byYear <- sapply(1:length(pres), function(x){
						sub <- findYear %>% filter(Node==x)
						if(nrow(sub)==0) return(rep(NA, length(SPR_seq)))
						if(nrow(sub)>0){
							## density at each node
							if(all(is.na(sub[,"SE"]))) sdinp <- 0.001
							if(all(is.na(sub[,"SE"]))==FALSE) sdinp <- sub[,"SE"]
							dxx <- dlnorm(SPR_seq, mean=log(sub[,"MLE"])-sdinp^2/2, sd=sdinp)
							return(dxx)						
						}
					})
					rownames(byYear) <- SPR_seq
					return(byYear)
				})
				density_results[[i]] <- SPR_dens
			}
		}
		names(density_results) <- vals
		# sum(density_results[[2]][[1]][,1])
		saveRDS(density_results, file.path(iterpath, paste0("results_density_", modname, "_", model, ".rds")))


		stack_interp <- list()
		stack_interp_dens <- list()
		run_nodes <- unique(node_results[[1]][,"Node"])
		unique_nodes <- matrix(unique(nodes[run_nodes,]), ncol=ncol(nodes))
		colnames(unique_nodes) <- colnames(nodes)
		Match <- match(unique_nodes[,1], nodes[,1])
		for(i in 1:length(vals)){
			if(vals[i]=="Depletion"){
				unique_dens <- lapply(1:length(density_results[["Depletion"]]), function(x){
					density_results[["Depletion"]][[x]][,Match]
				})
				SB_ysum2 <- lapply(1:length(unique_dens), function(x){
					byyr <- unique_dens[[x]]
					out <- sapply(1:nrow(byyr), function(z){
						subyr <- byyr[z,]
						bad <- which(is.finite(subyr)==FALSE)
						if(length(bad) > 0){
							newyr <- subyr[-bad]
							newnodes <- matrix(unique_nodes[-bad,], ncol=ncol(nodes))
						}
						if(length(bad) == 0){
							newyr <- subyr
							newnodes <- unique_nodes
						}
						newcov <- matrix(cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)], nrow=ncol(nodes), ncol=ncol(nodes))
						colnames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]
						rownames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]
						
						newmean <- mean[which(names(mean) %in% param)]
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Loo")) newmean[which(names(mean)=="Loo")] <-  exp(mean["Lm"])/exp(mean["Loo"])
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Lm")) newmean[which(names(mean)=="Lm")] <-  exp(mean["Lm"])/exp(mean["Loo"])

						Integrate(h_i=newyr,
							Nodes_i=newnodes,
							Mean=newmean, 
							Cov=newcov)
					})
					return(out)
				})
				SB_stack2 <- sapply(1:length(SB_ysum2), function(x) weighted.mean(SB_seq, w=SB_ysum2[[x]]))
				stack_interp[[i]] <- SB_stack2
				stack_interp_dens[[i]] <- SB_ysum2
			}
			if(vals[i]=="F"){
				unique_dens <- lapply(1:length(density_results[["F"]]), function(x){
					density_results[["F"]][[x]][,Match]
				})
				F_ysum2 <- lapply(1:length(unique_dens), function(x){
					byyr <- unique_dens[[x]]
					out <- sapply(1:nrow(byyr), function(z){
						subyr <- byyr[z,]
						bad <- which(is.finite(subyr)==FALSE)
						if(length(bad) > 0){
							newyr <- subyr[-bad]
							newnodes <- matrix(unique_nodes[-bad,], ncol=ncol(nodes))
						}
						if(length(bad) == 0){
							newyr <- subyr
							newnodes <- unique_nodes
						}
						newcov <- matrix(cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)], nrow=ncol(nodes), ncol=ncol(nodes))
						colnames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]
						rownames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]
						
						newmean <- mean[which(names(mean) %in% param)]
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Loo")) newmean[which(names(mean)=="Loo")] <-  exp(mean["Lm"])/exp(mean["Loo"])
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Lm")) newmean[which(names(mean)=="Lm")] <-  exp(mean["Lm"])/exp(mean["Loo"])

						Integrate(h_i=newyr,
							Nodes_i=newnodes,
							Mean=newmean, 
							Cov=newcov)
					})
					return(out)
				})
				F_stack2 <- sapply(1:length(F_ysum2), function(x) weighted.mean(F_seq, w=F_ysum2[[x]]))
				stack_interp[[i]] <- F_stack2
				stack_interp_dens[[i]] <- F_ysum2
			}
			if(vals[i]=="SPR"){
				unique_dens <- lapply(1:length(density_results[["SPR"]]), function(x){
					density_results[["SPR"]][[x]][,Match]
				})
				SPR_ysum2 <- lapply(1:length(unique_dens), function(x){
					byyr <- unique_dens[[x]]
					out <- sapply(1:nrow(byyr), function(z){
						subyr <- byyr[z,]
						bad <- which(is.finite(subyr)==FALSE)
						if(length(bad) > 0){
							newyr <- subyr[-bad]
							newnodes <- matrix(unique_nodes[-bad,], ncol=ncol(nodes))
						}
						if(length(bad) == 0){
							newyr <- subyr
							newnodes <- unique_nodes
						}
						newcov <- matrix(cov[which(rownames(cov) %in% param), which(colnames(cov) %in% param)], nrow=ncol(nodes), ncol=ncol(nodes))
						colnames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]
						rownames(newcov) <- colnames(cov)[which(colnames(cov) %in% param)]

						newmean <- mean[which(names(mean) %in% param)]
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Loo")) newmean[which(names(mean)=="Loo")] <-  exp(mean["Lm"])/exp(mean["Loo"])
						if(any(grepl("LmLoo", colnames(newnodes))) & any(param=="Lm")) newmean[which(names(mean)=="Lm")] <-  exp(mean["Lm"])/exp(mean["Loo"])

						Integrate(h_i=newyr,
							Nodes_i=newnodes,
							Mean=newmean, 
							Cov=newcov)
					})
					return(out)
				})
				SPR_stack2 <- sapply(1:length(SPR_ysum2), function(x) weighted.mean(SPR_seq, w=SPR_ysum2[[x]], na.rm=TRUE))
				stack_interp[[i]] <- SPR_stack2
				stack_interp_dens[[i]] <- SPR_ysum2
			}
	}
		names(stack_interp) <- vals
		names(stack_interp_dens) <- vals
		saveRDS(stack_interp, file.path(iterpath, paste0("results_stackinterpolation_", modname, "_", model, ".rds")))
		saveRDS(stack_interp_dens, file.path(iterpath, paste0("results_stackinterpolation_dens_", modname, "_", model, ".rds")))



	if(rewrite==TRUE | rewrite_summary==TRUE | file.exists(file.path(iterpath, paste0("results_summary_", modname, "_", model, ".rds")))==FALSE & is.null(iter)==FALSE){

		stack_interp <- readRDS(file.path(iterpath, paste0("results_stackinterpolation_",modname, "_", model, ".rds")))
		stack_interp_dens <- readRDS(file.path(iterpath, paste0("results_stackinterpolation_dens_", modname, "_", model, ".rds")))
		est <- lapply(1:length(vals), function(x){
			yrs <- 1:length(stack_interp[[x]])
				if(vals[[x]]=="Depletion"){
					if(model=="LIME"){
						if(all(is.null(iter))==FALSE){
							tval <- true$D_t

							ival <- ires$Report$D_t
							if(is.null(ival)==FALSE){
								sdrep <- summary(ires$Sdreport)
								isd <- sdrep[which(rownames(sdrep)=="lD_t"),2]								
							} else{
								isd <- NULL
							}
						} else{
							tval <- NULL
							ival <- NULL
							isd <- NULL
						}

						if(all(is.null(mres))==FALSE){
							mval <- mres$Report$D_t
							sdrep <- summary(mres$Sdreport)
							msd <- sdrep[which(rownames(sdrep)=="lD_t"),2]
						} else{
							mval <- NULL
							msd <- NULL
						}

						sval <- stack_interp[[x]]
						ssd <- sapply(1:length(stack_interp_dens[[x]]), function(y) sqrt(wtd.var(SB_seq, w=stack_interp_dens[[x]][[y]])))		
					}
				}
				if(vals[[x]]=="F"){
					if(model=="LIME"){
						if(all(is.null(iter))==FALSE){
							tval <- true$F_y

							ival <- ires$Report$F_y
							if(is.null(ival)==FALSE){
								sdrep <- summary(ires$Sdreport)
								isd <- sdrep[which(rownames(sdrep)=="lF_y"),2]								
							} else{
								isd <- NULL
							}
						} else{
							tval <- NULL
							ival <- NULL
							isd <- NULL
						}

						if(all(is.null(mres))==FALSE){
							mval <- mres$Report$F_y
							sdrep <- summary(mres$Sdreport)
							msd <- sdrep[which(rownames(sdrep)=="lF_y"),2]							
						} else{
							mval <- NULL
							msd <- NULL
						}


						sval <- stack_interp[[x]]
						ssd <- sapply(1:length(stack_interp_dens[[x]]), function(y) sqrt(wtd.var(F_seq, w=stack_interp_dens[[x]][[y]])))		
					}
				}


				if(vals[[x]]=="SPR"){
					if(model=="LBSPR"){
						if(all(is.null(iter))==FALSE){
							tval <- true$SPR_t

							ival <- ires@SPR
							isd <- sqrt(ires@Vars[,"SPR"])	
						} else {
							tval <- NULL
							ival <- NULL
							isd <- NULL
						}
						
						mval <- mres@SPR
						msd <- sqrt(mres@Vars[,"SPR"])		

						sval <- stack_interp[[x]]
						ssd <- sapply(1:length(stack_interp_dens[[x]]), function(y) sqrt(wtd.var(SPR_seq, w=stack_interp_dens[[x]][[y]])))		
					}
					if(model=="LIME"){
						if(all(is.null(iter))==FALSE){
							tval <- true$SPR_t

							ival <- ires$Report$SPR_t
							if(is.null(ival)==FALSE){
								sdrep <- summary(ires$Sdreport)
								isd <- sdrep[which(rownames(sdrep)=="SPR_t"),2]								
							} else{
								isd <- NULL
							}
						} else{
							tval <- NULL
							ival <- NULL
							isd <- NULL
						}

						if(all(is.null(mres))==FALSE){
							mval <- mres$Report$SPR_t
							sdrep <- summary(mres$Sdreport)
							msd <- sdrep[which(rownames(sdrep)=="SPR_t"),2]							
						} else{
							mval <- NULL
							msd <- NULL
						}

						sval <- stack_interp[[x]]
						ssd <- sapply(1:length(stack_interp_dens[[x]]), function(y) sqrt(wtd.var(SPR_seq, w=stack_interp_dens[[x]][[y]])))		
					}
				}

				if(all(is.null(iter))==FALSE){
					iup <- ival + 0.674*isd
					ilow <- ival - 0.674*isd
					icover <- ifelse(tval <= iup & tval >= ilow, 1, 0)
				} else {
					iup <- NULL
					ilow <- NULL
					icover <- NULL
				}

				mup <- mval + 0.674*msd
				mlow <- mval - 0.674*msd
				if(all(is.null(iter))==FALSE){
					mcover <- ifelse(tval <= mup & tval >= mlow, 1, 0)
				} else { mcover <- NULL}

				sup <- sval + 0.674*ssd
				slow <- sval - 0.674*ssd
				if(all(is.null(iter))==FALSE) {
					scover <- ifelse(tval <= sup & tval >= slow, 1, 0)
				} else { scover <- NULL }

				out <- data.frame("Year"=yrs, 
					"Variable"=vals[x], 
					"Type"=c(rep("BestCase",length(ival)), rep("Means",length(mval)), rep("Stacking",length(sval))), 
					"Value"=c(ival, mval, sval), 
					"SD"=c(isd, msd, ssd))
				if(all(is.null(iter))==FALSE){
					out <- cbind.data.frame(out, 
						"Cover50"=c(icover, mcover, scover), 
						"True"=tval)
					out <- out %>% mutate("RE"= (Value - True)/True)	
				} 
			return(out)
		})
		names(est) <- vals
		saveRDS(est, file.path(iterpath, paste0("results_summary_", modname, "_", model, ".rds")))
		return(est)
	}
}

}
