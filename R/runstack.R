runstack <- function(savedir, iter, tmax, nodes, param, mean, cov, modname, input_data, Fscenario, rewrite){
	iterpath <- file.path(savedir, itervec[iter])
	dir.create(iterpath, showWarnings=FALSE)

	## generate data
		## true values
		set.seed(iter)
		vbk_choose <- ifelse("K" %in% param, rlnorm(1, mean=mean["K"], sd=sqrt(cov["K","K"])), exp(mean["K"]))
		M_choose <- ifelse("M" %in% param, rlnorm(1, mean=mean["M"], sd=sqrt(cov["M","M"])), exp(mean["M"]))
		Linf_choose <- ifelse("Loo" %in% param, rlnorm(1, mean=mean["Loo"], sd=sqrt(cov["Loo","Loo"])), exp(mean["Loo"]))
		if(Fscenario=="equil"){
			SigmaF_inp <- 0.01
			SigmaR_inp <- 0.01
			rho_inp <- 0
			Fdynamics_inp <- "Constant"
		}
		if(Fscenario=="harvestdyn"){
			SigmaF_inp <- 0.2
			SigmaR_inp <- 0.6
			rho_inp <- 0.4
			Fdynamics_inp <- "Endogenous"
		}
		plist <- create_lh_list(linf=Linf_choose, vbk=vbk_choose, t0=-1.77,
								lwa=0.0076475, lwb=2.96,
								M=M_choose,
								M50=16.9, maturity_input="length",
								S50=16.9, S95=25, selex_input="length",
								SigmaF=SigmaF_inp, SigmaR=SigmaR_inp, rho=rho_inp,
								AgeMax=tmax)

	if(all(input_data==FALSE)){
		if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
			data <- generate_data(modpath=savedir, itervec=iter, 
							Fdynamics=Fdynamics_inp, Rdynamics="AR", 
							lh=plist, 
							Nyears=20, Nyears_comp=20, comp_sample=200,
							init_depl=c(0.10,0.90), 
							seed=rep(iter+1000,iter),
							nburn=50,
							rewrite=TRUE)
		}
		data <- readRDS(file.path(iterpath, "True.rds"))
		input_data <- list("years"=data$years, "LF"=data$LF)	
	}

		# plot_LCfits(Inputs=list("LF"=input_data$LF))

	## run at means from FishLife
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_FishLifeMeans.rds")))==FALSE){
				vbk_inp <- exp(mean["K"])
				M_inp <- exp(mean["M"])
				linf_inp <- exp(mean["Loo"])
		 		lhinp <- with(plist, 
		 				 create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=M50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=SigmaF, SigmaR=SigmaR,
									AgeMax=AgeMax))		
		 		out <- run_LIME(modpath=NULL, lh=lhinp, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE)	
		saveRDS(out, file.path(iterpath, paste0("res_FishLifeMeans.rds")))
	}

	## run at true values
	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_IterTrue.rds")))==FALSE){	
		out <- run_LIME(modpath=NULL, lh=plist, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE)	
		saveRDS(out, file.path(iterpath, paste0("res_IterTrue.rds")))
	}


	if(rewrite==TRUE | file.exists(file.path(iterpath, paste0("res_", modname, ".rds")))==FALSE){
		res <- lapply(1:nrow(nodes), function(x){

				vbk_inp <- ifelse("K" %in% param, exp(nodes[x,"K"]), exp(mean["K"]))
				M_inp <- ifelse("M" %in% param, exp(nodes[x,"M"]), exp(mean["M"]))
				linf_inp <- ifelse("Loo" %in% param, exp(nodes[x,"Loo"]), exp(mean["Loo"]))
		 		lhinp <- with(plist, 
		 				 create_lh_list(linf=linf_inp, vbk=vbk_inp, t0=t0,
									lwa=lwa, lwb=lwb,
									M=M_inp,
									M50=M50, maturity_input="length",
									S50=SL50, S95=SL95, selex_input="length",
									SigmaF=SigmaF, SigmaR=SigmaR,
									AgeMax=AgeMax))		
		 		out <- run_LIME(modpath=NULL, lh=lhinp, input_data=input_data, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE)	
			return(out)
		})
		saveRDS(res, file.path(iterpath, paste0("res_", modname, ".rds")))
	}
		return(paste0("Ran iter ", iter, " in ", savedir))

}
