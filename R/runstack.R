runstack <- function(savedir, iter, nodes, param, mean, cov, modname, rewrite){
	iterpath <- file.path(savedir, itervec[iter])
	dir.create(iterpath, showWarnings=FALSE)
	
	## true values
	set.seed(iter)
	vbk_choose <- ifelse("K" %in% param, exp(rnorm(1, mean=mean["K"], sd=sqrt(cov["K","K"]))), exp(mean["K"]))
	M_choose <- ifelse("M" %in% param, exp(rnorm(1, mean=mean["M"], sd=sqrt(cov["M","M"]))), exp(mean["M"]))
	Linf_choose <- ifelse("Loo" %in% param, exp(rnorm(1, mean=mean["Loo"], sd=sqrt(cov["Loo","Loo"]))), exp(mean["Loo"]))
	plist <- create_lh_list(linf=Linf_choose, vbk=vbk_choose, t0=-1.77,
							lwa=0.0076475, lwb=2.96,
							M=M_choose,
							M50=16.9, maturity_input="length",
							S50=16.9, S95=42.6, selex_input="length",
							SigmaF=0.2, SigmaR=0.737)

	if(rewrite==TRUE | file.exists(file.path(iterpath, "True.rds"))==FALSE){
		## generate data
		data <- generate_data(modpath=savedir, itervec=iter, 
						Fdynamics="Endogenous", Rdynamics="AR", 
						lh=plist, 
						Nyears=20, Nyears_comp=20, comp_sample=200,
						init_depl=c(0.10,0.90), 
						seed=rep(iter+1000,iter),
						nburn=50,
						rewrite=TRUE)
	}

		data <- readRDS(file.path(iterpath, "True.rds"))
		data_list <- list("years"=data$years, "LF"=data$LF)	
		# plot_LCfits(Inputs=list("LF"=data_list$LF))

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
									SigmaF=SigmaF, SigmaR=SigmaR))		
		 		out <- run_LIME(modpath=NULL, lh=lhinp, input_data=data_list, est_sigma="log_sigma_R", data_avail="LC", rewrite=TRUE)	
			return(out)
		})
		saveRDS(res, file.path(iterpath, paste0("res_", modname, ".rds")))
	}
		return(paste0("Ran iter ", iter, " in ", savedir))

}