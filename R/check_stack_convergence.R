check_stack_convergence <- function(savedir, itervec, modname, model, sim){

	## check for non-converged runs
	nc <- NULL
	not_run <- NULL
	for(i in 1:length(itervec)){
		choose_dir <- file.path(savedir, itervec[i])
		files <- list.files(choose_dir)
		if(any(grepl("pdHess", files)) | any(grepl("highgradient", files)) | any(grepl("modelNA", files))) nc <- c(nc, i)
		if(sim==TRUE) if(any(grepl(paste0(modname, "_res_stacking_", model, ".rds"), files))==FALSE | any(grepl(paste0(modname, "_res_FishLifeMeans_", model, ".rds"), files))==FALSE | any(grepl(paste0(modname, "_res_IterTrue_", model, ".rds"), files))==FALSE) not_run <- c(not_run, i)
		if(sim==FALSE) if(any(grepl(paste0(modname, "_res_stacking_", model, ".rds"), files))==FALSE | any(grepl(paste0(modname, "_res_FishLifeMeans_", model, ".rds"), files))==FALSE) not_run <- c(not_run, i)
	}
	out <- NULL
	out$nc <- nc
	out$not_run <- not_run	
	## while there are non-converged runs
	return(out)
}