multilevelannotationstep2 <- function(list_number, outloc1) {
    
    # ,adduct_weights=NA,max.time.diff=NA,filter.by=c('M+H'),max_isp=100,numnodes=2,
    # MplusH.abundance.ratio.check=FALSE,mass_defect_window=0.01,mass_defect_mode='pos'){
    
    setwd(outloc1)
    
    load("step1_results.Rda")
    load("global_cor.Rda")
    
    library(enviPat)
    data(isotopes)
    
    unlink("allmatches_with_isotopes.txt")
    
    outloc <- outloc1
    if (is.na(max.rt.diff) == FALSE) {
        max_diff_rt <- max.rt.diff
    }
    
    if (list_number > length(chemids_split)) {
        #list_number <- length(chemids_split)
        stop(paste('Invalid list number. Must be less than or equal to ',length(chemids_split),sep=''))
        return(0)
    }
    
    if (list_number > num_sets) {
        #list_number <- length(chemids_split)
        stop(paste('Invalid list number. Must be less than or equal to ',num_sets,sep=''))
        return(0)
    }
    
    outloc1 <- paste(outloc, "/stage2/", sep = "")
    suppressWarnings(dir.create(outloc1))
    setwd(outloc1)
    
    # if (is.na(adduct_weights) == TRUE) {
    #     data(adduct_weights)
    #     adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
    #     adduct_weights1[1, ] <- c("M+H", 1)
    #     adduct_weights1[2, ] <- c("M-H", 1)
    #     adduct_weights <- as.data.frame(adduct_weights1)
    #     colnames(adduct_weights) <- c("Adduct", "Weight")
    # }
    
    if(!exists("adduct_weights") == T){
      #treats all adducts equally
      data("adducts_enviPat.rda")
      adduct_weights = adduct_table
    }
    
    
    #for each set of formulae (as defined by chemids_split and list_number) get the list of associated 
    #    features based on mass defect and retention time
    
    chem_score <- lapply(as.list(list_number), function(j) {
      
      #get 'Formula_X' name
      chemid <- chemids[j[1]]
      print(as.character(chemid))
      
      chemscoremat <- {}
      curmchemdata <- mchemdata[which(mchemdata$chemical_ID == chemid), ]
      curmchemdata$mz <- as.numeric(as.character(curmchemdata$mz))
      curmchemdata$time <- as.numeric(as.character(curmchemdata$time))
      curmchemdata <- as.data.frame(curmchemdata)
      curmchemdata$Module_RTclust <- gsub(curmchemdata$Module_RTclust, pattern = "_[0-9]*", replacement = "")
      isop_res_md$Module_RTclust <- gsub(isop_res_md$Module_RTclust, pattern = "_[0-9]*", replacement = "")
      isp_masses_mz_data <- {}
      
      isp_masses_mz_data <- lapply(1:length(curmchemdata$mz),function(m) {
        
        #determine which group the mz value is part of (as defined by mass defect, retention time and mz value)
        #i.e. Module_RTclust with ISgroup info also indicated
        isp_group <- as.character(curmchemdata$ISgroup[m])
        module_rt_group <- as.character(curmchemdata$Module_RTclust[m])
        module_rt_group <- gsub(module_rt_group, pattern = "_[0-9]*", replacement = "")
        
        #determine mass defect for mz value
        query_md <- curmchemdata$mz[m] - round(curmchemdata$mz[m])
        
        #isop_res_md contains information on which mz are part of which mass defect group and which Module_RTclust
        #hence, this step captures all mz with similar mass defect and retention time to those of the input mz (i.e. possible isotope peaks)
        put_isp_masses_curmz_data <- isop_res_md[which(abs((isop_res_md$MD) - (query_md)) < mass_defect_window & isop_res_md$Module_RTclust == module_rt_group), ]
        put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
        return(put_isp_masses_curmz_data)
        
      })
      
      #combine list of valid isotope candidates in to a data frame
      isp_masses_mz_data <- ldply(isp_masses_mz_data, rbind)
      
      if (is.na(mass_defect_mode) == TRUE) {
        mass_defect_mode = "pos"
      }
      
      out = get_chemscorev1.6.73_custom(chemicalid = chemid,
                                        mchemicaldata = curmchemdata,
                                        corthresh = corthresh,
                                        global_cor = global_cor,
                                        mzid = mzid,
                                        max_diff_rt = max_diff_rt,
                                        level_module_isop_annot = isp_masses_mz_data,
                                        adduct_table = adduct_table,
                                        adduct_weights = adduct_table,
                                        filter.by=c("M+H"),
                                        max_isp = max_isp,
                                        MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
                                        mass_defect_window = mass_defect_window,
                                        mass_defect_mode = mass_defect_mode,
                                        outlocorig = outloc,
                                        iso_ppm_tol = iso_ppm_tol,
                                        iso_int_tol = iso_int_tol,
                                        list_number = list_number)
      
      return(out)

    })
    
    
    
    cur_fname <- paste("chem_score", list_number, "_a.Rda", sep = "")
    chem_score2 <- chem_score[which(chem_score != "NULL")]
    curchemscoremat <- ldply(chem_score2, rbind)
    curchemscoremat = curchemscoremat[,-c(1)]
   
    cur_fname <- paste("chem_score", list_number, ".Rda", sep = "")
    save(curchemscoremat, file = cur_fname)

    Sys.sleep(1)

    #clean up cluster
    rm(chem_score)
    rm(chem_score2)
    rm( "mchemdata", "chemids", "adduct_table", 
        "global_cor", "mzid", "max_diff_rt", "isop_res_md", 
        "corthresh", "level_module_isop_annot", "chemids_split", 
        "corthresh", "max.mz.diff", "outloc", "num_sets", 
        "db_name", "num_nodes", "num_sets", "adduct_weights", 
        "filter.by")
    return(curchemscoremat)
    
    
}
