get_confidence_stage2 <- function(curdata, adduct_weights = NA, max_diff_rt = 10) {
  
  curdata <- curdata[order(curdata$Adduct), ]
  cur_adducts_with_isotopes <- curdata$Adduct
    
  #a vector of adducts, including updated isotope formats
  cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
  data(adducts_enviPat)
  
  if(nrow(adduct_weights) == 0){
    adduct_weights = adduct_table
  }
  
  adduct_table <- adduct_table[order(adduct_table$Adduct), ]
  
  
  #if no valid adducts found, confidence level = None
  if (length(which(adduct_table$Adduct %in% cur_adducts)) < 1) {
    chemscoremat_conf_levels <- "None"
    score_level = 0
    return(score_level)
  }else {
    chemscoremat_conf_levels <- "High" #by default, if adducts found, score_level is started at 'High'
  }
  
  
  #merge data and adduct information, then drop column #1 (mz?)
  curdata <- cbind(curdata, cur_adducts)
  curdata <- merge(curdata, adduct_table, by.x = "cur_adducts", by.y = "Adduct")
  curdata <- curdata[, -c(1)]
  
  #get formula
  formula_vec <- curdata$Formula
  

  ########################################################################
  ####### PERFORM CHECK OF RETENTION TIME SPREAD IN RTclust MODULE #######
  ########################################################################
  
  #extract information on number of Module_RTclust exists for this input data
  table_modules <- table(curdata$Module_RTclust)
  module_names <- names(table_modules[which(table_modules > 0)])
  
  #determine the difference in largest and smallest feature retention time 
  curdata$time <- as.numeric(as.character(curdata$time))
  delta_rt <- max(curdata$time) - min(curdata$time)
  delta_rt <- round(delta_rt)
  
  #if the retention time range was larger than the max_diff_rt, but only one module found, we have medium confidence as chromatography may have been suboptimal
  # else, if more-than one module found, we have no confidence in the assignment
  # finally, if delta_rt less-than max_diff_rt and more-than module found, we have low confidence.
  if ((delta_rt > max_diff_rt) && length(module_names) == 1) {
    chemscoremat_conf_levels <- "Medium"
  } else if ((delta_rt > max_diff_rt) && length(module_names) > 1) {
    chemscoremat_conf_levels <- "None"
  } else if (length(module_names) > 1) {
    chemscoremat_conf_levels <- "Low"
  } else{
    chemscoremat_conf_levels <- "High"
  }
  
  
  ##################################################################
  ######### CHECK CHARGE STATE OF ION ##############################
  ##################################################################
  
  #check charge-state of the annotated ion forms
  min_charge <- min(curdata$charge)
  min_charge_ind <- which(curdata$charge == min_charge)
  max_int_min_charge <- max(curdata$AvgIntensity)
  
  #if the lowest charge-state is 3+, then set the confidence level to 'Medium'
  if (min_charge > 2) {
    
    #this step is needed to ensure the above check of retention times is not over-written during charge-state checking
    if(chemscoremat_conf_levels != "Low" & chemscoremat_conf_levels != "None"){
      chemscoremat_conf_levels <- "Medium"  
    }
  }

  
  ###################################################################
  ####### CHECK WHICH X-MER FORMS WERE ANNOTATED ####################
  ###################################################################
  
  #determine the lowest number of molecules among the assigned adduct forms and then determine which row corresponds to this assignment
  min_nummol <- min(curdata$num_molecules)
  min_nummol_ind <- which(curdata$num_molecules == min_nummol)
    
  # num molecules check
  if (min_nummol > 1) {
    
    # if X-mer found (from 2 to 9):
    check1 <- gregexpr(text = cur_adducts, pattern = "([2-9]+M)")
    
    #if X-mer found, then:
    if (length(check1) > 0) {
    
      #get the intensity for the adduct form with lowest-number of monomeric units
      max_int_min_mol <- max(curdata$AvgIntensity[min_nummol_ind])
      
      #generate object for recording which assignments are possibly bad annotations
      bad_ind_status <- rep(0, length(check1))
                  
      # check pattern of each adduct
      for (a1 in 1:length(check1)) {
            
        strlength <- attr(check1[[a1]], "match.length")
          
        #for di-, tri-....poly-meric ions:
        if (strlength[1] > (-1)) {
            
          #if the X-mer ion form has higher intensity than the most-intense monomer-based ion, then assume it's a bad annotation
          abund_ratio <- curdata$AvgIntensity[a1]/max_int_min_mol
            
          #if the abundance_ratio is better for the smaller X-mer ion, then leave the chemscore_conf_level at its earlier setting, else adjust to 'low'
          if (is.na(abund_ratio) == FALSE & abund_ratio > 1){
            chemscoremat_conf_levels <- "Low"
            bad_ind_status[a1] <- 1
          } 
        }
      }
    }
      
    #if in the above step there were annotations that were deemed 'bad' based on intensity check, record the bad adducts
    if (length(which(bad_ind_status == 1)) > 0) {
      bad_adducts <- cur_adducts[which(bad_ind_status == 1)]
    }
  }

  
  
  
  #######################################################################
  #### ISOTOPE CHECK - CONFIRM M+ FORM SEEN IF M+1 or M+2 ANNOTATED #####
  #######################################################################
    
  # isotope based check for +1 and +2 to make sure that the monoisotopic form is present
  
  #determine which of the adducts in curdata are isotopes
  check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
  
  #for each adduct in curdata$adduct
  for (a1 in 1:length(check2)) {
      
    #confirm which adducts are isotopes
    strlength <- attr(check2[[a1]], "match.length")
      
    #for isotopes, count the number of times the associated adduct was found (should be at least 2, one for the monoisotopic feature and one for the isotope peak)
    if (strlength[1] > (-1)) {
      
      count_abundant_form <- length(which(cur_adducts %in% cur_adducts[a1]))
      
      if (count_abundant_form < 2) {
        chemscoremat_conf_levels <- "None"
      }
      
    }
  }

  
  ################################################################
  ##### CHECK OF FORMULA ELEMENT COUNTS ##########################
  ################################################################
  
  curformula <- as.character(formula_vec[1])
  curformula <- gsub(curformula, pattern = "Ca", replacement = "")
  curformula <- gsub(curformula, pattern = "Cl", replacement = "")
  curformula <- gsub(curformula, pattern = "Cd", replacement = "")
  curformula <- gsub(curformula, pattern = "Cr", replacement = "")
  curformula <- gsub(curformula, pattern = "Se", replacement = "")
  
  numcarbons <- check_element(curformula, "C")
  
  if (numcarbons < 1) {
    chemscoremat_conf_levels <- "None"
  }
    
  
  ################################################################
  ######### SUMMARISE SCORE LEVEL CHECK ##########################
  ################################################################
    
  if (nrow(curdata) < 1) {
    score_level <- 0
  } else {
    # 3: High 2: medium 1: Low 0: None
    score_level <- as.character(chemscoremat_conf_levels)
    score_level <- gsub(score_level, pattern = "High", replacement = "3")
    score_level <- gsub(score_level, pattern = "Medium", replacement = "2")
    score_level <- gsub(score_level, pattern = "Low", replacement = "1")
    score_level <- gsub(score_level, pattern = "None", replacement = "0")
  }
    
  # print('Score is') print(score_level)
  return(score_level)
    
}
