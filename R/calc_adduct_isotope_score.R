#script in which the adduct score is calculated and, if isotopes are present, is multiplied by a factor of 10


calc_adduct_isotope_score = function(k_power = 1, adduct_weights = adduct_weights, mchemicaldata = mchemicaldata, topquant_cor = NA,
                                     calc_iso_score = F, adduct_scoring_method = 'simple'){
  
  #from mchemicaldata object, extract assigned adducts
  cur_adducts_with_isotopes<-mchemicaldata$Adduct
  
  #fix formating of adduct for isotope (usually looks like M+H_[+1] and becomes M+H)
  cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-9]*\\])",replacement="")
  
  #determine how many adducts are valid according to user-specified adduct_weights information
  good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
  
  if(good_adducts_len>0){
    
    if(!is.na(topquant_cor)){
      #determine the score for this feature, based on the number of associated adducts
      chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
      conf_level = 2 #adduct info and correlation information recorded, so confidence level is higher
    }else{
      chemical_score = 1
      conf_level = 1 #only adduct info was available, therefore confidence level is just 1
    }
    
    if(adduct_scoring_method == 'simple'){
      chemical_score<-sum(chemical_score*(10^max(adduct_weights[which(adduct_weights$Adduct%in%cur_adducts_with_isotopes),]$Weight)))
    }else{
      chemical_score<-'Insert custom scoring approach here'
      #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
      #print(chemical_score)
    }
    chemical_score<-chemical_score[1]
    
    #### IMPORTANT STEP THAT FACTORS IN ISOTOPE INFORMATION, ON TOP OF THE ADDUCT-BASED SCORE
    # here, increase the chemical score if the +1 or +2 isotopes were found
    
    if(calc_iso_score==TRUE){
      check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-9]\\])")
      if(length(which(check2>0))>0){
        chemical_score<-100*chemical_score
        conf_level = 3
      }
    }
    
  }else{
    chemical_score = 0
    conf_level = 0
  }
  
  #get name of formula associated with the current chemid
  names(chemical_score)<-mchemicaldata$chemical_ID[1]
  fname<-paste(mchemicaldata$chemical_ID[1],"score.txt",sep="_")
  
  return(list("chemical_score" = chemical_score,
              "fname" = fname,
              'conf_level' = conf_level))
  
}
