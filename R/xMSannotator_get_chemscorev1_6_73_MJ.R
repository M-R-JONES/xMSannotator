get_chemscorev1.6.73_custom = function(chemicalid = chemid,
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
                                       iso_ppm_tol = 5,
                                       iso_int_tol = 0.15){
  
  #mchemicaldata = all data from the original dataA dataframe that is associated with a specific chemical formula from the DB
  #level_module_isop_annot = a dataframe providing details of the mz, retention time,ISgroup, RT clust module, intensity and mass defect of ALL features
  #   the ISgroup is essentially a combination of RTclust (clustering by reten. time and mz) and grouping by mass defect magnitude
  #   that is, ISgroup groups features that are similar in retention time and mass defect.
   
  #THIS FUNCTION NOW CALCULATES ISOTOPES CORRECTLY, INCLUDING FOR X-MER AND MULTIPLY CHARGED FORMS
  
  #set the working directory to the output directory defined by the user
  setwd(outlocorig)
  #set up a subdirectory for storing isotope annotation information
  outloc1<-paste(outlocorig,"/stage2/",sep="")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1) 

  #set default value for chemical scoring
  chemical_score<-(-100)
  fname<-paste(chemicalid,"data.txt",sep="_")

    
  #if(length(which(duplicated(mchemicaldata$mz)==TRUE))>0){
    #  mchemicaldata<-mchemicaldata[-which(duplicated(mchemicaldata$mz)==TRUE),]
  #}
  
  #simplify the Module_RTclust descriptor (from initial clustering procedure)
  mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
  mchemicaldata_orig<-mchemicaldata
    
  #if the adduct form associated with this chemID is valid, according to the user-defined adduct_weights obj,
  #  then retain all features associated with this adduct for further evaluation.
  if(length(which(mchemicaldata_orig$Adduct%in%as.character(adduct_weights$Adduct)))>0){
    
    #for formula currently considered, collect all valid adduct forms
    rel_adduct_data<-mchemicaldata[which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct)),]
    #for each adduct form, refine the Module_RTclust label (as defined during clustering in stage 1)
    rel_adduct_module<-gsub(rel_adduct_data$Module_RTclust,pattern="_[0-9]*",replacement="")
    #from mchemicaldata, extract data corresponding to each adduct form
    module_rt_group<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
    mchemicaldata<-mchemicaldata[which(module_rt_group%in%rel_adduct_module),]
    
  }
  
  #get all unique features (i.e. remove duplicated rows)
  mchemicaldata<-unique(mchemicaldata)
  
  #get the formula
  curformula<-as.character(mchemicaldata$Formula)
  
  #historical - used only for crude check of elements validity 
  formula_check<-getMolecule(as.character(curformula[1]))
  exp_isp<-which(formula_check$isotopes[[1]][2,]>=0.001)
  abund_ratio_vec<-formula_check$isotopes[[1]][2,exp_isp]
  
  #check appropriate number of oxygens in formula
  numoxygen<-check_element(curformula,"O")
  
  
  water_adducts<-c("M+H-H2O","M+H-2H2O","M-H2O-H")
  water_adduct_ind<-which(mchemicaldata$Adduct%in%water_adducts)
  
  #historical - I think the calculation is wrong as the elements included in the formula do not include the adduct
  #if there are no oxygens in the formula, remove associated water adduct assignments
  if(numoxygen<1){
    if(length(water_adduct_ind)>0){
      mchemicaldata<-mchemicaldata[-water_adduct_ind,]
    }
  }
  
  #for features that fail the oxygen and water adduct evaluations, return terrible score value
  if(length(mchemicaldata$mz)<1){
    chemical_score<-(-100)
    return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
  }
  
  #get set of unique adduct types associated with chemid and determine which are valid by comparing to adduct_weights table
  mchemicaldata$Adduct<-as.character(mchemicaldata$Adduct)
  mchemicaldata_goodadducts_index<-which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct))
  #legacy version:
  #uniq_adducts<-unique(mchemicaldata$Adduct)
  #good_adducts<-uniq_adducts[which(uniq_adducts%in%as.character(adduct_weights$Adduct))]
  
  #determine frequency of occurrence of each Module_RTclust in mchemicaldata object for this chemid
  #   keep those with frequency >0 and then reverse order, from highest to lowest frequency
  table_mod<-table(mchemicaldata[mchemicaldata_goodadducts_index,]$Module_RTclust)
  table_mod<-table_mod[table_mod>0]
  table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
  top_mod<-names(table_mod)
  
  #[NOT USED BELOW] determine frequency of occurrence of each Module_RTclust in the mchemicaldata object for this chemid
  #table_iso<-table(mchemicaldata$ISgroup)
  
  #store the original chemical data info for this mass-defect and tr set in mchemicaldata_origA
  mchemicaldata_origA<-mchemicaldata
  
  bool_check<-0
  
  #tidy-up the cluster module identifier in the level_module_isop_annot obj.
  level_module_isop_annot$Module_RTclust<-gsub(level_module_isop_annot$Module_RTclust,pattern="_[0-9]*",replacement="")
  
  #add new column to mchemicaldata object for capturing isotope element information
  mchemicaldata$isotope_elements = NA
  
  #object in to which all information relating to valid features will later be stored, including isotope info.
  final_isp_annot_res_all<-mchemicaldata
  
  
  ################################################################################
  ######## FIND AND ASSIGN ISOTOPE INFORMATION WITHIN THIS IF STATEMENT ##########
  ################################################################################
  
  if(length(mchemicaldata_goodadducts_index)>0){
    
    final_isp_annot_res_isp<-lapply(1:length(mchemicaldata_goodadducts_index),assign_isotopes,
                                    mchemicaldata = mchemicaldata,
                                    adduct_weights = adduct_weights, 
                                    level_module_isop_annot = level_module_isop_annot,
                                    max_diff_rt = max_diff_rt, 
                                    mass_defect_window = mass_defect_window, 
                                    mass_defect_mode = mass_defect_mode,
                                    max_isp = max_isp, 
                                    iso_int_tol = iso_int_tol, 
                                    iso_ppm_tol = iso_ppm_tol)
    
    rm(level_module_isop_annot) 
    final_isp_annot_res2<-ldply(final_isp_annot_res_isp,rbind)
    
    #NOT USED BELOW - getting total number of features in each isotope group
    #isp_group_check<-table(final_isp_annot_res2$mz)
    #good_groups<-which(isp_group_check>=max(isp_group_check)) #change to > 1 makes sense (i.e. an adduct + 1 isotope, minimally)
    #group_name<-names(isp_group_check)[good_groups]
    
    #remove the isotope grouping information
    final_isp_annot_res2<-final_isp_annot_res2[,-c(1)]
    
    rm(final_isp_annot_res_isp)
    
    #finally, add annotated isotope information to the mchemicaldata matrix
    mchemicaldata<-rbind(final_isp_annot_res_all,final_isp_annot_res2)  #[,-c(12)]
    
    ##################################################################################################
    ############################# FINISHED ASSIGNING ISOTOPE FEATURES ################################
    ##################################################################################################

  }
  
  #removed rows (corresponding to monisotopic peaks) that were inserted during addition of isotope information
  mchemicaldata<-unique(mchemicaldata)

  bad_rows<-which(is.na(mchemicaldata$mz)==TRUE)
  if(length(bad_rows)>0){
    mchemicaldata<-mchemicaldata[-c(bad_rows),]
  }
  
  mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
  
  #check which clustering modules are in the dataframe
  mod_names<-unique(mchemicaldata$Module_RTclust)
  
  #generate unique identifier for each feature in mchemicaldata
  #mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")
  
  #not used - this provides a route to extracting the correlation information for monoisotopic and isotopic peaks
  #temp_global_cor<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
  
  #generates plots showing where, chromatographically, each formula was assigned and with what frequency
  #higher frequency can mean either a single formula was assigned to multiple mz in that region, or that multiple isotopes were assigned to a single formula
  diffmatB<-lapply(1:length(mod_names),function(i){
    
    groupA_num<-mod_names[i]
    subdata<-mchemicaldata[which(mchemicaldata$Module_RTclust==groupA_num),]
    subdata<-subdata[order(subdata$time),]
    
    if(nrow(subdata)>0){
      groupB<-group_by_rt_histv2(subdata,time_step=1,max_diff_rt=10,groupnum=groupA_num)
    }else{
      groupB<-subdata
    }
    rownames(groupB)<-NULL
    return(groupB)
    
  })

  #a re-ordered form of the mchemicaldata, based on RT clust module and then by increasing retention times
  mchemicaldata<-ldply(diffmatB,rbind)

  #ensure only unique entries are stored in mchemicaldata
  mchemicaldata<-unique(mchemicaldata)
  
  #clean up enviroment
  rm(diffmatB)
  rm(final_isp_annot_res)
  
  #STORE INFORMATION RELATING TO ISOTOPE ASSIGNMENTS!
  write.table(mchemicaldata,file="../Stage2_withisotopes.txt",append=TRUE,sep="\t",col.names=FALSE)

  
  
  
  #################################################################################################
  ##  CALCULATE SCORE BASED ON ADDUCT AND ISOTOPES
  #################################################################################################
  
  #summarise which modules (from stage1A clustering procedure) remain in the mchemicaldata object
  #mchemicaldata object: contains info. relating to mz values associated with module
  # in this second stage of checking modules, the isotopes assigned above will contribute to the frequencies
  table_mod<-table(mchemicaldata$Module_RTclust)
  table_mod<-table_mod[table_mod>0]
  table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
  top_mod<-names(table_mod) #modules associated with the greatest number of compound / isotope annotations

  #NOT USED AFTER THIS POINT:
  #which IS groups (which combines the RTclust modules and the mass defect clustering) exist
  #table_iso<-table(mchemicaldata$ISgroup)

  #create backup of mchemicaldata object
  mchemicaldata_orig<-mchemicaldata
  
  #prepare for scoring
  bool_check<-0
  topquant_cor<-0
  best_conf_level<-(-100)
  k_power<-1
  
  if(length(which(table_mod>=1))>0){
    
    #by default, chemical score is negative
    best_chemical_score<-(-100) #benchmark against which chemical scores are compared
    
    #for each module associated with one-or-more annotations:
    for(i in 1:length(which(table_mod>=1))){
      
      dup_add<-{}
      chemical_score<-(-99999) #base score, against which chemical scores are compared and updated
      conf_level<-0
      
      #extract information corresponding to the module & re-order by mz value
      mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
      mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]

      
      ###### FINE TUNE LATER FOR REMOVAL OF X-MER + X-charge adducts
      
      mzid_cur = paste(mchemicaldata$mz, mchemicaldata$time, sep = '_')
      
      #get associated adducts
      cur_adducts_with_isotopes<-mchemicaldata$Adduct
      cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
      
      #remove duplicated features from the mchemicaldata dataframe (i.e. those associated with M+H, 2M+2H and 3M+3H)
      dup_features<-which(duplicated(mzid_cur)==TRUE)
      if(length(dup_features)>0){
        mchemicaldata<-mchemicaldata[-c(dup_features),]
        cur_adducts_with_isotopes = cur_adducts_with_isotopes[-dup_features]
        mzid_cur = mzid_cur[-dup_features]
      }
      
      
      
      ##############################################################################################
      ######### ASSIGN SCORE TO ADDUCT WITH THE HEIGHEST WEIGHT IN ADDUCT_WEIGHTS OBJ. #############
      ##############################################################################################
      # the adduct with the heighest weight will be used (as defined in adduct_weights$Weight)
      # contribution equals 10^[adduct weight]
      
      if(length(mchemicaldata$mz)>=1){ 
        
        #assign chemical score based on the best-weighted adduct form (here, ignores isotope peaks and peak correlations)
        assigning.score = calc_adduct_isotope_score(k_power = 1, adduct_weights = adduct_weights,mchemicaldata = mchemdata, 
                                  topquant_cor = NA, calc_iso_score = F, adduct_scoring_method = 'simple')
        
        chemical_score = as.numeric(assigning.score$chemical_score)
        fname = as.character(assigning.score$fname)
        temp_best_conf_level = as.numeric(assigning.score$conf_level)
        
        if(chemical_score>best_chemical_score){
          
          best_chemical_score<- chemical_score
          best_conf_level<-temp_best_conf_level
          best_mod_ind<-i
          best_data<-mchemicaldata
          
        }else if (chemical_score==best_chemical_score){
          
          best_chemical_score<-chemical_score
          best_conf_level<-temp_best_conf_level
          best_mod_ind<-c(i,best_mod_ind)
          best_data<-rbind(best_data,mchemicaldata)
          
        }
      }

      ############### END OF ADDUCT SCORING REGION ###################
      
      ############### BEGIN: EXTRACTION OF FEATURES CORRELATED WITH CURRENTLY-CONSIDERED FEATURE ############
      
      #check whether row in mchemicaldata corresponds to an assigned isotope (values below 0 equals no match)
      check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])")
      
      #generate identifiers for the features present in the currently-considered Module_RTclust group
      mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")

      ## extract the correlation coefficients between peaks in the module ####
      # following lines are original approach - should work fine if retention times are not rounded in earlier steps
      cor_mz<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
      cor_mz<-round(cor_mz,1)
      
      #if more than one feautre in currently-considered Module_RTclust, then matrix > 1 row.
      
      if(length(cor_mz)>1){
        
        #get the mz-time identifier for each feature extracted from global_cor
        corrownamesA<-rownames(cor_mz)
     
        #create a list with each item containing [1] mz value, [2] retention time
        mat_rownames<-strsplit(as.character(corrownamesA),split="_")

        #generate a data frame detailing the mz, retention time and correlations between features
        m1 = ldply(mat_rownames, function(x) { data.frame('mz' = x[1], 'time' = x[2] )})
        cor_mz2 = cbind(m1, cor_mz)
        cor_mz2<-cor_mz2[order(cor_mz2$mz),]
        cor_mz<-cor_mz2[,-c(1:2)]
        
      }
      
      ###############################
      # 
      # 1) for each peak, determine how many others have a correlation co-efficient with it that > corthresh (check_cor)
      # 2) determine the highest correlation value (topquent_cor)
      # 3) determine the proportion of features the feature is associated with (# of correlated features / total # peaks)
      # 4) none-relevant adduct forms are set to have 0 correlated peaks, even if correlation existed (we ignore them)
      # 5) for peaks with correlated peaks, determine which had the highest intensity (hub_mz_int)
      # 6) layer_one_associated obj. = determines which features are linked by intensity correlation and retention times
      # 7) second_level_associations obj. = features that are correlated only based on intensity (retention times are too far apart)
      
      if(length(cor_mz)>1){
        
        #determine how many mz are correlated with each peak in the correlation matrix
        check_cor<-sapply(1:dim(cor_mz)[1],function(k){
          count_mz<-length(which(cor_mz[k,]>=corthresh))-1 #-1 here is to remove the count for the feature itself
          return(count_mz)
        })
        
        #determine the maximum correlation coefficient from those above the correlation threshold
        #only considers positive correlations - this is right, we only want features that are +ve correlated with one another
        topquant_cor<-max(cor_mz[upper.tri(cor_mz)]) 
        
        if(is.na(topquant_cor)==TRUE){
          topquant_cor=0
        }
        
        #
        check_cor2<-check_cor/length(check_cor)
        
        #original form appears to set correlation of isotopes to 0
        #if(length((cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))>0){
          #check_cor[check2>0 | (cur_adducts_with_isotopes%in%adduct_weights$Adduct==FALSE)]<-0
        #}
        
        #new version - only assigns non-relevant adducts to have 0 correlation coefficients (ignores isotopes)
        if(length((cur_adducts_with_isotopes%in%adduct_weights$Adduct==FALSE))>0){
          check_cor[(cur_adducts_with_isotopes%in%adduct_weights$Adduct==FALSE)]<-0
        }
        
        #for those adduct forms with correlation coefficients greater than the cutoff:
        if(length(which(check_cor>0)==TRUE)>0){
          
          #check whether one of the assigned adducts is in the filter.by list
          hub_mz_list<-which(check_cor>0  & (cur_adducts_with_isotopes%in%filter.by==TRUE))
          
          #if no adducts match the filter.by option, check if an adduct form matches to the valid adducts defined by user
          if(length(hub_mz_list)<1){
            #the original version, below, seemed to ignore isotopes when checking which peaks were correlated
            #hub_mz_list<-which(check_cor>0 & check2<0 & (cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))
            hub_mz_list<-which(check_cor>0 & (cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))
          } 
          
          #if adduct is not in the adduct_weights list, just take the mz
          if(length(hub_mz_list)<1){
            #original form, below, seemed to ignore isotope peaks when checking number of correlated features
            #hub_mz_list<-which(check_cor>0 & check2<0)
            hub_mz_list<-which(check_cor>0)
          }
          
          #which of the correlated features has the highest intensity (if a large feature, could be the M+1 peak)
          hub_mz_int<-hub_mz_list[which(mchemicaldata$AvgIntensity[hub_mz_list]==max(mchemicaldata$AvgIntensity[hub_mz_list]))[1]]
          
          #of the annotateed features, which has the greatest number of neighbours (based on retention time)
          max_time_neighbors<-0
          best_hub_time_mz<-hub_mz_int #by default, peak with highest intensity is the hub_mz, unless something is has more neighbours!

          #determine which of the correlated features has the greatest number of correlated features
          for(h1 in hub_mz_list){

            mz_name_hub<-paste(mchemicaldata$mz[h1],mchemicaldata$time[h1],sep="_")
            
            hub_rt<-mchemicaldata$time[h1]
            
            #calculate the difference between the input feature (h1) and the other peaks
            diff_rt_hubmz<-apply(mchemicaldata,1,function(k){
              curtime = as.numeric(k[which(names(k) == 'time')])
              #curtime<-as.numeric(as.character(k[3])) #this was taking the wrong value for the retention time!
              return(abs(hub_rt-curtime))
            })
            
            #adjusted - need to remove the feature itself from the calculation of the number of time neighbours
            num_time_neighbors<-length(which(diff_rt_hubmz<=max_diff_rt)) - 1
            
            if(num_time_neighbors>max_time_neighbors){
              best_hub_time_mz<-h1
              max_time_neighbors<-num_time_neighbors
            }
          }
          
          #get mz and retention time for the hub feature
          hub_mz<-best_hub_time_mz
          hub_rt<-mchemicaldata$time[hub_mz]
          
          mz_name_hub<-paste(mchemicaldata$mz[hub_mz],hub_rt,sep="_")
          
          diff_rt_hubmz<-apply(mchemicaldata,1,function(k){
            curtime = as.numeric(k[which(names(k) == 'time')])
            #curtime<-as.numeric(as.character(k[3])) #this was taking the wrong value for the retention time!
            return(abs(hub_rt-curtime))
          })

          #layer_one_associations = features whose intensities are correlated and retention times are sufficiently close
          #stop("sf")
          if(MplusH.abundance.ratio.check==TRUE){
            #might need to be <= for mchemicaldata$AvgIntensity < mchemicaldata$AvgIntensity[hub_mz]
            layer_one_associations<-which(cor_mz[hub_mz,]>=corthresh & 
                                            mchemicaldata$AvgIntensity<mchemicaldata$AvgIntensity[hub_mz] & 
                                            diff_rt_hubmz<=max_diff_rt)
          }else{
            layer_one_associations<-which(cor_mz[hub_mz,]>=corthresh & diff_rt_hubmz<=max_diff_rt)
          }
          
          
          #NOT USED LATER IN SCRIPT
          #object for storing information on which features are associated ONLY BY CORRELATIONS (i.e. not by rt and correlation of intensities
          #original code (below) only subsetted those features in layer_one_associations, all of which were > corthresh anyway
          #I believe the correct subsetting here should be to extract features from cor_mz (not layer_one_associations) > 0.7
          #for(l in layer_one_associations){
          #  second_level_associations<-c(which(cor_mz[l,]>=0.7))
          #}
          
          #NOT USED LATER IN SCRIPT:
          #matrix showing features correlated with the hub_mz feature (i.e. the one with highest intensity and with greater number of tr neighbours)
          #sub_cor_mz<-cor_mz[hub_mz,which(cor_mz[hub_mz,]>=corthresh)]

          
          #subset the mchemicaldata object to keep only those features that are associated with one another!
          selected_mz<-c(hub_mz,layer_one_associations) #intersect(layer_one_associations,second_level_associations))
          selected_mz<-unique(selected_mz)
          
          #if more-than one feature is correlated with at least one other feature
          if(length(which(check_cor2>= 1/length(cor_mz)))>0){
            mchemicaldata<-mchemicaldata[selected_mz,]
          }else{
            mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh),]
          }
          
          #remove empty rows (ignore isotope_elements column during checks)
          keep.indx = which(complete.cases(mchemicaldata[,c(1:(ncol(mchemicaldata)-1))])) #check which rows are not empty
          mchemicaldata = mchemicaldata[keep.indx,]
          rm(keep.indx)
          
          #remove rows from mchemicaldata if no time data is recorded (surely this should have happened at start of processing?)
          if(length(which(is.na(mchemicaldata$time))==TRUE)>0){
            mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$time)==TRUE),]
          }
          
          #if fewer-than 2 rows in mchemicaldata, no need to consider feature further as nothing is associated with it
          #  as such, no score can be calculated
          if(length(mchemicaldata$mz)<2){
            next #this escapes the top_mod for-loop, with the best_chemical_score and top_mod_index having already been recorded above (i.e. just adduct-based scoring)
          }
          
          #calculate difference between min. and max. retention times in mchemicaldata features
          diff_rt<-round(abs(min(as.numeric(mchemicaldata$time))-max(as.numeric(mchemicaldata$time))))

          
          #store information on duplicate rows
          dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
          if(length(dup_add)>0){
            dup_data<-mchemicaldata[c(dup_add),]
            #mchemicaldata<-mchemicaldata[-c(dup_add),]
          }
          
          dup_mz_check<-dup_add
          
          #check if diff_rt less-than user-defined max_diff_rt value
          if(diff_rt<=max_diff_rt){
            
            
            #as there are > 1 row in the mchemicaldata dataframe, check to see if isotope was assigned
            calc_iso_score = T
            
            assigning.score = calc_adduct_isotope_score(k_power = 1, adduct_weights=adduct_weights, 
                                                        mchemicaldata = mchemicaldata, topquant_cor = topquant_cor, 
                                                        calc_iso_score = calc_iso_score, adduct_scoring_method = 'simple')
              
            chemical_score = as.numeric(assigning.score$chemical_score)
            fname = as.character(assigning.score$fname)
            temp_best_conf_level = as.numeric(assigning.score$conf_level)
            
            if(chemical_score>best_chemical_score){
              
              best_chemical_score<- chemical_score
              best_conf_level<-temp_best_conf_level
              best_mod_ind<-i
              best_data<-mchemicaldata
              
            }else if (chemical_score==best_chemical_score){
              
              best_chemical_score<-chemical_score
              best_conf_level<-temp_best_conf_level
              best_mod_ind<-c(i,best_mod_ind)
              best_data<-rbind(best_data,mchemicaldata)
              
            }
     
            #return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
          } else{

            #print("fixing time")
          
            mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
            
            #updated approach, which includes a column for storing the isotope_elements information
            mchemicaldata<-cbind(mchemicaldata[,c(2:11)],mchemicaldata[,1],mchemicaldata[,c(12:15)])
            colnames(mchemicaldata)<-c("mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","AvgIntensity", "MD", "isotope_elements")
            
            #get the identifier of the RTclust module (top_mod[i] under consideration)
            groupnumA<-unique(mchemicaldata$Module_RTclust)
            
            #determine the density of a given compound/formula annotation along chromatographic axis (and plot)
            mchemicaldata<-group_by_rt_histv2(mchemicaldata,time_step=1,max_diff_rt=max_diff_rt,groupnum=groupnumA)
            
            #determine frequency of occurrence of each RTclust module in mchemicaldata obj. and extract their names
            top_mod_sub<-table(mchemicaldata$Module_RTclust)
            top_mod_sub_names<-names(top_mod_sub)
            
            #determine which RTclust module occurred most frequently (most number of associated features)
            #  and subset mchemicaldata to keep only this module. Finally, order by mz
            max_top_mod<-which(top_mod_sub==max(top_mod_sub))[1]
            mchemicaldata<-mchemicaldata[which(mchemicaldata$Module_RTclust==top_mod_sub_names[max_top_mod]),]
            mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]

            ######### HERE: adjust retention time tolerance if features within RTclust module are separated by more-than max_diff_rt variable
            
            #preferred approach would be to hierarchically cluster the data and split according to max_diff_rt (more than one group generated if distance is more than max_diff_rt)
            #cuts = cutree(hclust(dist(mchemicaldata$time)), h = max_diff_rt)
            #time_cor_groups = as.list(split(mchemicaldata, cuts))
            
            
            ### effectively performing outlier check here, with features above or below IQR*1.5 being considered outliers
            ##   caution is applied, as if IQR is too large (due to skewed rt values &/or outliers) the iqr is refined
            
            #summarise distribution of data and extract the interquartile range (iqr1)
            s1<-summary(mchemicaldata$time)
            iqr1<-s1[5]-s1[2]
            
            # adjust retention time of features with lowest and highest tr in group (expand using 1.5* IQR as expansion factor)
            min_val<-s1[2]-(1.5*iqr1)
            max_val<-s1[5]+(1.5*iqr1)
            
            #check whether BOTH min_val and max_val have been refined further than their original values 
            # if so, it's likely that the data is skewed (should not happen for normally-distributed values)
            # When this happens the min_ and max_ vals are refined once more, this time using a potentially smaller iqr1 value (accounts for skewed data & outliers)
            if(min_val<s1[1] && max_val>s1[6]){
              iqr1<-min(abs(s1[3]-s1[2]),abs(s1[3]-s1[5]))
              min_val<-s1[2]-(1.5*iqr1)
              max_val<-s1[5]+(1.5*iqr1)
            }
            
            # if the min_val was less-than the original min_val in the feature set, then reset to original min_val
            if(min_val<s1[1]){
              min_val<-s1[1]
            }
          
            # if the max_val was greater-than the original max-val in the feature set, then reset to original max_val
            if(max_val>s1[6]){
              max_val<-s1[6]
            }
            
            #determine range of times for features in RTclust
            diff_rt<-abs(max_val - min_val) #used to be: abs(max(mchemicaldata$time) - min(mchemicaldata$time))
            
            #re-define the iqr1 as the smaller of either median - 1st quartile, or 3rd quartile - median
            #   presumably this is done to account for skewed data, i.e. taking the least spread side of the distribution as 0.5* iqr
            iqr1<-min(abs(s1[3]-s1[2]),abs(s1[3]-s1[5]))
            
            #if the newly-established iqr1 is larger than the original max_diff_rt, update, otherwise keep max_diff_rt
            iqr1<-max(2*iqr1,max_diff_rt)
            
            print(paste('New max_diff_rt value set at: ', iqr1, sep = ''))
            
            if(nrow(mchemicaldata)<1){
              next
            }
            
            #if the range is more-than 2*newly-established iqr1, then split data in to chunks
            #  splitting done from min_val - max_diff_rt, to max_val + max_diff_rt, splitting in to chunks of size iqr1
            #  likely leads to sub-optimal splitting of features
            if(diff_rt>=iqr1){

              #earlier version: if(diff_rt>2*max_diff_rt){
              
              #it is not clear why iqr1 is used for binning features here.
              time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                split(x,cut(mchemicaldata$time,breaks=seq(min_val-iqr1,max_val+iqr1,iqr1)))
                #split(x,cut(mchemicaldata$time,breaks=seq(min_val-max_diff_rt,max_val+max_diff_rt,iqr1))) #earlier version
              })
              
            }else{
              
              #for features that elute closely, this part applies
              max_val<-max(mchemicaldata$time)
              min_val<-min(mchemicaldata$time)
              
              #given that diff_rt is less-than max_diff_rt, all peaks are grouped together in a single set
              time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                split(x,cut(mchemicaldata$time,breaks=seq(min_val-iqr1,max_val+iqr1,iqr1)))
              })
            }

            #determine how many features are in each retention time bin/group
            group_sizes<-sapply(time_cor_groups,function(x){dim(as.data.frame(x))[1]})

            #an empty container that is used for storing the index of peak groups that have at least one valid adduct
            good_temp<-c()
            group_ind_size = 0

            #for each subset of features (in a group) determine whether valid adducts were found (as defined by adduct_weights obj.)
            for(g1 in 1:length(group_sizes)){
              #if the number of adducts is greater than 0, append to good_temp object
              tempdata<-time_cor_groups[[g1]]
              check_reladd<-which(tempdata$Adduct%in%as.character(adduct_weights$Adduct))
              if(length(check_reladd)>0 & nrow(tempdata)>0){
                good_temp<-c(good_temp,g1)
                if(nrow(tempdata)>group_ind_size){
                  #if the number of features in the group exceeds 1, add g1 identifier to the group_ind_val
                  group_ind_val<-g1
                  group_ind_size<-nrow(tempdata)
                }
              }
            }
            
            # if there are two features in the group
            if(length(good_temp) > 0){
              
              temp_best_score<-(-100)
              temp_best_data<-{}
              
              #for each group
              for(g2 in good_temp){
                
                mchemicaldata<-time_cor_groups[[g2]]
                
                diff_rt<-round(abs(max(as.numeric(mchemicaldata$time)) - min(as.numeric(mchemicaldata$time))))

                #if at least two rows in the group of features
                if(dim(mchemicaldata)[1]>1){
                  calc_iso_score = T
                }else{
                  calc_iso_score = F
                }
                  
                assigning.score = calc_adduct_isotope_score(k_power = 1, adduct_weights = adduct_weights,
                                                            mchemicaldata = mchemicaldata, topquant_cor = topquant_cor, 
                                                            calc_iso_score = calc_iso_score, adduct_scoring_method = 'simple')
                
                chemical_score = as.numeric(assigning.score$chemical_score)
                fname = as.character(assigning.score$fname)
                temp_best_conf_level = as.numeric(assigning.score$conf_level)
                
                if(chemical_score>best_chemical_score){
                  
                  best_chemical_score<- chemical_score
                  best_conf_level<-temp_best_conf_level
                  best_mod_ind<-i
                  best_data<-mchemicaldata
                  
                }else if (chemical_score==best_chemical_score){
                  
                  best_chemical_score<-chemical_score
                  best_conf_level<-temp_best_conf_level
                  best_mod_ind<-c(i,best_mod_ind)
                  best_data<-rbind(best_data,mchemicaldata)
                  
                }
              }
            }else{
              next
            }
          }
        }
      }
    }
  }
  
  ################################################################
  ######## FINAL CLEANUP BEFORE RETURNING CHEMICAL SCORE #########
  ################################################################
  
  # if(chemical_score>min_chemical_score){
  #   chemical_score<-chemical_score*(conf_level^conf_level)
  # }else{
  #   chemical_score<-0
  # }
  # 
  # if(length(dup_add)>0){
  #   mchemicaldata<-rbind(mchemicaldata,dup_data)
  # }
  # 
  # if(is.na(conf_level)==TRUE){
  #   conf_level<-0
  # }
  # 
  # if(is.na(chemical_score)==TRUE){
  #   chemical_score<-0
  # }
  
  mchemicaldata = best_data
  score_level = get_confidence_stage2(curdata = mchemicaldata, adduct_weights=adduct_weights)
  
  mchemicaldata = cbind(mchemicaldata, 'conf_level' = score_level, 'chemical_score' = best_chemical_score)
  
  
  
  
  rm("mzid","global_cor","temp_global_cor")
  return(list("chemical_score"=best_chemical_score,"filtdata"=mchemicaldata))
}