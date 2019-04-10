get_chemscorev1.6.71_custom = function(chemicalid = chemid,
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
  
  #THIS FUNCTION NOW CALCULATES ISOTOPES CORRECTLY, INCLUDING FOR X-MER AND MULTIPLY CHARGED FORMS
  
  #set the working directory to the output directory defined by the user
  setwd(outlocorig)
  #set up a subdirectory for storing isotope annotation information
  outloc1<-paste(outlocorig,"/stage2/",sep="")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1) 
  
  #this is the object that is submitted for isotope annotation
  level_module_isop_annot$mz<-as.numeric(as.character(level_module_isop_annot$mz))
  level_module_isop_annot$time<-as.numeric(as.character(level_module_isop_annot$time))
  
  #set default value for chemical scoring
  chemical_score<-(-100)
  fname<-paste(chemicalid,"data.txt",sep="_")
  
  
  if(length(mchemicaldata$mz)>0){
    
    if(length(which(duplicated(mchemicaldata$mz)==TRUE))>0){
      #  mchemicaldata<-mchemicaldata[-which(duplicated(mchemicaldata$mz)==TRUE),]
    }
    
    #simplify the Module_RTclust descriptor (from initial clustering procedure)
    mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
    mchemicaldata_orig<-mchemicaldata
    
    
    if(length(which(mchemicaldata_orig$Adduct%in%as.character(adduct_weights$Adduct)))>0){
      
      #for formula currently considered, collect all valid adduct forms
      rel_adduct_data<-mchemicaldata[which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct)),]
      #for each adduct form, refine the Module_RTclust label (as defined during clustering in stage 1)
      rel_adduct_module<-gsub(rel_adduct_data$Module_RTclust,pattern="_[0-9]*",replacement="")
      #from mchemicaldata, extract data corresponding to each adduct form
      module_rt_group<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
      mchemicaldata<-mchemicaldata[which(module_rt_group%in%rel_adduct_module),]
      
    }
    
    #get all unique features
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
    if(nrow(mchemicaldata)<1){
      chemical_score<-(-100)
      return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
    }
    
    #get set of unique adduct types associated with chemid
    mchemicaldata$Adduct<-as.character(mchemicaldata$Adduct)
    uniq_adducts<-unique(mchemicaldata$Adduct)
    
    #keep only adducts that were defined as valid adducts in the adduct_weights object
    good_adducts<-uniq_adducts[which(uniq_adducts%in%as.character(adduct_weights$Adduct))]
    
    #from the full list of chemical annotations, extract those corresponding to the set of formulae in this loop
    table_mod<-table(mchemicaldata$Module_RTclust)
    
    #get the IS group for this group of formulae
    table_iso<-table(mchemicaldata$ISgroup)
    
    #keep only those Module_RTclusts which have value greater than 0
    table_mod<-table_mod[table_mod>0]
    table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
    
    #store the original chemical data info for this mass-defect and tr set in mchemicaldata_origA
    mchemicaldata_origA<-mchemicaldata
    
    #corthresh<-0.5
    top_mod<-names(table_mod)
    
    bool_check<-0
    
    #tidy-up the cluster module identifier
    level_module_isop_annot$Module_RTclust<-gsub(level_module_isop_annot$Module_RTclust,pattern="_[0-9]*",replacement="")
    
    #determine of the adducts in the mchemicaldata object are in the list of valid adducts
    mchemicaldata_goodadducts_index<-which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct))
    
    #add new column for capturing isotope element information
    mchemicaldata$isotope_elements = NA
    
    #object in to which isotope information will be stored
    final_isp_annot_res_all<-mchemicaldata
    
    #temporary object used to capture isotope annotation outputs
    final_isp_annot_res_isp<-{}
    
    if(length(mchemicaldata_goodadducts_index)>0){
      
      #for all good adducts, try to find the associated isotopes
      final_isp_annot_res_isp<-lapply(1:length(mchemicaldata_goodadducts_index),function(i){
        
        m<-mchemicaldata_goodadducts_index[i]
        
        #from mchemicaldata, extract the row corresponding to the input adduct index
        final_isp_annot_res<-cbind(paste("group",i,sep=""),mchemicaldata[m,])
        #get the ISgroup corresponding to the input adduct
        isp_group<-as.character(mchemicaldata$ISgroup[m])
        #get the cluster (based on mz and retention) for currently-considered mz
        module_rt_group<-as.character(mchemicaldata$Module_RTclust[m])
        module_rt_group<-gsub(module_rt_group,pattern="_[0-9]*",replacement="")
        
        ##########################################################################################
        #from the formula, calculate isotope profile
        
        ##### IMPORTANT - NEED TO FIND A WAY OF GENERATING ISOTOPIC INFORMATION FOR DIMERS / TRIMERS etc
        
        #here, set form to the formula for the assigned suspect
        form = final_isp_annot_res$Formula
        
        #get reference to row in adducts_weights table corresponding to current adduct
        indx_add_table = which(adduct_weights$Adduct == as.character(mchemicaldata[m,]$Adduct))
        print('found index')
        
        #multiply formula weight by number of monomers assigned in adduct
        nummols = as.numeric(adduct_weights[indx_add_table,]$num_molecules)
        
        #get the charge state associated with the adduct form
        z = as.numeric(adduct_weights[indx_add_table,]$charge)
        
        #use enviPat multiform function to calculate total formula of X-mer
        #form = enviPat::multiform(curformula[1], fact = nummols)
        form = enviPat::multiform(check_chemform(isotopes, curformula[1])$new_formula, fact = nummols)
        
        print(paste('Formula for compound: ', form, sep = ''))
        print(paste('adduct for compound: ', as.character(adduct_weights$Adduct[indx_add_table]), sep =''))
        print(paste('charge for compound: ', z, sep = ''))
        
        #correct formula to take in to account adduct elements
        if(!is.na(adduct_weights[indx_add_table,]$Merge_add)){
          form = enviPat::mergeform(form, adduct_weights[indx_add_table,]$Merge_add)
        }
        
        if(!is.na(adduct_weights[indx_add_table,]$Merge_sub)){
          form = enviPat::subform(form, adduct_weights[indx_add_table,]$Merge_sub)
        }
        
        print(paste('Sum formula for compound and associated adduct is: ', form, sep = ''))
        
        ##########################################################################################
        
        #from enviPat package, read in isotopes table
        data(isotopes)
        
        max_isp_count = max_isp
        
        #calculate the fine isotope structure for the formula, taking in to account charge-state and adduct elements
        mol_isos = enviPat::isopattern(isotopes, chemforms = form, charge = z, threshold = 0.1)
        
        ### As the adduct table has factored in the mass of electrons, correction of theoretical mz is not required
        
        #correct for slight shift in mass of isotope m/z values due to the charge
        #if(z>0 & adduct_weights[indx_add_table,]$Mode == 'positive'){
        #  
        #  mol_isos[[1]][,1] = mol_isos[[1]][,1] - (z * 0.00054858)
        #  
        #}else if(z>0 & adduct_weights[indx_add_table,]$Mode == 'negative'){
        #  
        #  mol_isos[[1]][,1] = mol_isos[[1]][,1] + (z * 0.00054858)
        #  
        #}
        
        #keep only isotope peaks, not the monoisotopic peak!
        top_isotopes = mol_isos[[1]][order(mol_isos[[1]][,2], decreasing = T)[1:max_isp+1],] #the +1 is to account for the monoisotopic peak
        
        top_isotopes = data.frame(top_isotopes, check.names = F)
        
        #get the string identifier for each top isotope
        top_isotopes$name = apply(top_isotopes, 1, function(x){
          
          colnm = names(x)
          stng = list()
          
          for(c in 1:length(colnm)){
            if(x[c] > 0){
              stng[[c]] = paste(colnm[c], x[c], sep = '')
            }
          }
          
          paste(unlist(stng[3:length(stng)]), collapse = '_')
          
        })
        
        isp_mat_module_rt_group<-as.character(level_module_isop_annot$Module_RTclust)
        
        query_md<-mchemicaldata$mz[m]-round(mchemicaldata$mz[m])
        query_rt<-mchemicaldata$time[m]
        
        #to account for fact that [M+1] can have a higher intensity than the monoisotopic peak,
        #first determine if [M+1] expected abundance exceeds the monoisotopic abundance and if so
        #calculate the ratio between the [M+1] and M+ peak abundances and multiply by the M+ peak intensity
        #add a further 'cushion' to ensure more-abundant features are collected during the subset procedure
        
        max_isotope_abundance = max(top_isotopes$abundance)
        if(max_isotope_abundance < 100){
          query_int<-1*(mchemicaldata$AvgIntensity[m])
        }else{
          #multiply the compound intensity by the ratio of its intensity to that of its largest isotope
          ratio = max(top_isotopes$abundance) / mchemicaldata$AvgIntensity[m]
          query_int = ratio * mchemdata$AvgIntensity[m]
          #this is a tolerance factor that adds 10% further to the intensity value to ensure largest isotope captured
          query_int = query_int + (0.1 * query_int)
        }
        
        #subset procedure
        #A VERY IMPORTANT STEP WHERE ALL POSSIBLE ISOTOPE CANDIDATES ARE EXTRACTED FROM THE level_module_isop_annot OBJECT
        put_isp_masses_curmz_data<-level_module_isop_annot[which(abs(level_module_isop_annot$time-query_rt)<max_diff_rt & 
                                                                   abs((level_module_isop_annot$MD)-(query_md))<mass_defect_window & 
                                                                   isp_mat_module_rt_group==module_rt_group & 
                                                                   level_module_isop_annot$AvgIntensity<=query_int),]
        
        put_isp_masses_curmz_data<-as.data.frame(put_isp_masses_curmz_data)
        put_isp_masses_curmz_data$mz<-as.numeric(as.character(put_isp_masses_curmz_data$mz))
        put_isp_masses_curmz_data$time<-as.numeric(as.character(put_isp_masses_curmz_data$time))
        mchemicaldata<-as.data.frame(mchemicaldata)
        
        put_isp_masses_curmz_data<-unique(put_isp_masses_curmz_data)
        put_isp_masses_curmz_data<-put_isp_masses_curmz_data[order(put_isp_masses_curmz_data$mz),]
        
        if(nrow(put_isp_masses_curmz_data)>0){
          
          #calculate the relative intensities of peaks versus the monoisotopic mass
          int_vec<-put_isp_masses_curmz_data$AvgIntensity
          int_vec<-int_vec/mchemicaldata$AvgIntensity[m] 
          
          ###no idea why these elements are excluded###
          #curformula<-gsub(curformula,pattern="Sr",replacement="")
          #curformula<-gsub(curformula,pattern="Sn",replacement="")
          #curformula<-gsub(curformula,pattern="Se",replacement="")
          #curformula<-gsub(curformula,pattern="Sc",replacement="")
          #curformula<-gsub(curformula,pattern="Sm",replacement="")
          
          #if(FALSE){
          #  numchlorine<-check_element(curformula,"Cl")
          #  numsulphur<-check_element(curformula,"S")
          #  numbromine<-check_element(curformula,"Br")
          #  max_isp_count<-max(numchlorine,numsulphur,numbromine,max_isp)
          #}
          
          #max_isp_count<-max(exp_isp)
          
          if(is.na(max_isp_count)==TRUE){
            max_isp_count=1
          }
          
          #previous version
          #ischeck<-which(int_vec<=max(abund_ratio_vec[-c(1)]+0.10))
          
          #ADDED HERE - earlier part was based on incorrect Rdisop calculations without adduct and X-mer
          abund_ratio_vec = top_isotopes$abundance / 100
          ischeck<-which(int_vec<=max(abund_ratio_vec+0.20))
          
          put_isp_masses_curmz_data$time<-as.numeric(as.character(put_isp_masses_curmz_data$time))
          
          if(length(ischeck)>0){
            
            for(rnum in 1:length(ischeck)){
              temp_var<-{}
              bool_check<-1
              
              isp_v<-ischeck[rnum]
              
              isp_v<-as.numeric(as.character(isp_v))
              #print(isp_v)
              #print(put_isp_masses_curmz_data[isp_v,])
              
              # print(mchemicaldata[m,2])
              diff_rt<-abs(put_isp_masses_curmz_data$time[isp_v]-mchemicaldata$time[m])
              
              #### MODIFICATION
              
              #if the candidate isotope is indeed an isotope, this defined the +/- value associated with t
              #e.g. if isnum is 1, this is the M+1 isotope
              isnum = round((put_isp_masses_curmz_data$mz[isp_v]-mchemicaldata$mz[m])*z)
              
              #determines whether peak is +[integer] form of isotope, e.g. [M+1]
              #isnum<-(round(put_isp_masses_curmz_data[isp_v,1])-round(mchemicaldata[m,1]))
              bool_check<-1
              #print(mass_defect_mode)
              
              if(mass_defect_mode=="neg" | mass_defect_mode=="both"){
                isnum<-abs(isnum)
              }else{
                if(isnum>0){
                  bool_check<-1
                }else{
                  bool_check<-0
                }
              }
              
              #if X in M+X is less than the maximum number of allowed isotopes, perform assignment procedure
              if(diff_rt<max_diff_rt & isnum<=max_isp){
                
                max_isp_count=max_isp
                
                if(max_isp_count>0 && isnum<=max_isp_count && bool_check>0){
                  
                  
                  ###### CHANGED HERE!!!!!!#######
                  #isnum2<-(round(put_isp_masses_curmz_data[isp_v,1])-round(mchemicaldata$MonoisotopicMass[m]))
                  #isnum2<-(round(put_isp_masses_curmz_data[isp_v,1])-round(mchemicaldata$mz[m]))
                  #isnum2<-round(isnum2)
                  
                  #ADDED MULTIPLICATION TO ACCOUNT FOR CHARGE STATE OF ADDUCT
                  #doubly-charged ions that are 1 Da apart are +2 isotope.
                  isnum2 = round((put_isp_masses_curmz_data$mz[isp_v]-mchemicaldata$mz[m])*z)
                  
                  if(isnum2<=0){
                    isp_sign<-"-"
                  }else{
                    isp_sign<-"+"
                  }
                  
                  isnum2<-abs(isnum2)
                  if(isnum2<=max_isp){
                    
                    #form_name<-as.character(paste(mchemicaldata[m,7],"_[+",(isnum),"]",sep=""))
                    
                    top_isotopes
                    
                    #iso_ppm_tol = 5
                    
                    
                    #### FIRST CHECK PPM FROM POSSIBLE ISOTOPE ASSIGNMENTS
                    iso_hits = ((put_isp_masses_curmz_data$mz[isp_v] - top_isotopes$`m/z`) / top_isotopes$`m/z`)* 1000000
                    iso_hits = top_isotopes[which(abs(iso_hits) < iso_ppm_tol),]
                    
                    #IF > 1 ISOTOPE WITHIN PPM TOLERANCE, CHECK IF ANY WITHIN ISO_INT_TOL
                    #start by calculating ratio of candidate isotope peak to monoisotopic peak
                    percent = (put_isp_masses_curmz_data$AvgIntensity[isp_v] / mchemicaldata$AvgIntensity[m])*100
                    
                    ##### IMPORTANT PARAMETER ###########
                    #iso_int_tol = 0.05 #5% tolerance allowed
                    
                    #collect isotopes that occur within +/- the intensity tolerance
                    iso_hits = subset(iso_hits, 
                                      iso_hits$abundance > (percent - (percent * iso_int_tol)) &
                                        iso_hits$abundance < (percent + (percent * iso_int_tol)))
                    
                    #FORM NAME IS FOR NEUTRAL MOLECULE AS DEFINED IN THE USER-SELECTED DATABASE!
                    if(nrow(iso_hits) > 0){
                      
                      #provides details of the elements contributing to the isotope
                      isotope_elements = paste(iso_hits$name, collapse='//', sep='')
                      
                      #provides the CcHhNnOoPpSs_[+X] isotope form (neutral)
                      form_name = as.character(paste(mchemicaldata$Formula[m],"_[",isp_sign,(isnum2),"]",sep=""))
                      
                      #COMMENTED THIS LINE!!!!
                      #form_name<-as.character(paste(mchemicaldata[m,7],"_[",isp_sign,(isnum2),"]",sep=""))
                      
                      #end of modification
                      ###########################################################
                      
                      other_inf<-cbind(rep("-",7))
                      temp_var<-cbind(put_isp_masses_curmz_data[isp_v,c(1:2)],
                                      t(other_inf),
                                      put_isp_masses_curmz_data[isp_v,c(3:4)],
                                      put_isp_masses_curmz_data[isp_v,c(2,5:6)],
                                      isotope_elements) #ADDED - provides details of which elements contributed to the putative isotope annotation
                      
                      temp_var<-as.data.frame(temp_var)
                      
                      colnames(temp_var)<-colnames(mchemicaldata)
                      temp_var$Formula<-form_name
                      temp_var$Name<-as.character(mchemicaldata$Name[m])
                      temp_var$chemical_ID<-as.character(mchemicaldata$chemical_ID[m])
                      
                      #temp_var$Adduct<-paste(mchemicaldata[m,9],"_[+",isnum,"]",sep="")
                      
                      #get the adduct associated with the assigned isotopes
                      adductname=mchemicaldata$Adduct[m] #queryadductlist[adnum]
                      #adductmass=adduct_table[as.character(adductname),4]
                      
                      #temp_var$Adduct<-paste("M","_[",isp_sign,(abs(isnum2)),"]",sep="")
                      
                      ###### MODIFIED THIS - incorrectly reference 'isnum' versus 'isnum2'
                      temp_var$Adduct<-paste(mchemicaldata$Adduct[m],"_[",isp_sign,(abs(isnum2)),"]",sep="")
                      
                      temp_var<-as.data.frame(temp_var)
                      temp_var<-cbind(paste("group",i,sep=""),temp_var)
                      final_isp_annot_res<-as.data.frame(final_isp_annot_res)
                      
                      if(nrow(temp_var)>0){
                        
                        check_mz<-which(temp_var$mz%in%final_isp_annot_res$mz)
                        
                        if(length(check_mz)>0){
                          temp_var<-temp_var[-c(check_mz),]
                        }
                        
                        if(nrow(temp_var)>0){
                          
                          final_isp_annot_res<-rbind(final_isp_annot_res,temp_var)
                        }
                      
                        return(final_isp_annot_res)
                      
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        #return(final_isp_annot_res)
        
      })
      
      
      rm(level_module_isop_annot) 
      final_isp_annot_res2<-ldply(final_isp_annot_res_isp,rbind)
      
      #get total number of features in each isotope group
      isp_group_check<-table(final_isp_annot_res2[,1])
      
      #earlier form
      #good_groups<-which(isp_group_check==max(isp_group_check))
      good_groups<-which(isp_group_check>=max(isp_group_check))
      
      group_name<-names(isp_group_check)[good_groups]
      final_isp_annot_res2<-as.data.frame(final_isp_annot_res2)
      
      #remove the isotope grouping information
      final_isp_annot_res2<-final_isp_annot_res2[,-c(1)]
      final_isp_annot_res2<-as.data.frame(final_isp_annot_res2)
      
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
    mod_names<-mchemicaldata$Module_RTclust
    mod_names<-unique(mod_names)

    #mchemicaldata$mz = round(mchemicaldata$mz, 5)
    #mchemicaldata$time = round(mchemicaldata$time, 1)
    
    mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")
    
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

    #a re-ordered form of the mchemicaldata, based on retention times
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
    table_mod<-table(mchemicaldata$Module_RTclust)
    table_mod<-table_mod[table_mod>0]
    table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
    top_mod<-names(table_mod) #modules associated with the greatest number of compound / isotope annotations

    #which IS groups (which combines the RTclust modules and the mass defect clustering) exist
    table_iso<-table(mchemicaldata$ISgroup)

    #create backup of mchemicaldata object
    mchemicaldata_orig<-mchemicaldata
    
    #################
    # restart from here if required
    #load("~/temp/outputs/env_isotope_assigned_anabaenopeptin_finishScoreSection.rdata")
    #load("~/temp/outputs/global_cor.rda")
    #load("~/temp/outputs/xMSannotator_levelA_modules.rda")
    
    #fix the mzid rounding from earlier steps - TEMPORARY
    #levelA_res$mz = round(levelA_res$mz, 5)
    #levelA_res$time = round(levelA_res$time, 1)
    #mzid = paste(levelA_res$mz, levelA_res$time, sep = '_')
    
    
    #prepare for scoring
    bool_check<-0
    topquant_cor<-0
    best_conf_level<-(-100)
    k_power<-1
    
    
    if(length(which(table_mod>=1))>0){
      
      #by default, chemical score is negative
      best_chemical_score<-(-100)
      
      #for each module associated with one-or-more annotations:
      for(i in 1:length(which(table_mod>=1))){
        
        dup_add<-{}
        chemical_score<-(-99999) #benchmark against which chemical scores are compared
        conf_level<-0
        
        top_mod_cleaned = strsplit(top_mod[i], '_')[[1]][1] #added
        
        #extract information corresponding to the module & re-order by mz value
        #mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),] #earlier approach
        mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
        mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
        
        #get associated adducts
        cur_adducts_with_isotopes<-mchemicaldata$Adduct
        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
        
        #if only 1 row associated with the annotation
        #if(nrow(mchemicaldata)<2){ # this is the earlier check - makes no sense as it only considered best scores when a single adduct assigned
        
        
        ############ SECTION IN WHICH ADDUCT CONTRIBUTION TO OVERALL CHEMICAL SCORE IS DEFINED #############
        # the adduct with the heighest weight will be used (as defined in adduct_weights$Weight)
        # contribution equals 10^[adduct weight]
        
        if(length(mchemicaldata$mz)>=1){ 
          
          #confirm adduct existed in the valid adducts defined in the adduct_weights object
          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
          
          #go through the adducts and determine which offers the best score
          if(good_adducts_len>0){
            
            #calculate the chemical score for this single annotation
            chemical_score<-(1*(10^max(adduct_weights[which(adduct_weights$Adduct%in%cur_adducts_with_isotopes),]$Weight )))
            chemical_score<-chemical_score[1]
            
            if(chemical_score>best_chemical_score){
              
              best_chemical_score<-chemical_score
              best_conf_level<-1
              best_mod_ind<-i
              best_data<-mchemicaldata
              
            }else if (chemical_score==best_chemical_score){
              
              best_chemical_score<-chemical_score
              best_conf_level<-1
              best_mod_ind<-c(i,best_mod_ind)
              best_data<-rbind(best_data,mchemicaldata)
              
            }
          }
        }

        ############### END OF ADDUCT SCORING REGION ###################
        
        ############### BEGIN: EXTRACTION OF FEATURES CORRELATED WITH CURRENTLY-CONSIDERED FEATURE ############
        
        #
        check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])")
        
        #generate identifiers for the features present in the currently-considered Module_RTclust group
        mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")
        
        #remove duplicated features from the mchemicaldata dataframe (i.e. those associated with M+H, 2M+2H and 3M+3H)
        dup_features<-which(duplicated(mzid_cur)==TRUE)
        if(length(dup_features)>0){
          mchemicaldata<-mchemicaldata[-c(dup_features),]
          mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
        }
        
        
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
          #only considers positive correlations - this is right, we only want features that are correlared with one another
          topquant_cor<-max(cor_mz[upper.tri(cor_mz)]) 
          
          #
          check_cor2<-check_cor/length(check_cor)
          
          #set none-relevant adduct forms to have a 0 correlated features
          if(length((cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))>0){
            check_cor[check2>0 | (cur_adducts_with_isotopes%in%adduct_weights$Adduct==FALSE)]<-0
          }
          
          #for those adduct forms with correlation coefficients greater than the cutoff:
          if(length(which(check_cor>0)==TRUE)>0){
            
            #check whether one of the assigned adducts is in the filter.by list
            hub_mz_list<-which(check_cor>0  & (cur_adducts_with_isotopes%in%filter.by==TRUE))
            
            #if no adducts match the filter.by option, check if an adduct form matches to the valid adducts defined by user
            if(length(hub_mz_list)<1){
              hub_mz_list<-which(check_cor>0 & check2<0 & (cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))
            }
            
            #replicate of the above procedure
            # if(length(hub_mz_list)<1){
            #   hub_mz_list<-which(check_cor>0 & check2<0 & (cur_adducts_with_isotopes%in%adduct_weights$Adduct==TRUE))
            # }
            
            #if adduct is not in the adduct_weights list, just take the mz
            if(length(hub_mz_list)<1){
              hub_mz_list<-which(check_cor>0 & check2<0)
            }
            
            #which of the correlated features has the highest intensity
            hub_mz_int<-hub_mz_list[which(mchemicaldata$AvgIntensity[hub_mz_list]==max(mchemicaldata$AvgIntensity[hub_mz_list]))[1]]
            
            #of the annotateed features, which has the greatest number of neighbours (based on retention time)
            max_time_neighbors<-0
            best_hub_time_mz<-hub_mz_int

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
            
            mz_name_hub<-paste(mchemicaldata$mz[hub_mz],mchemicaldata$time[hub_mz],sep="_")
            
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
            
            #object for storing information on which features are associated ONLY BY CORRELATIONS (i.e. not by rt and correlation of intensities)
            second_level_associations<-{}

            for(l in layer_one_associations){
              second_level_associations<-c(which(cor_mz[l,]>=0.7))
            }

            
            selected_mz<-c(hub_mz,layer_one_associations) #intersect(layer_one_associations,second_level_associations))
            selected_mz<-unique(selected_mz)

            #matrix showing features correlated with the hub_mz feature (i.e. the one with highest intensity and with greater number of tr neighbours)
            sub_cor_mz<-cor_mz[hub_mz,which(cor_mz[hub_mz,]>=corthresh)]
            
            if(is.na(topquant_cor)==TRUE){
              topquant_cor=0
            }
            
            
            #subset the mchemicaldata object to keep only those features that are associated with one another!
            if(length(which(check_cor2>=0.1))>0){
              mchemicaldata<-mchemicaldata[selected_mz,]		#[which(cor_mz[hub_mz,]>=corthresh & check_cor2>=0.1),]
            }else{
              mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh),]
            }
            
            #remove empty rows 
            
            #mchemicaldata<-na.omit(mchemicaldata) #original version
            keep.indx = which(complete.cases(mchemicaldata[,c(1:(ncol(mchemicaldata)-1))])) #check which rows are not empty
            mchemicaldata = mchemicaldata[keep.indx,]
            rm(keep.indx)
            #if fewer-than 2 rows in mchemicaldata, no need to consider feature further as nothing is associated with it
            if(length(mchemicaldata$mz)<2){
              next
            }
            
            diff_rt<-abs(min(as.numeric(mchemicaldata$time))-max(as.numeric(mchemicaldata$time)))
            
            diff_rt<-round(diff_rt)

            #remove rows from mchemicaldata if no time data is recorded (surely this should have happened at start of processing?)
            if(length(which(is.na(mchemicaldata$time))==TRUE)>0){
              mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$time)==TRUE),]
            }
            
            if(length(mchemicaldata$time)<2){
              next
            }
            
            #if retention times of associated features differ by less than user-defined max_diff_rt value
            if(diff_rt<=max_diff_rt){
              
              #store information on duplicate rows
              dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
              if(length(dup_add)>0){
                dup_data<-mchemicaldata[c(dup_add),]
                #mchemicaldata<-mchemicaldata[-c(dup_add),]
              }
              
              dup_mz_check<-dup_add
              
              if(dim(mchemicaldata)[1]>1){
                #chemical_score<-3
                
                k_power<-1
                
                #associated adducts
                cur_adducts_with_isotopes<-mchemicaldata$Adduct
                
                #fix formating of adduct for isotope (usually looks like M+H_[+1] and becomes M+H)
                cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                
                #determine how many adducts are valid according to user-specified adduct_weights information
                good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
                
                #determine the score for this feature, based on the number of associated adducts
                chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                
                #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
                #print(chemical_score)
                
                if(good_adducts_len>0){
                  
                  chemical_score<-sum(chemical_score*(10^max(adduct_weights[which(adduct_weights$Adduct%in%cur_adducts_with_isotopes),2])))
                  chemical_score<-chemical_score[1]
                }
                
                #### IMPORTANT STEP THAT FACTORS IN ISOTOPE INFORMATION, ON TOP OF THE ADDUCT-BASED SCORE
                # here, increase the chemical score if the +1 or +2 isotopes were found
                check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]\\])")
                if(length(which(check2>0))>0){
                  chemical_score<-100*chemical_score
                }
                
                names(chemical_score)<-chemicalid[1]
                
                #print(chemical_score)
                
                if(length(which(mchemicaldata$Adduct%in%adduct_weights$Adduct))>0){
                  fname<-paste(chemicalid,"score.txt",sep="_")
                }
              }
              
              #return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
            } else{

              ############ 11 April 2019 ######

              #print("fixing time")
              
              #return(mchemicaldata)
              
              #mchemicaldata<-chem_score
              mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
              #mchemicaldata<-cbind(mchemicaldata[,c(2:11)],mchemicaldata[,1],mchemicaldata[,c(12:14)])
              #colnames(mchemicaldata)<-c("mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","AvgIntensity", "MD")
              #updated to include isotope_elements column
              mchemicaldata<-cbind(mchemicaldata[,c(2:11)],mchemicaldata[,1],mchemicaldata[,c(12:15)])
              colnames(mchemicaldata)<-c("mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","AvgIntensity", "MD", "isotope_elements")
              mchemicaldata<-as.data.frame(mchemicaldata)
              mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
              
              #get the identifier of the RTclust module
              groupnumA<-unique(mchemicaldata$Module_RTclust)
              
              #determine the density (and plot) of a given compound along chromatographic axis
              mchemicaldata<-group_by_rt_histv2(mchemicaldata,time_step=1,max_diff_rt=max_diff_rt,groupnum=groupnumA)
              
              #determine frequency of occurrence of each RTclust modules in mchemicaldata obj. and assign names to each in resulting table
              top_mod_sub<-table(mchemicaldata$Module_RTclust)
              top_mod_sub_names<-names(top_mod_sub)
              
              #determine which RTclust module occurred most frequently and subset mchemicaldata to keep only this, then order by mz
              max_top_mod<-which(top_mod_sub==max(top_mod_sub))[1]
              mchemicaldata<-mchemicaldata[which(mchemicaldata$Module_RTclust==top_mod_sub_names[max_top_mod]),]
              mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
              
              #determine which adducts were detected (includes isotopes assigned as adduct, e.g. M+H_[+1] which becomes M+H)
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")

              ## it looks as though interquartile range is calculated in two ways, therefore commented-out first approach
              
              ##determine the interquartile range of retention times within the RTclust module
              #d1<-density(mchemicaldata$time,bw=max_diff_rt,from=min(mchemicaldata$time)-0.001,to=(0.01+max(mchemicaldata$time)),na.rm=TRUE)
              #s1<-summary(d1$x) #mchemicaldata$time)
              
              ##iqr1<-s1[5]-s1[2] #previous approach
              #iqr1<-s1[which(names(s1)=='3rd Qu.')]-s1[which(names(s1)=='1st Qu.')]
              
              ##determine whether the interquartile range was more than half the allowed difference in retention times
              #if(iqr1>max_diff_rt/2){
              #  iqr1<-max_diff_rt/2
              #}
    
              #min_val<-s1[2] #-(1.5*iqr1)
              #max_val<-s1[5] #+(1.5*iqr1)
                
              ######### HERE: adjust retention time tolerance if features within RTclust module are separated by more-than max_diff_rt variable
              print(paste('working on Module_RTclust_', unique(mchemicaldata$Module_RTclust), sep =''))
              print(paste('max_diff_rt was set at: ', max_diff_rt, sep = ''))
              print(paste('attempting to adjust...'))
              s1<-summary(mchemicaldata$time)
              iqr1<-s1[5]-s1[2]
              #print(s1)
              min_val<-s1[2]-(1.5*iqr1)
              max_val<-s1[5]+(1.5*iqr1)

              if(min_val<s1[1] && max_val>s1[6]){
                iqr1<-min(abs(s1[3]-s1[2]),abs(s1[3]-s1[5]))
                min_val<-s1[2]-(1.5*iqr1)
                max_val<-s1[5]+(1.5*iqr1)
                #max_val<-s1[6]
              }
              
              #print(mchemicaldata[hub_mz_ind,])
              if(min_val<s1[1]){
                min_val<-s1[1]
              }
              if(max_val>s1[6]){
                max_val<-s1[6]
              }
              
              iqr1<-min(abs(s1[3]-s1[2]),abs(s1[3]-s1[5]))
              iqr1<-max(iqr1,max_diff_rt)
              
              print(paste('New max_diff_rt value set at: ', iqr1, sep = ''))
              
              #determine range of times for features in RTclust
              diff_rt<-abs(max(mchemicaldata$time)-min(mchemicaldata$time))
              
              if(nrow(mchemicaldata)<1){
                next
              }
              
              #divide data in to time chunks, with the bin size being the iqr1 value, and the range being from below min and above max
              if(diff_rt>2*max_diff_rt){
                time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                  split(x,cut(mchemicaldata$time,breaks=seq(min_val-max_diff_rt,max_val+max_diff_rt,iqr1)))
                })
              }else{
                max_val<-max(mchemicaldata$time)
                min_val<-min(mchemicaldata$time)
                diff_rt<-abs(max(mchemicaldata$time)-min(mchemicaldata$time))
              }
              
              
              if(min_val<max_diff_rt){
                #time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x) split(x,cut(mchemicaldata$time,breaks=seq(0,max_val,1*max_diff_rt))))
                time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                  split(x,cut(mchemicaldata$time,breaks=c(0,max_val+1)))
                })
              }else if(diff_rt<max_diff_rt){
                time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                  split(x,cut(mchemicaldata$time,breaks=c(min_val-diff_rt,max_val+diff_rt,diff_rt*2)))
                })
              }else{
                  time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x){
                    split(x,cut(mchemicaldata$time,breaks=seq(min_val-diff_rt,max_val+diff_rt,1*max_diff_rt)))
                  })
              }

              #determine how many features are in each retention time bin
              group_sizes<-sapply(time_cor_groups,function(x){dim(as.data.frame(x))[1]})
              
              #set up objects for capturing processing outputs
              group_ind_val<-1
              group_ind_size<-1
              check_reladd<-{}
              check_data<-{}
              
              #determine which group had the greatest number of features in it and record the index of that group
              if(length(which(group_sizes>1))>0){
                group_ind_val<-which(group_sizes==max(group_sizes)[1])
              }
              
              #record the number of features in the largest group
              group_ind_size<-max(group_sizes)[1]
              
              #if fewer than 2 features in the largest group
              if(group_ind_size<2){
                
                #as before, now process chemical score
                
                k_power<-1.25
                
                mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
                mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
                  
                #generate list of adducts, including those associated with isotopes (corrected e.g. M+H_[+1] becomes M+H)
                cur_adducts_with_isotopes<-mchemicaldata$Adduct
                cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                  
                good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                # chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))*(1/((diff_rt*0.1)+1)^k_power)
                
                #change
                chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                if(good_adducts_len>0){
                  chemical_score<-sum(chemical_score*(as.numeric(adduct_weights[which(adduct_weights[,1]%in%cur_adducts),2])))
                }
                
              }
                
              good_temp<-{}
              
              for(g1 in 1:length(group_sizes)){
                tempdata<-time_cor_groups[[g1]]
                check_reladd<-which(tempdata$Adduct%in%as.character(adduct_weights$Adduct))
                if(length(check_reladd)>0){
                  good_temp<-c(good_temp,g1)
                  if(nrow(tempdata)>group_ind_size){
                    group_ind_val<-g1
                    group_ind_size<-nrow(tempdata)
                  }
                }
              }
            
              if(length(which(group_sizes>1))>0){
                
                temp_best_score<-(-100)
                temp_best_data<-{}
                
                for(g2 in 1:length(group_sizes)){
                  
                  if(g2%in%good_temp){
                    
                    mchemicaldata<-{}
                    mchemicaldata<-rbind(mchemicaldata,time_cor_groups[[g2]])
                    mchemicaldata<-as.data.frame(mchemicaldata)
                    
                    cur_adducts_with_isotopes<-mchemicaldata$Adduct
                    cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                    
                    if(length(mchemicaldata)<1){
                      next
                    }
                    
                    if(nrow(mchemicaldata)<1){
                      next
                    }

                    diff_rt<-abs(min(as.numeric(mchemicaldata$time))-max(as.numeric(mchemicaldata$time)))
                    diff_rt<-round(diff_rt)

                    dup_add<-{} #which(duplicated(mchemicaldata$Adduct)==TRUE)
                    if(length(dup_add)>0){
                      dup_data<-mchemicaldata[c(dup_add),]
                      #mchemicaldata<-mchemicaldata[-c(dup_add),]
                    }
                    
                    if(length(mchemicaldata)<1){
                      next
                    }
                    if(nrow(mchemicaldata)<1){
                      next
                    }
                    
                    if(diff_rt<=max_diff_rt){

                      if(dim(mchemicaldata)[1]>1){
                        k_power<-1
                        
                        mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
                        cur_adducts_with_isotopes<-mchemicaldata$Adduct
                        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                        
                        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                        #print("score 3")
                        #chemical_score<-2*2^length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                        
                        #change
                        chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                        #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
                        if(good_adducts_len>0){
                          #chemical_score<-chemical_score*(2^adduct_weights[which(adduct_weights[,1]%in%cur_adducts),2])
                          chemical_score<-sum(chemical_score*(10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))))
                          chemical_score<-chemical_score[1]
                        }
                        
                        check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])")
                        
                        if(length(which(check2>0))>0){
                          chemical_score<-100*chemical_score
                        }

                      }else{
                        chemical_score<-0
                      }
                      
                      names(chemical_score)<-chemicalid[1]
                      
                    }else{
                      
                      d1<-density(mchemicaldata$time,bw=max_diff_rt,from=min(mchemicaldata$time)-0.001,to=(0.01+max(mchemicaldata$time)),na.rm=TRUE)
                      s1<-summary(d1$x)
                      iqr1<-s1[5]-s1[2]
                      
                      if(iqr1>max_diff_rt/2){
                        iqr1<-max_diff_rt/2
                      }
                      
                      min_val<-s1[2] #-(1.5*iqr1)
                      max_val<-s1[5] #+(1.5*iqr1)

                      if(length(which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val))>1){
                        mchemicaldata<-mchemicaldata[which(mchemicaldata$time>=(min_val-1) & mchemicaldata$time<=(max_val-1)),]
                        cur_adducts_with_isotopes<-mchemicaldata$Adduct
                        
                        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                        if(dim(mchemicaldata)[1]>1){
                          k_power<-1
                          mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
                          cur_adducts_with_isotopes<-mchemicaldata$Adduct
                          cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                          
                          # print("Setting score 3")
                          #print(mchemicaldata)
                          # chemical_score<-100*2^length(which(cur_adducts%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power) #*length(unique(cur_adducts))
                          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                          
                          #print("score 4")
                          #chemical_score<-2*2^length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                          
                          #change
                          chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                          
                          #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
                          if(good_adducts_len>0){
                            
                            #chemical_score<-chemical_score*(2^adduct_weights[which(adduct_weights[,1]%in%cur_adducts),2])
                            chemical_score<-sum(chemical_score*(10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))))
                            chemical_score<-chemical_score[1]
                          }
                          cur_adducts_with_isotopes<-mchemicaldata$Adduct
                          cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                          check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]\\])")
                          
                          if(length(which(check2>0))>0){
                            chemical_score<-100*chemical_score
                          }

                          #*(1/((diff_rt*0.1)+1)^k_power))
                        }
                        
                        names(chemical_score)<-chemicalid[1]
                        
                      }else{
                        
                        # print("Setting score 4")
                        cur_adducts_with_isotopes<-mchemicaldata$Adduct
                        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                        
                        #print("score 5")
                        #chemical_score<-0
                        #					chemical_score<-length(which(cur_adducts%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power)
                        chemical_score<-length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights$Adduct))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))
                        
                        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
                        
                        #change
                        chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
                        
                        if(good_adducts_len>0){
                          chemical_score<-sum(chemical_score*(as.numeric(adduct_weights[which(adduct_weights$Adduct%in%cur_adducts),2])))
                        }
                        names(chemical_score)<-chemicalid[1]
                        
                      }
                    }
                    
                    if(chemical_score>temp_best_score){
                      temp_best_data<-mchemicaldata
                      mchemicaldata<-temp_best_data
                      temp_best_score<-chemical_score
                    }
                  }
                }
                mchemicaldata<-temp_best_data
                chemical_score<-temp_best_score
              }
              else{
                
                #v1.4.1 removed duplicate
                #dup_mz_check<-which(duplicated(mchemicaldata$Adduct)==TRUE)
                
                #		if(length(dup_mz_check)>0){
                #			mchemicaldata<-mchemicaldata[-c(dup_mz_check),]
                
                #chemical_score<-0
                #chemical_score<-length(which(cur_adducts%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power)
                #names(chemical_score)<-chemicalid[1]
                #					}
                if(dim(mchemicaldata)[1]>1){
                  
                  diff_rt<-abs(min(as.numeric(mchemicaldata$time))-max(as.numeric(mchemicaldata$time)))
                  
                  
                  
                  d1<-density(mchemicaldata$time,bw=max_diff_rt,from=min(mchemicaldata$time)-0.001,to=(0.01+max(mchemicaldata$time)),na.rm=TRUE)
                  s1<-summary(d1$x) #mchemicaldata$time)
                  iqr1<-s1[5]-s1[2]
                  
                  if(iqr1>max_diff_rt/2){
                    iqr1<-max_diff_rt/2
                  }
                  min_val<-s1[2] #-(1.5*iqr1)
                  max_val<-s1[5] #+(1.5*iqr1)
                  
                  if(length(which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val))>1){
                    mchemicaldata<-mchemicaldata[which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val),]
                    dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
                    if(length(dup_add)>0){
                      dup_data<-mchemicaldata[c(dup_add),]
                      #mchemicaldata<-mchemicaldata[-c(dup_add),]
                    }
                    
                    if(dim(mchemicaldata)[1]>1){

                      k_power<-1
                      mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
                      dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
                      
                      if(length(dup_add)>0){
                        dup_data<-mchemicaldata[c(dup_add),]
                        #mchemicaldata<-mchemicaldata[-c(dup_add),]
                      }
                      if(length(mchemicaldata)>0){
                        if(nrow(mchemicaldata)>1){
                          cur_adducts_with_isotopes<-mchemicaldata$Adduct
                          
                          cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                          
                          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                          
                          #print("score 6")
                          # chemical_score<-2*2^length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                          #change
                          chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                          
                          #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
                          if(good_adducts_len>0){
                            chemical_score<-sum(chemical_score*(10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))))
                            chemical_score<-chemical_score[1]
                          }
                          
                          cur_adducts_with_isotopes<-mchemicaldata$Adduct
                          cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                          check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])")
                          
                          if(length(which(check2>0))>0){
                            chemical_score<-100*chemical_score
                          }
                          
                        } else{
                          chemical_score<-0
                        }
                      } else{
                        chemical_score<-0
                      }
                    }
                    
                    names(chemical_score)<-chemicalid[1]
                    
                  } else{
                    dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
                    if(length(dup_add)>0){
                      dup_data<-mchemicaldata[c(dup_add),]
                      #mchemicaldata<-mchemicaldata[-c(dup_add),]
                    }
                    
                    if(length(mchemicaldata)>0){
                      if(nrow(mchemicaldata)>1){
                        cur_adducts_with_isotopes<-mchemicaldata$Adduct
                        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                        if(diff_rt>max_diff_rt){
                          k_power<-3
                          #	print("score 7")
                          #chemical_score<-0
                          #length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights[,1]))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))
                          chemical_score<-length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights[,1]))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))
                          
                          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                          
                          #change
                          chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                          
                          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                          
                          if(good_adducts_len>0){
                            chemical_score<-sum(chemical_score*(as.numeric(adduct_weights[which(adduct_weights[,1]%in%cur_adducts),2])))
                          }
                          
                        }
                      }else{
                        chemical_score<-0
                        #chemical_score<-length(which(cur_adducts%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power)
                      }
                    }else{
                      chemical_score<-0
                    }
                    
                  }
                  #chemical_score<-(topquant_cor)*(1/(diff_rt+1)^k_power)*length(unique(mchemicaldata$Adduct))
                  
                  #print("setting score 6")
                  
                  
                  
                }
                names(chemical_score)<-chemicalid[1]
                
              }
              names(chemical_score)<-chemicalid[1]
              mchemicaldata<-mchemicaldata #[,c(1:11)]
              mchemicaldata<-na.omit(mchemicaldata)
              
              if(length(which(mchemicaldata$Adduct%in%adduct_weights[,1]))>0){
                #   return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
              }
              
            }
          }
          else{
            #no correlation between putative adducts
            #print("setting score 7")
            cur_adducts_with_isotopes<-mchemicaldata$Adduct
            good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights$Adduct))
            #print(mchemicaldata)
            if(good_adducts_len>0){
              #acceptable adducts found
              
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
              
              
              max_adduct_weight<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))[1]
              chemical_score<-((10^max_adduct_weight))
              #print("score 7")
              chemical_score<-chemical_score[1]
              good_adduct_index<-which(adduct_weights[,2]==max_adduct_weight)
              chemical_score<-chemical_score[1]
              mchemicaldata<-mchemicaldata[which(cur_adducts_with_isotopes%in%adduct_weights[good_adduct_index,1]),]
              
              if(chemical_score>best_chemical_score){
                
                best_chemical_score<-chemical_score
                best_conf_level<-1
                best_mod_ind<-i
                best_data<-mchemicaldata
                
              }

            }else{
              #no good adducts found
              
              chemical_score<-0
              mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
              mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
              
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
              #next;
            }
            
          }
          
        }
        
        cur_adducts_with_isotopes<-mchemicaldata$Adduct
        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
        
        check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])")

        if(length(check2)>0){
          
          for(a1 in 1:length(check2)){
            strlength<-attr(check2[[a1]],"match.length")
            
            if(strlength[1]>(-1)){
              
              count_abundant_form<-length(which(cur_adducts%in%cur_adducts[a1]))

              if(count_abundant_form<2){
                chemscoremat_conf_levels<-"None"
                mchemicaldata<-mchemicaldata[-a1,]
              }
            }
          }
          
          
        }
        #if(nrow(mchemicaldata)<1){

        conf_level<-0
        if(nrow(mchemicaldata)<1){
          conf_level<-0
        }else if(nrow(mchemicaldata)>0){
            conf_level<-get_confidence_stage2(curdata=mchemicaldata,adduct_weights=adduct_weights)
            conf_level<-as.numeric(as.character(conf_level))
        } else{
            conf_level<-0
        }
        
        if(nrow(mchemicaldata)>0){

          if(nrow(mchemicaldata)>1){
            
            diff_rt<-max(mchemicaldata$time)-min(mchemicaldata$time)
            k_power=1
            
            if(diff_rt>max_diff_rt){
              k_power=10
            }

            chemical_score<-chemical_score*(1/((diff_rt*0.1)+1)^k_power)
          }
          
        } else{
          chemical_score<-0
        }

        min_chemical_score<-100*2*(1*(corthresh))*(1/((max_diff_rt*0.1)+1)^3)
        
        if(chemical_score>min_chemical_score){
          chemical_score<-chemical_score*(conf_level^conf_level)
        }else{
          chemical_score<-0
        }
        
        if(length(dup_add)>0){
          mchemicaldata<-rbind(mchemicaldata,dup_data)
        }
        
        if(is.na(conf_level)==TRUE){
          conf_level<-0
        }
        
        if(is.na(chemical_score)==TRUE){
          chemical_score<-0
        }
        
        if(chemical_score>best_chemical_score & conf_level>0){  #| conf_level>=best_conf_level)
          
          best_chemical_score<-chemical_score
          best_conf_level<-conf_level
          best_mod_ind<-i
          best_data<-mchemicaldata

        }else (chemical_score==best_chemical_score){
            
            best_chemical_score<-chemical_score
            best_conf_level<-1
            best_mod_ind<-c(i,best_mod_ind)
            best_data<-rbind(best_data,mchemicaldata)
            
        }
      
      print("i is")
      print(i)
      print(mchemicaldata)
      print(top_mod[i])
      print("score is")
      print(chemical_score)
      print(best_chemical_score)
      print("conf level")
      print(conf_level)
      
      
      
      ##for loop complete
      
      if(best_chemical_score>0){
        chemical_score<-best_chemical_score
        best_mod<-best_mod_ind
        mchemicaldata<-best_data
        names(chemical_score)<-chemicalid[1]
        #return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
      }else{
        chemical_score<-0
      }
    }
    #######add code for only correlation criteria here

    
    if(chemical_score<=1){
      
      # || length(which(mchemicaldata$Adduct%in%as.character(filter.by)))<1)
      
      mchemicaldata<-mchemicaldata_orig
      cur_adducts_with_isotopes<-mchemicaldata$Adduct
      cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
      
      good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
      
      if(good_adducts_len>0){
        
        max_adduct_weight<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))[1]
        chemical_score<-((10^max_adduct_weight))
        chemical_score<-chemical_score[1]-1
        good_adduct_index<-which(adduct_weights[,2]==max_adduct_weight)
        chemical_score<-chemical_score[1]
        mchemicaldata<-mchemicaldata[which(cur_adducts_with_isotopes%in%adduct_weights[good_adduct_index,1]),]
      }
    }

    
    if(nrow(mchemicaldata)>0){
      mchemicaldata<-unique(mchemicaldata)
      mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
      
      #dweight1<-degree_weights[which(mzid%in%mzid_cur),]
    }else{
      dweight1<-c(0)
      chemical_score<-0
    }
    
    rm("mzid","global_cor","temp_global_cor")
    return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
  }
  }
  