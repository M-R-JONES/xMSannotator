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
                                       iso_int_tol = 0.05){
  
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
      
      rel_adduct_data<-mchemicaldata[which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct)),]
      rel_adduct_module<-gsub(rel_adduct_data$Module_RTclust,pattern="_[0-9]*",replacement="")
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
    final_isp_annot_res_all<-mchemicaldata
    level_module_isop_annot$Module_RTclust<-gsub(level_module_isop_annot$Module_RTclust,pattern="_[0-9]*",replacement="")
    
    mchemicaldata_goodadducts_index<-which(mchemicaldata$Adduct%in%as.character(adduct_weights$Adduct))
    final_isp_annot_res_isp<-{}
    
    if(length(mchemicaldata_goodadducts_index)>0){
      
      final_isp_annot_res_isp<-lapply(1:length(mchemicaldata_goodadducts_index),function(i){
        
        m<-mchemicaldata_goodadducts_index[i]
        final_isp_annot_res<-cbind(paste("group",i,sep=""),mchemicaldata[m,])
        isp_group<-as.character(mchemicaldata$ISgroup[m])
        module_rt_group<-as.character(mchemicaldata$Module_RTclust[m])
        module_rt_group<-gsub(module_rt_group,pattern="_[0-9]*",replacement="")
        
        ##########################################################################################
        #from the formula, calculate isotope profile
        
        ##### IMPORTANT - NEED TO FIND A WAY OF GENERATING ISOTOPIC INFORMATION FOR DIMERS / TRIMERS etc
        
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
        
        #correct formula to take in to account adduct elements
        if(!is.na(adduct_weights[indx_add_table,]$Merge_add)){
          form = enviPat::mergeform(form, adduct_weights[indx_add_table,]$Merge_add)
        }
        
        if(!is.na(adduct_weights[indx_add_table,]$Merge_sub)){
          form = enviPat::subform(form, adduct_weights[indx_add_table,]$Merge_sub)
        }
        
        print('Sum formula for compound and associated adduct is: ', form, sep = '')
        
        ##########################################################################################
        
        #from enviPat package, read in isotopes table
        data(isotopes)
        
        #calculate the fine isotope structure for the formula, taking in to account charge-state and adduct elements
        mol_isos = enviPat::isopattern(isotopes, chemforms = form, charge = z, threshold = 0.1)
        
        temp = mol_isos[[1]][,1]
        
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
          #multiple the compound intensity by the ratio of its intensity to that of its largest isotope
          ratio = max(top_isotopes$abundance) / mchemicaldata$AvgIntensity[m]
          query_int = ratio * mchemdata$AvgIntensity[m]
          #this is a tolerance factor that adds 10% further to the intensity value to ensure largest isotope captured
          query_int = query int + (0.1 * query_int)
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
        
        if(length(put_isp_masses_curmz_data)>0){
          
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
          abund_ratio_vec = top_isotopes$abundance
          ischeck<-which(int_vec<=max(abund_ratio_vec+0.20))
          
          put_isp_masses_curmz_data[,2]<-as.numeric(as.character(put_isp_masses_curmz_data[,2]))
          
          if(length(ischeck)>0){
            
            for(rnum in 1:length(ischeck)){
              temp_var<-{}
              bool_check<-1
              
              isp_v<-ischeck[rnum]
              
              isp_v<-as.numeric(as.character(isp_v))
              #print(isp_v)
              #print(put_isp_masses_curmz_data[isp_v,])
              
              # print(mchemicaldata[m,2])
              diff_rt<-abs(put_isp_masses_curmz_data[isp_v,2]-mchemicaldata[m,2])
              
              #### MODIFICATION
              
              isnum = round((put_isp_masses_curmz_data[isp_v,1]-mchemicaldata$mz[m])*z)
              
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
              
              if(diff_rt<max_diff_rt & isnum<=max_isp){
                
                max_isp_count=max_isp
                
                if(max_isp_count>0 && isnum<=max_isp_count && bool_check>0){
                  
                  
                  ###### CHANGED HERE!!!!!!#######
                  #isnum2<-(round(put_isp_masses_curmz_data[isp_v,1])-round(mchemicaldata$MonoisotopicMass[m]))
                  #isnum2<-(round(put_isp_masses_curmz_data[isp_v,1])-round(mchemicaldata$mz[m]))
                  #isnum2<-round(isnum2)
                  
                  #ADDED MULTIPLICATION TO ACCOUNT FOR CHARGE STATE OF ADDUCT
                  #doubly-charged ions that are 1 Da apart are +2 isotope.
                  isnum2 = round((put_isp_masses_curmz_data[isp_v,1]-mchemicaldata$mz[m])*z)
                  
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
                    iso_hits = ((put_isp_masses_curmz_data[isp_v,1] - top_isotopes$`m/z`) / top_isotopes$`m/z`)* 1000000
                    iso_hits = top_isotopes[which(iso_hits < iso_ppm_tol),]
                    
                    #IF > 1 ISOTOPE WITHIN PPM TOLERANCE, CHECK IF ANY WITHIN ISO_INT_TOL
                    #start by calculating ratio of candidate isotope peak to monoisotopic peak
                    percent = (put_isp_masses_curmz_data[isp_v,]$AvgIntensity / mchemicaldata[m,]$mean_int_vec)*100
                    
                    ##### IMPORTANT PARAMETER ###########
                    #iso_int_tol = 0.05 #5% tolerance allowed
                    
                    #collect isotopes that occur within +/- the intensity tolerance
                    iso_hits = subset(iso_hits, 
                                      iso_hits$abundance > (percent - (percent * iso_int_tol)) &
                                        iso_hits$abundance < (percent + (percent * iso_int_tol)))
                    
                    if(nrow(iso_hits) > 1){
                      form_name_alt = paste('_[', paste(iso_hits$name, collapse='//'), ']', sep='')
                      form_name = form_name_alt
                    }else{
                      form_name = as.character(paste(mchemicaldata[m,7],"_[",isp_sign,(isnum2),"]",sep=""))
                    }
                    
                    #COMMENTED THIS LINE!!!!
                    #form_name<-as.character(paste(mchemicaldata[m,7],"_[",isp_sign,(isnum2),"]",sep=""))
                    
                    #end of modification
                    ###########################################################
                    
                    other_inf<-cbind(rep("-",7))
                    temp_var<-cbind(put_isp_masses_curmz_data[isp_v,c(1:2)],t(other_inf),put_isp_masses_curmz_data[isp_v,c(3:4)],put_isp_masses_curmz_data[isp_v,c(2,5:6)])
                    
                    temp_var<-as.data.frame(temp_var)
                    
                    colnames(temp_var)<-colnames(mchemicaldata)
                    
                    temp_var$Formula<-form_name
                    temp_var$Name<-as.character(mchemicaldata[m,6])
                    temp_var$chemical_ID<-as.character(mchemicaldata[m,5])
                    
                    #temp_var$Adduct<-paste(mchemicaldata[m,9],"_[+",isnum,"]",sep="")
                    
                    adductname=mchemicaldata[m,9] #queryadductlist[adnum]
                    #adductmass=adduct_table[as.character(adductname),4]
                    
                    #temp_var$Adduct<-paste("M","_[",isp_sign,(abs(isnum2)),"]",sep="")
                    
                    ###### MODIFIED THIS - incorrectly reference 'isnum' versus 'isnum2'
                    temp_var$Adduct<-paste(mchemicaldata[m,9],"_[",isp_sign,(abs(isnum2)),"]",sep="")
                    
                    temp_var<-as.data.frame(temp_var)
                    temp_var<-cbind(paste("group",i,sep=""),temp_var)
                    final_isp_annot_res<-as.data.frame(final_isp_annot_res)
                    
                    if(nrow(temp_var)>0){
                      
                      check_mz<-which(temp_var$mz%in%final_isp_annot_res)
                      
                      if(length(check_mz)>0){
                        temp_var<-temp_var[-c(check_mz),]
                      }
                      
                      if(nrow(temp_var)>0){
                        final_isp_annot_res<-rbind(final_isp_annot_res,temp_var)
                      }
                      
                      
                      
                    }
                  }
                }
              }
            }
          }
        }
        
        return(final_isp_annot_res)
        
      })
      
      
      rm(level_module_isop_annot) 
      final_isp_annot_res2<-ldply(final_isp_annot_res_isp,rbind)
      
      
      isp_group_check<-table(final_isp_annot_res2[,1])
      good_groups<-which(isp_group_check==max(isp_group_check))
      
      group_name<-names(isp_group_check)[good_groups]
      final_isp_annot_res2<-as.data.frame(final_isp_annot_res2)
      
      final_isp_annot_res2<-final_isp_annot_res2 #[which(final_isp_annot_res2[,1]%in%group_name),]
      
      final_isp_annot_res2<-final_isp_annot_res2[,-c(1)]
      final_isp_annot_res2<-as.data.frame(final_isp_annot_res2)
      
      rm(final_isp_annot_res_isp)
      
      
      #write.table(temp_var,file="finaltempvar.txt",sep="\t",row.names=FALSE)
      
      mchemicaldata<-rbind(final_isp_annot_res_all,final_isp_annot_res2)  #[,-c(12)]
      
    }
    mchemicaldata<-unique(mchemicaldata)
    
    
    
    bad_rows<-which(is.na(mchemicaldata$mz)==TRUE)
    if(length(bad_rows)>0){
      mchemicaldata<-mchemicaldata[-c(bad_rows),]
    }
    
    mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
    
    
    mod_names<-mchemicaldata$Module_RTclust
    mod_names<-unique(mod_names)
    
    
    
    mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")
    
    
    
    temp_global_cor<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
    
    #save(temp_global_cor,file="temp_global_cor.Rda")
    #rm(global_cor) 
    
    #rm(mzid)
    
    #mzid<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
    
    #if(FALSE)
    {
      diffmatB<-{}
      diffmatB<-lapply(1:length(mod_names),function(i){
        
        groupA_num<-mod_names[i]
        
        subdata<-mchemicaldata[which(mchemicaldata$Module_RTclust==groupA_num),]
        subdata<-subdata[order(subdata$time),]
        
        
        
        if(nrow(subdata)>0){
          
          #print(subdata)
          #print(length(subdata))
          groupB<-group_by_rt_histv2(subdata,time_step=1,max_diff_rt=10,groupnum=groupA_num)
          
        }else{
          groupB<-subdata
        }
        rownames(groupB)<-NULL
        
        
        #print(groupB)
        return(groupB)
      })
      
      
      
      mchemicaldata<-ldply(diffmatB,rbind)
    }
    mchemicaldata<-unique(mchemicaldata)
    
    rm(diffmatB)
    dupmz<-{} #which(duplicated(mchemicaldata$mz)==TRUE)
    
    if(length(dupmz)>0){
      mchemicaldata<-mchemicaldata[-c(dupmz),]
      
    }
    
    rm(final_isp_annot_res)
    
    
    
    write.table(mchemicaldata,file="../Stage2_withisotopes.txt",append=TRUE,sep="\t",col.names=FALSE)
    
    #write.csv(mchemicaldata,file="../Stage2_withisotopes.csv",append=TRUE,col.names=FALSE)
    
    table_mod<-table(mchemicaldata$Module_RTclust)
    
    table_iso<-table(mchemicaldata$ISgroup)
    
    table_mod<-table_mod[table_mod>0]
    
    table_mod<-table_mod[order(table_mod,decreasing=TRUE)]
    
    mchemicaldata_orig<-mchemicaldata
    
    top_mod<-names(table_mod)
    
    bool_check<-0
    
    
    
    topquant_cor<-0
    
    best_conf_level<-(-100)
    
    k_power<-1
    if(length(which(table_mod>=1))>0)
    {
      best_chemical_score<-(-100)
      
      for(i in 1:length(which(table_mod>=1)))
      {
        
        
        dup_add<-{}
        
        chemical_score<-(-99999)
        conf_level<-0
        mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
        
        
        mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
        
        cur_adducts_with_isotopes<-mchemicaldata$Adduct
        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
        
        if(nrow(mchemicaldata)<2)
        {
          
          
          good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
          if(good_adducts_len>0){
            
            
            chemical_score<-(1*(10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))))
            chemical_score<-chemical_score[1]
            if(chemical_score>best_chemical_score)
            {
              
              best_chemical_score<-chemical_score
              
              best_conf_level<-1
              
              
              best_mod_ind<-i
              best_data<-mchemicaldata
              
              
            }else{
              
              if(chemical_score==best_chemical_score)
              {
                
                best_chemical_score<-chemical_score
                
                best_conf_level<-1
                
                
                best_mod_ind<-c(i,best_mod_ind)
                best_data<-rbind(best_data,mchemicaldata)
                
                
              }
              
              
              
            }
          }
          #next;
        }
        check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])")
        mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_") #mzid_cur<-paste(curmchemdata$mz,curmchemdata$time,sep="_") #mzid_cur<-paste(chem_score$filtdata$mz,chem_score$filtdata$time,sep="_")
        
        dup_features<-which(duplicated(mzid_cur)==TRUE)
        if(length(dup_features)>0){
          
          mchemicaldata<-mchemicaldata[-c(dup_features),]
          mzid_cur<-paste(mchemicaldata$mz,mchemicaldata$time,sep="_")
        }
        
        
        
        cor_mz<-global_cor[which(mzid%in%mzid_cur),which(mzid%in%mzid_cur)]
        
        
        cor_mz<-round(cor_mz,1)
        
        #return(cor_mz)
        if(length(cor_mz)>1){
          corrownamesA<-rownames(cor_mz)
          
          
          
          mat_rownames<-strsplit(as.character(corrownamesA),split="_")
          
          #  print(mchemicaldata)
          #print(cor_mz)
          m1<-{}
          for(i in 1:length(mat_rownames)){m1<-rbind(m1,cbind(mat_rownames[[i]] [1],mat_rownames[[i]][2]))}
          m1<-as.data.frame(m1)
          colnames(m1)<-c("mz","time")
          cor_mz2<-cbind(m1,cor_mz)
          
          mz_order<-order(cor_mz2$mz)
          cor_mz2<-cor_mz2[c(mz_order),]
          
          cor_mz<-cor_mz2[,-c(1:2)]
          cor_mz<-cor_mz[,mz_order]
        }
        
        # print(corthresh)
        if(length(cor_mz)>1){
          check_cor<-sapply(1:dim(cor_mz)[1],function(k){
            
            count_mz<-length(which(cor_mz[k,]>=corthresh))-1
            return(count_mz)
          })
          
          topquant_cor<-max(cor_mz[upper.tri(cor_mz)])
          
          check_cor2<-check_cor/length(check_cor)
          
          if(length((cur_adducts_with_isotopes%in%adduct_weights[,1]==TRUE))>0){
            check_cor[check2>0 | (cur_adducts_with_isotopes%in%adduct_weights[,1]==FALSE)]<-0
            
          }
          
          
          
          
          #return(cor_mz)
          #at least score of 2
          if(length(which(check_cor>0)==TRUE)>0)
          {
            
            #hub_mz<-which(check_cor==max(check_cor)[1] & (cur_adducts_with_isotopes%in%adduct_weights[,1]==TRUE))
            
            hub_mz_list<-which(check_cor>0  & (cur_adducts_with_isotopes%in%filter.by==TRUE))
            
            if(length(hub_mz_list)<1){
              hub_mz_list<-which(check_cor>0 & check2<0 & (cur_adducts_with_isotopes%in%adduct_weights[,1]==TRUE))
            }
            
            
            if(length(hub_mz_list)<1){
              hub_mz_list<-which(check_cor>0 & check2<0 & (cur_adducts_with_isotopes%in%adduct_weights[,1]==TRUE))
            }
            
            if(length(hub_mz_list)<1){
              hub_mz_list<-which(check_cor>0 & check2<0)
            }
            
            hub_mz_int<-hub_mz_list[which(mchemicaldata$mean_int_vec[hub_mz_list]==max(mchemicaldata$mean_int_vec[hub_mz_list]))[1]]
            
            max_time_neighbors<-0
            best_hub_time_mz<-hub_mz_int
            
            for(h1 in hub_mz_list){
              
              
              mz_name_hub<-paste(mchemicaldata$mz[h1],mchemicaldata$time[h1],sep="_")
              
              hub_rt<-mchemicaldata$time[h1]
              
              
              
              diff_rt_hubmz<-apply(mchemicaldata,1,function(k){curtime<-as.numeric(as.character(k[3]));return(abs(hub_rt-curtime))})
              
              num_time_neighbors<-length(which(diff_rt_hubmz<=max_diff_rt))
              
              
              if(num_time_neighbors>max_time_neighbors){
                
                best_hub_time_mz<-h1
                max_time_neighbors<-num_time_neighbors
              }
            }
            
            
            #print(best_hub_time_mz)
            #		print(mchemicaldata)
            
            hub_mz<-best_hub_time_mz
            
            mz_name_hub<-paste(mchemicaldata$mz[hub_mz],mchemicaldata$time[hub_mz],sep="_")
            
            hub_rt<-mchemicaldata$time[hub_mz]
            
            diff_rt_hubmz<-apply(mchemicaldata,1,function(k){curtime<-as.numeric(as.character(k[3]));return(abs(hub_rt-curtime))})
            #print(diff_rt_hubmz)
            #print(cor_mz)
            #print(hub_mz)
            #print(mchemicaldata$mean_int_vec)
            #print(cor_mz[hub_mz,])
            
            #stop("sf")
            if(MplusH.abundance.ratio.check==TRUE){
              
              layer_one_associations<-which(cor_mz[hub_mz,]>=corthresh & mchemicaldata$mean_int_vec<mchemicaldata$mean_int_vec[hub_mz] & diff_rt_hubmz<=max_diff_rt)
            }else{
              
              layer_one_associations<-which(cor_mz[hub_mz,]>=corthresh & diff_rt_hubmz<=max_diff_rt)
            }
            second_level_associations<-{}
            
            for(l in layer_one_associations){
              second_level_associations<-c(second_level_associations,which(cor_mz[l,]>=0.7))
              
            }
            #print("layer one")
            #  print(layer_one_associations)
            selected_mz<-c(hub_mz,layer_one_associations) #intersect(layer_one_associations,second_level_associations))
            
            selected_mz<-unique(selected_mz)
            
            
            sub_cor_mz<-cor_mz[hub_mz,which(cor_mz[hub_mz,]>=corthresh)]
            
            #topquant_cor<-quantile(sub_cor_mz,0.9,na.rm=TRUE)[1]
            
            if(is.na(topquant_cor)==TRUE){
              topquant_cor=0
            }
            
            
            if(length(which(check_cor2>=0.1))>0){
              
              mchemicaldata<-mchemicaldata[selected_mz,]		#[which(cor_mz[hub_mz,]>=corthresh & check_cor2>=0.1),]
              #print(mchemicaldata)
              
              
            }else{
              
              mchemicaldata<-mchemicaldata[which(cor_mz[hub_mz,]>=corthresh),]
            }
            
            #	print(mchemicaldata)
            #return(mchemicaldata)
            mchemicaldata<-na.omit(mchemicaldata)
            if(nrow(mchemicaldata)<2){
              next
            }
            #print(mchemicaldata)
            #print(nrow(mchemicaldata))
            mchemicaldata<-as.data.frame(mchemicaldata)
            
            diff_rt<-abs(min(as.numeric(mchemicaldata$time))-max(as.numeric(mchemicaldata$time)))
            
            diff_rt<-round(diff_rt)
            
            
            
            
            if(length(which(is.na(mchemicaldata$time))==TRUE)>0){
              mchemicaldata<-mchemicaldata[-which(is.na(mchemicaldata$time)==TRUE),]
            }
            
            if(nrow(mchemicaldata)<2){
              next
            }
            
            
            
            
            
            if(diff_rt<=max_diff_rt)
            {
              
              #print("DOING THIS")
              dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
              if(length(dup_add)>0){
                dup_data<-mchemicaldata[c(dup_add),]
                #mchemicaldata<-mchemicaldata[-c(dup_add),]
              }
              
              dup_mz_check<-dup_add
              #print(dup_add)
              #print(mchemicaldata)
              if(dim(mchemicaldata)[1]>1){
                #chemical_score<-3
                
                k_power<-1
                mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
                
                cur_adducts_with_isotopes<-mchemicaldata$Adduct
                cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
                
                good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                #print("score 2")
                #chemical_score<-2*2^length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                
                #change
                
                chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                
                #chemical_score<-chemical_score/sqrt(log2((0.5*diff_rt)+1))
                #print(chemical_score)
                if(good_adducts_len>0){
                  
                  chemical_score<-sum(chemical_score*(10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))))
                  chemical_score<-chemical_score[1]
                }
                
                
                check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]\\])")
                
                #print(chemical_score)
                if(length(which(check2>0))>0){
                  
                  chemical_score<-100*chemical_score
                }
                names(chemical_score)<-chemicalid[1]
                
                #print(chemical_score)
                
                if(length(which(mchemicaldata$Adduct%in%adduct_weights[,1]))>0){
                  
                  fname<-paste(chemicalid,"score.txt",sep="_")
                  
                  
                }
              }
              
              #return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
            }
            else{
              
              #print("fixing time")
              
              #return(mchemicaldata)
              
              #mchemicaldata<-chem_score
              mchemicaldata$Module_RTclust<-gsub(mchemicaldata$Module_RTclust,pattern="_[0-9]*",replacement="")
              
              mchemicaldata<-cbind(mchemicaldata[,c(2:11)],mchemicaldata[,1],mchemicaldata[,c(12:14)])
              
              colnames(mchemicaldata)<-c("mz","time","MatchCategory","theoretical.mz","chemical_ID","Name","Formula","MonoisotopicMass","Adduct","ISgroup","Module_RTclust","time.y","mean_int_vec", "MD")
              
              #return(mchemicaldata)
              
              mchemicaldata<-as.data.frame(mchemicaldata)
              
              mchemicaldata$time<-as.numeric(as.character(mchemicaldata$time))
              
              #return(mchemicaldata)
              
              groupnumA<-unique(mchemicaldata$Module_RTclust)
              
              mchemicaldata<-group_by_rt_histv2(mchemicaldata,time_step=1,max_diff_rt=max_diff_rt,groupnum=groupnumA)
              
              #print("done")
              top_mod_sub<-table(mchemicaldata$Module_RTclust)
              
              top_mod_sub_names<-names(top_mod_sub)
              
              max_top_mod<-which(top_mod_sub==max(top_mod_sub))[1]
              
              mchemicaldata<-mchemicaldata[which(mchemicaldata$Module_RTclust==top_mod_sub_names[max_top_mod]),]
              
              
              mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
              
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
              
              #print(mchemicaldata)
              
              
              d1<-density(mchemicaldata$time,bw=max_diff_rt,from=min(mchemicaldata$time)-0.001,to=(0.01+max(mchemicaldata$time)),na.rm=TRUE)
              s1<-summary(d1$x) #mchemicaldata$time)
              iqr1<-s1[5]-s1[2]
              
              if(iqr1>max_diff_rt/2){
                iqr1<-max_diff_rt/2
              }
              min_val<-s1[2] #-(1.5*iqr1)
              max_val<-s1[5] #+(1.5*iqr1)
              
              
              {
                
                s1<-summary(mchemicaldata$time)
                iqr1<-s1[5]-s1[2]
                #print(s1)
                min_val<-s1[2]-(1.5*iqr1)
                max_val<-s1[5]+(1.5*iqr1)
                
                #print(min_val)
                #print(max_val)
                if(min_val<s1[1] && max_val>s1[6]){
                  #min_val<-s1[1]
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
                
                
                diff_rt<-abs(max(mchemicaldata$time)-min(mchemicaldata$time))
                
                #	print("fixed time")
                #				print(diff_rt)
                #				print(mchemicaldata)
                #return(mchemicaldata)
                
                if(nrow(mchemicaldata)<1){
                  next
                }
                
                if(diff_rt>2*max_diff_rt){
                  time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x)  split(x,cut(mchemicaldata$time,breaks=seq(min_val-max_diff_rt,max_val+max_diff_rt,iqr1))))
                  
                }else{
                  
                  max_val<-max(mchemicaldata$time)
                  min_val<-min(mchemicaldata$time)
                  diff_rt<-abs(max(mchemicaldata$time)-min(mchemicaldata$time))
                  
                  if(min_val<max_diff_rt){
                    
                    #time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x) split(x,cut(mchemicaldata$time,breaks=seq(0,max_val,1*max_diff_rt))))
                    
                    time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x) split(x,cut(mchemicaldata$time,breaks=c(0,max_val+1))))
                  }else{
                    
                    if(diff_rt<max_diff_rt){
                      
                      time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x) split(x,cut(mchemicaldata$time,breaks=c(min_val-diff_rt,max_val+diff_rt,diff_rt*2))))
                    }else{
                      time_cor_groups<-sapply(list(myData1=mchemicaldata),function(x)  split(x,cut(mchemicaldata$time,breaks=seq(min_val-diff_rt,max_val+diff_rt,1*max_diff_rt))))
                    }
                  }
                }
                
                
                
              }
              
              #return(time_cor_groups)
              
              group_sizes<-sapply(time_cor_groups,function(x){dim(as.data.frame(x))[1]})
              
              
              
              group_ind_val<-1
              group_ind_size<-1
              check_reladd<-{}
              check_data<-{}
              if(length(which(group_sizes>1))>0){
                group_ind_val<-which(group_sizes==max(group_sizes)[1])
                
              }
              group_ind_size<-max(group_sizes)[1]
              
              if(group_ind_size<2){
                
                
                
                k_power<-1.25
                mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
                
                
                mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
                
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
              for(g1 in 1:length(group_sizes))
              {
                tempdata<-time_cor_groups[[g1]]
                
                
                check_reladd<-which(tempdata$Adduct%in%as.character(adduct_weights[,1]))
                if(length(check_reladd)>0){
                  good_temp<-c(good_temp,g1)
                  if(nrow(tempdata)>group_ind_size){
                    group_ind_val<-g1
                    group_ind_size<-nrow(tempdata)
                    
                  }
                  
                  
                  
                }
                
              }
              
              
              if(length(which(group_sizes>1))>0)
              {
                
                
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
                    #print("diff rT")
                    #print(diff_rt)
                    #print(mchemicaldata)
                    
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
                    if(diff_rt<=max_diff_rt)
                    {
                      
                      
                      if(dim(mchemicaldata)[1]>1)
                      {
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
                      
                      
                      if(length(which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val))>1)
                      {
                        mchemicaldata<-mchemicaldata[which(mchemicaldata$time>=(min_val-1) & mchemicaldata$time<=(max_val-1)),]
                        cur_adducts_with_isotopes<-mchemicaldata$Adduct
                        
                        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
                        if(dim(mchemicaldata)[1]>1)
                        {
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
                        chemical_score<-length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights[,1]))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))
                        
                        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                        
                        #change
                        chemical_score<-length(unique(cur_adducts))*good_adducts_len*(1*(topquant_cor))
                        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
                        
                        if(good_adducts_len>0){
                          
                          chemical_score<-sum(chemical_score*(as.numeric(adduct_weights[which(adduct_weights[,1]%in%cur_adducts),2])))
                          
                        }
                        
                        names(chemical_score)<-chemicalid[1]
                        
                      }
                    }
                    
                    if(chemical_score>temp_best_score)
                    {
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
                  
                  if(length(which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val))>1)
                  {
                    mchemicaldata<-mchemicaldata[which(mchemicaldata$time>=min_val & mchemicaldata$time<=max_val),]
                    dup_add<-which(duplicated(mchemicaldata$Adduct)==TRUE)
                    if(length(dup_add)>0){
                      dup_data<-mchemicaldata[c(dup_add),]
                      #mchemicaldata<-mchemicaldata[-c(dup_add),]
                    }
                    
                    if(dim(mchemicaldata)[1]>1)
                    {
                      
                      
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
                          
                        }else{
                          
                          
                          chemical_score<-0
                          
                        }
                      }else{
                        chemical_score<-0
                        
                      }
                      
                      
                    }
                    names(chemical_score)<-chemicalid[1]
                    
                  }else
                  {
                    
                    
                    
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
            good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
            #print(mchemicaldata)
            if(good_adducts_len>0)
            {
              
              
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
              
              
              max_adduct_weight<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))[1]
              chemical_score<-((10^max_adduct_weight))
              #print("score 7")
              chemical_score<-chemical_score[1]
              good_adduct_index<-which(adduct_weights[,2]==max_adduct_weight)
              chemical_score<-chemical_score[1]
              mchemicaldata<-mchemicaldata[which(cur_adducts_with_isotopes%in%adduct_weights[good_adduct_index,1]),]
              
              
              
              if(FALSE){
                if(chemical_score>best_chemical_score)
                {
                  
                  best_chemical_score<-chemical_score
                  
                  best_conf_level<-1
                  
                  
                  best_mod_ind<-i
                  best_data<-mchemicaldata
                  
                  
                  
                  
                }
              }
            }else{
              #
              
              chemical_score<-0
              mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
              
              
              mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]
              
              cur_adducts_with_isotopes<-mchemicaldata$Adduct
              cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
              #next;
            }
            
          }
          
        }
        
        
        #mchemicaldata<-mchemicaldata[order(mchemicaldata$Adduct),]
        
        cur_adducts_with_isotopes<-mchemicaldata$Adduct
        
        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")
        
        check2<-gregexpr(text=cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])")
        
        
        if(length(check2)>0){
          
          for(a1 in 1:length(check2)){
            strlength<-attr(check2[[a1]],"match.length")
            
            if(strlength[1]>(-1)){
              
              count_abundant_form<-length(which(cur_adducts%in%cur_adducts[a1]))
              
              
              if(count_abundant_form<2)
              {
                
                chemscoremat_conf_levels<-"None"
                
                mchemicaldata<-mchemicaldata[-a1,]
                
              }
              
              
            }
          }
          
          
        }
        #if(nrow(mchemicaldata)<1){
        
        
        conf_level<-0
        if(length(mchemicaldata)<1){
          
          conf_level<-0
        }else{
          
          if(nrow(mchemicaldata)>0){
            conf_level<-get_confidence_stage2(curdata=mchemicaldata,adduct_weights=adduct_weights)
            
            
            conf_level<-as.numeric(as.character(conf_level))
          }else{
            
            conf_level<-0
          }
        }
        
        if(length(mchemicaldata)>0){
          if(nrow(mchemicaldata)>1){
            
            diff_rt<-max(mchemicaldata$time)-min(mchemicaldata$time)
            
            k_power=1
            if(diff_rt>max_diff_rt)
            {
              k_power=10
            }
            chemical_score<-chemical_score*(1/((diff_rt*0.1)+1)^k_power)
            
            
            
          }
        }else{
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
        
        if(chemical_score>best_chemical_score & conf_level>0)  #| conf_level>=best_conf_level)
        {
          
          
          
          best_chemical_score<-chemical_score
          
          best_conf_level<-conf_level
          
          
          best_mod_ind<-i
          best_data<-mchemicaldata
          
          
        }else{
          
          
          if(chemical_score==best_chemical_score)
          {
            
            best_chemical_score<-chemical_score
            
            best_conf_level<-1
            
            
            best_mod_ind<-c(i,best_mod_ind)
            best_data<-rbind(best_data,mchemicaldata)
            
            
          }
          
          
          
          
        }
        
        if(FALSE)
        {
          print("i is")
          print(i)
          print(mchemicaldata)
          print(top_mod[i])
          print("score is")
          print(chemical_score)
          print(best_chemical_score)
          print("conf level")
          print(conf_level)
        }
      }
      
      
      
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
    
    
    #if(FALSE)
    {
      #if(chemical_score<=0)
      #print("length is")
      #print(length(which(mchemicaldata$Adduct%in%as.character(filter.by))))
      if(chemical_score<=1) # || length(which(mchemicaldata$Adduct%in%as.character(filter.by)))<1)
      {
        mchemicaldata<-mchemicaldata_orig
        cur_adducts_with_isotopes<-mchemicaldata$Adduct
        cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[1-2]*\\])",replacement="")
        
        good_adducts_len<-length(which(cur_adducts_with_isotopes%in%adduct_weights[,1]))
        
        if(good_adducts_len>0)
        {
          
          max_adduct_weight<-max(as.numeric(as.character(adduct_weights[which(adduct_weights[,1]%in%cur_adducts_with_isotopes),2])))[1]
          chemical_score<-((10^max_adduct_weight))
          chemical_score<-chemical_score[1]-1
          good_adduct_index<-which(adduct_weights[,2]==max_adduct_weight)
          chemical_score<-chemical_score[1]
          mchemicaldata<-mchemicaldata[which(cur_adducts_with_isotopes%in%adduct_weights[good_adduct_index,1]),]
          
          
        }
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
    
    #
    #print(best_chemical_score)
    rm("mzid","global_cor","temp_global_cor")
    
    return(list("chemical_score"=chemical_score,"filtdata"=mchemicaldata))
  }
}