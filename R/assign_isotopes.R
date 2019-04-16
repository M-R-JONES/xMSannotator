#' assign_isotopes
#' 
#' This function expands putative adduct assignments with associated isotopes.
#' 
#' @param add_indx Integer indicating the index of the adduct in the mchemicaldata dataframe
#' @param mchemicaldata Dataframe with columns mz (numeric), time (numeric), MatchCategory (character; assigned as either "Unique"
#' or "Multiple"), theoretical.mz (numeric; monoisotopic mass + adduct mass),
#' chemical_ID (character), Name (character; name of suspect), Formula (character; molecular formula without adducts), 
#' MonoisotopicMass (numeric; monoisotopic mass excluding adducts), Adduct (character; see adduct_weights)
#' @param adduct_weights
#' @param level_module_isop_annot
#' @param max_diff_rt
#' @param mass_defect_window
#' @param mass_defect_mode
#' @param max_isp Integer indicating the top 'n' isotope peaks to consider during assignment
#' @param iso_int_tol Numeric. The valid percent by which the putative isotope feature may differ from theoretical isotope intensity
#' @param iso_ppm_tol Numeric. The tolerance applied along the m/z axis during assigning isotopes to putative isotopic features
#' @import "enviPat"
#' @import "Rdisop"


assign_isotopes = function(add_indx,
                           mchemicaldata,
                           adduct_weights, 
                           level_module_isop_annot,
                           max_diff_rt = 10, 
                           mass_defect_window = 0.1, 
                           mass_defect_mode = "pos",
                           max_isp = 5, 
                           iso_int_tol = 0.3, 
                           iso_ppm_tol = 5,
                           i = ''){
  
  print(paste('Adduct index', add_indx, sep = ''))
  
  #from mchemicaldata, extract the row corresponding to the input adduct index
  #final_isp_annot_res<-cbind(paste("group",i,sep=""),mchemicaldata[add_indx,])#used before parLapply
  final_isp_annot_res<-cbind(paste("group",i,sep=""),mchemicaldata[add_indx,])
  
  #get the ISgroup corresponding to the input adduct
  isp_group<-as.character(mchemicaldata$ISgroup[add_indx])
  
  #get the cluster (based on mz and retention) for currently-considered mz
  module_rt_group<-as.character(mchemicaldata$Module_RTclust[add_indx])
  module_rt_group<-gsub(module_rt_group,pattern="_[0-9]*",replacement="")
  
  ##########################################################################################
  #from the formula, calculate isotope profile
  
  #here, set form to the formula for the assigned suspect
  form = final_isp_annot_res$Formula
  
  #get reference to row in adducts_weights table corresponding to current adduct
  indx_add_table = which(adduct_weights$Adduct == as.character(mchemicaldata[add_indx,]$Adduct))
  print('found index')
  
  #multiply formula weight by number of monomers assigned in adduct
  nummols = as.numeric(adduct_weights[indx_add_table,]$num_molecules)
  
  #get the charge state associated with the adduct form
  z = as.numeric(adduct_weights[indx_add_table,]$charge)
  
  #use enviPat multiform function to calculate total formula of X-mer
  #form = enviPat::multiform(curformula[1], fact = nummols)
  form = enviPat::multiform(check_chemform(isotopes, form)$new_formula, fact = nummols)
  
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
  
  query_md<-mchemicaldata$mz[add_indx]-round(mchemicaldata$mz[add_indx])
  query_rt<-mchemicaldata$time[add_indx]
  
  #to account for fact that [M+1] can have a higher intensity than the monoisotopic peak,
  #first determine if [M+1] expected abundance exceeds the monoisotopic abundance and if so
  #calculate the ratio between the [M+1] and M+ peak abundances and multiply by the M+ peak intensity
  #add a further 'cushion' to ensure more-abundant features are collected during the subset procedure
  
  max_isotope_abundance = max(top_isotopes$abundance)
  if(max_isotope_abundance < 100){
    query_int<-1*(mchemicaldata$AvgIntensity[add_indx])
  }else{
    #multiply the compound intensity by the ratio of its intensity to that of its largest isotope
    ratio = max(top_isotopes$abundance) / mchemicaldata$AvgIntensity[add_indx]
    query_int = ratio * mchemdata$AvgIntensity[add_indx]
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
    
    #read in the earlier-saved Obj_for_iso_check.csv object, which was a cleaned-up (non-filled) version of dataA
    iso.check = as.data.frame(read.csv(file.path( '../Obj_for_iso_check.csv'), header = T, as.is = T, check.names = F))
    iso.check = iso.check[,-c(1)]
    
    #from the iso.check obj, get the row corresponding to the currently-considered adduct form
    iso.check$mzid = paste(round(iso.check$mz,5), iso.check$time, sep = '_')
    ref_peak_row = iso.check[which(iso.check$mzid == paste(mchemicaldata$mz[add_indx], 
                                                      mchemicaldata$time[add_indx], sep = '_')),]

    #for each isotope candidate, recover the intensity information from ....
    iso.candidates.mzid = as.list(paste(round(put_isp_masses_curmz_data$mz,5), put_isp_masses_curmz_data$time, sep='_'))
    
    iso.candidates.info = ldply(iso.candidates.mzid, function(id, iso.check, ref_peak_row){
      
      iso.candidate = iso.check[iso.check$mzid == id,]
      temp.df = rbind(ref_peak_row, iso.candidate)
      
      int.ratios = apply(temp.df[,-c(1:2, ncol(temp.df))], 2, function(col.check) {
        
        if(length(!is.na(col.check)) >= 2){
          
          int.ratio = col.check[2] / col.check[1]
          return(int.ratio)
          
        }

      })
      
      int.ratios.avrg = median(int.ratios, na.rm = T)
      
      return(data.frame('mzid' = id, 'int_ratios_avrg' = int.ratios.avrg))

    }, iso.check = iso.check, ref_peak_row)
    
    
    int_vec = iso.candidates.info$int_ratios_avrg    
 
    
    #calculate the relative intensities of peaks versus the monoisotopic mass
    #int_vec<-put_isp_masses_curmz_data$AvgIntensity
    #int_vec<-int_vec/mchemicaldata$AvgIntensity[add_indx] 
    
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
        diff_rt<-abs(put_isp_masses_curmz_data$time[isp_v]-mchemicaldata$time[add_indx])
        
        #### MODIFICATION
        
        #if the candidate isotope is indeed an isotope, this defined the +/- value associated with t
        #e.g. if isnum is 1, this is the M+1 isotope
        isnum = round((put_isp_masses_curmz_data$mz[isp_v]-mchemicaldata$mz[add_indx])*z)
        
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
            isnum2 = round((put_isp_masses_curmz_data$mz[isp_v]-mchemicaldata$mz[add_indx])*z)
            
            if(isnum2<=0){
              isp_sign<-"-"
            }else{
              isp_sign<-"+"
            }
            
            isnum2<-abs(isnum2)
            if(isnum2<=max_isp){
              
              #form_name<-as.character(paste(mchemicaldata[m,7],"_[+",(isnum),"]",sep=""))
              
              #define the adduct and isotope names (relative to the currently-considered adduct form and charge)
              colnames(top_isotopes) = gsub('m/z', 'mz', colnames(top_isotopes))
              
              top_isotopes = ddply(top_isotopes, ~mz, function(row.in, mchemicaldata.sub, z, isp_sign){
                
                iso.mz = row.in['mz'] 
                
                adduct.mz = mchemicaldata.sub$mz
                
                id = round((iso.mz-adduct.mz)*z)
                
                add_name = as.character(paste(mchemicaldata.sub$Adduct,"_[",isp_sign,(id),"]",sep=""))
                
                form_name = as.character(paste(mchemicaldata.sub$Formula,"_[",isp_sign, (id),"]",sep=""))
                
                out.df = data.frame(row.in, 'iso_name' = form_name, 'adduct_name' = add_name, check.names = F)
                
              }, mchemicaldata.sub = mchemicaldata[add_indx,], z, isp_sign)
              
               
              #iso_ppm_tol = 5
              
              
              #### FIRST CHECK PPM FROM POSSIBLE ISOTOPE ASSIGNMENTS
              iso_hits = ((put_isp_masses_curmz_data$mz[isp_v] - top_isotopes$mz) / top_isotopes$mz)* 1000000
              iso_hits = top_isotopes[which(abs(iso_hits) < iso_ppm_tol),]
              
              #IF > 1 ISOTOPE WITHIN PPM TOLERANCE, CHECK IF ANY WITHIN ISO_INT_TOL
              #start by calculating ratio of candidate isotope peak to monoisotopic peak
              
              #percent = (iso.candidates.info$int_ratios_avrg[isp_v] / mchemicaldata$AvgIntensity[add_indx])*100
              
              percent = (iso.candidates.info$int_ratios_avrg[isp_v])*100
              
              
              ##### IMPORTANT PARAMETER ###########
              #iso_int_tol = 0.05 #5% tolerance allowed
              
              #collect isotopes that occur within +/- the intensity tolerance
              iso_hits = subset(iso_hits, 
                                iso_hits$abundance > (percent - (percent * iso_int_tol)) &
                                  iso_hits$abundance < (percent + (percent * iso_int_tol)))
              
              #temp = rbind( data.frame('mz' =top_isotopes$mz, 'int'=top_isotopes$abundance, 'type' = 'ref'), 
              #              data.frame('mz' = put_isp_masses_curmz_data$mz, 'int' = (put_isp_masses_curmz_data$AvgIntensity / mchemicaldata$AvgIntensity[add_indx])*100, 
              #                         'type' = 'exp')   )
              #rownames(temp) = paste(temp$mz, temp$type, sep = '_')
              #temp
              #temp = hclust(dist(temp[,c(1:2)], diag = T))
              #cutree(temp, 5)
              
              #FORM NAME IS FOR NEUTRAL MOLECULE AS DEFINED IN THE USER-SELECTED DATABASE!
              if(nrow(iso_hits) > 0){
                
                #provides details of the elements contributing to the isotope
                isotope_elements = paste(iso_hits$name, collapse='//', sep='')
                isotope_elements = cbind(as.list(strsplit(isotope_elements, '//')[[1]]))
                
                #provides the CcHhNnOoPpSs_[+X] isotope form (neutral)
                #form_name = as.character(paste(mchemicaldata$Formula[add_indx],"_[",isp_sign,(isnum2),"]",sep=""))
                
                #COMMENTED THIS LINE!!!!
                #form_name<-as.character(paste(mchemicaldata[m,7],"_[",isp_sign,(isnum2),"]",sep=""))
                
                #end of modification
                ###########################################################
                
                temp_var = mchemicaldata[add_indx,]
                temp_var$mz = put_isp_masses_curmz_data$mz[isp_v]
                temp_var$time = put_isp_masses_curmz_data$time[isp_v]
                temp_var$time.y = put_isp_masses_curmz_data$time[isp_v]
                temp_var$AvgIntensity = put_isp_masses_curmz_data$AvgIntensity[isp_v]
                temp_var$MD = put_isp_masses_curmz_data$MD[isp_v]
                temp_var$Formula = '-'
                temp_var$isotope_elements = NULL
                
                temp_var = data.frame(temp_var, 'isotope_elements' = as.character(isotope_elements), check.names = F)
                temp_var$theoretical.mz = iso_hits$mz
                temp_var$Adduct = as.character(iso_hits$adduct_name)
                temp_var$Formula = as.character(iso_hits$iso_name)

                #this earlier approach relied on a sub-optimal approach to mapping columns
                #other_inf<-cbind(rep("-",7))
                #temp_var<-cbind(put_isp_masses_curmz_data[isp_v,c(1:2)],
                #                t(other_inf),
                #                put_isp_masses_curmz_data[isp_v,c(3:4)],
                #                put_isp_masses_curmz_data[isp_v,c(5:6)],
                #                isotope_elements) #ADDED - provides details of which elements contributed to the putative isotope annotation

                #colnames(temp_var)<-colnames(mchemicaldata)
                #temp_var$Formula<-form_name
                #temp_var$Name<-as.character(mchemicaldata$Name[add_indx])
                #temp_var$chemical_ID<-as.character(mchemicaldata$chemical_ID[add_indx])
                
                #temp_var$Adduct<-paste(mchemicaldata[m,9],"_[+",isnum,"]",sep="")
                
                #get the adduct associated with the assigned isotopes
                #adductname=mchemicaldata$Adduct[add_indx] #queryadductlist[adnum]
                
                
                #temp_var$Adduct<-paste("M","_[",isp_sign,(abs(isnum2)),"]",sep="")
                
                ###### MODIFIED THIS - incorrectly reference 'isnum' versus 'isnum2'
                #temp_var$Adduct<-paste(mchemicaldata$Adduct[add_indx],"_[",isp_sign,(abs(isnum2)),"]",sep="")
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
                  
                  
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(final_isp_annot_res)
  
}