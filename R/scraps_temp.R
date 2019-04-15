#Scraps from end of xMSannotator_get_chemscorev1_6_73_MJ.R

k_power<-3
#length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights[,1]))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))
#chemical_score<-length(unique(cur_adducts))*length(which(cur_adducts%in%adduct_weights[,1]))*(1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power))

#chemical_score<-length(which(cur_adducts%in%adduct_weights[,1]))+1*(topquant_cor)*(1/((diff_rt*0.1)+1)^k_power)
#chemical_score<-(topquant_cor)*(1/(diff_rt+1)^k_power)*length(unique(mchemicaldata$Adduct))



#no good adducts found

chemical_score<-0
mchemicaldata<-mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust==top_mod[i]),]
mchemicaldata<-mchemicaldata[order(mchemicaldata$mz),]

cur_adducts_with_isotopes<-mchemicaldata$Adduct
cur_adducts<-gsub(cur_adducts_with_isotopes,pattern="(_\\[(\\+|\\-)[0-9]*\\])",replacement="")

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