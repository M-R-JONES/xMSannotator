#' get_mz_by_monoisotopicmass
#' 
#' Generate list of expected m/z for a specific monoisotopic mass
#' 
#' 
#' @param monoisotopicmass Monoisotopic mass. e.g.: 149.051
#' @param dbid Database or user-defined ID. e.g.: "M001"
#' @param name Metabolite name. e.g.: "Methionine"
#' @param formula Chemical formula. e.g.: "C5H11NO2S"
#' @param queryadductlist List of adducts to be used for searching.  eg:
#' c("M+H","M+Na","M+K"), c("positive") for positive adducts, c("negative") for
#' negative adducts c("all") for all adducts
#' @param syssleep Wait time between queries to prevent overloading the KEGG
#' REST interface. e.g.: 0.1
#' @return Returns an R object with a list of expected m/z for the input
#' monoisotopic mass.
#' @author Karan Uppal & Martin Jones
get_mz_by_monoisotopicmass <- function(monoisotopicmass, 
    dbid = NA, name = NA, formula = NA, queryadductlist = NA, 
    syssleep = 0.01, adduct_table = NA) {
    cnames <- c("mz", "ID", "Name", "Formula", "MonoisotopicMass", 
        "Adduct", "AdductMass")
    
    
    if (is.na(adduct_table) == TRUE) {
      load('adducts_enviPat.rda')
    }
    
    if(is.na(queryadductlist) == TRUE) {
      queryadductlist = c("M+H")
      queryadductlist = as.list(queryadductlist)
    }

    adductmass <- adduct_table$adductMass
    mult_charge <- adduct_table$charge
    num_mol <- adduct_table$num_molecules
    
    names(adductlist) <- as.character(adduct_names)
    names(mult_charge) <- as.character(adduct_names)
    names(num_mol) <- as.character(adduct_names)
    
    alladducts <- adduct_names
    
    
    if (queryadductlist[1] == "positive") {
        queryadductlist <- adduct_names[which(adduct_table$Mode == "positive")]
    } else if (queryadductlist[1] == "negative") {
      queryadductlist <- adduct_names[which(adduct_table$Mode == "negative")]
    } else if (queryadductlist[1] == "all") {
      queryadductlist <- alladducts
    } else if (length(which(queryadductlist %in% alladducts == FALSE)) > 0) {
        errormsg <- paste("Adduct should be one of:", sep = "")
        for (i in alladducts) {
          errormsg <- paste(errormsg, i, sep = " ; ")
        }
        stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)", 
          "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)", 
          "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)", 
          "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)")
    }
    
    if (is.na(dbid) == TRUE) {
        dbid = "-"
    }
    
    if (is.na(name) == TRUE) {
        name = "-"
    }
    if (is.na(formula) == TRUE) {
        formula = "-"
    }

    map_res = ldply(as.list(queryadductlist), function(adname, adduct_table, monoisotopicmass) {
      
      if(!is.na(adname)){
        
        #extract adduct information from the adduct table
        adduct.info = adduct_table[which(adduct_table$Adduct == adname),]
        
        #extract mass, charge and no. of molecules associated with adduct form
        adductmass = adduct.info$adductMass
        adductcharge = adduct.info$charge
        adductnmol = adduct.info$num_molecules
        
        #calculate mass of the adduct form
        mz = ((monoisotopicmass * as.numeric(adductnmol)) + (as.numeric(adductmass)))/as.numeric(adductcharge)
        
        # delta_ppm=(max.mz.diff)*(mz/1000000)
        # min_mz=round((mz-delta_ppm),5)
        # max_mz=round((mz+delta_ppm),5)
        
        #return data frame
        data.frame('mz' = mz, 'ID' = dbid, 'Name' = name, 'Formula' = formula, 'Adduct' = adname, 
                  'AdductMass' = adductmass, 'MonoisotopicMass' = monoisotopicmass, stringsAsFactors = F)
        
      }
    }, adduct_table = adduct_table, monoisotopicmass = monoisotopicmass)
    
    return(map_res)
}
