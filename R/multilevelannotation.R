#' multilevelannotation
#' 
#' The function uses a multi-level scoring algorithm to annotate features using
#' HMDB, KEGG, T3DB, and LipidMaps.
#' 
#' Multistage clustering algorithm based on intensity profiles, retention time
#' characteristics, mass defect, isotope/adduct patterns and correlation with
#' signals for metabolic precursors and products. The algorithm uses
#' high-resolution mass spectrometry data for a series of samples with common
#' properties and publicly available chemical, metabolic and environmental
#' databases to assign confidence levels to annotation results.
#' 
#' @param dataA Peak intensity matrix. The first two columns should be "mz" and
#' "time" followed by sample intensities.
#' @param max.mz.diff Mass tolerance in ppm for database matching. e.g.: 10
#' @param max.rt.diff Retention time (s) tolerance for finding peaks derived
#' from the same parent metabolite. e.g.: 10
#' @param cormethod Method for correlation. e.g.: "pearson" or "spearman". The
#' "pearson" implementation is computationally faster.
#' @param clustermethod Method for clustering of features across samples. Valid
#' options are 'WGCNA' (default), 'hc' (hierarchical clustering) and 'graph'.
#' @param num_nodes Number of computing cores to be used for parallel
#' processing. e.g.:2
#' @param queryadductlist Adduct list to be used for database matching.  e.g.
#' c("all") for all possible positive or negative adducts or
#' c("M+H","M+Na","M+ACN+H") for specifying subset of adducts.  Run
#' data(adduct_table) for list of all adducts.
#' @param mode Ionization mode. e.g.:"pos" or "neg"
#' @param num_sets How many subsets should the total number of metabolites in a
#' database should be divided into to faciliate parallel processing or prevent
#' memory overload? e.g.: 1000
#' @param outloc Output folder location.  e.g.: "C:\Documents\ProjectX\"
#' @param db_name Database to be used for annotation. e.g.: "HMDB", "KEGG",
#' "T3DB", "LipidMaps"
#' @param adduct_weights Adduct weight matrix. Run data(adduct_weights) to see
#' an example adduct weight matrix.
#' @param allsteps If FALSE, only step 1 that involves module and retention
#' time based clustering is performed. e.g.: TRUE
#' @param corthresh Minimum correlation threshold between peaks to qualify as
#' adducts/isotopes of the same metabolite.
#' @param NOPS_check Should elemental ratio checks be performed as outlined in
#' Fiehn 2007? e.g. TRUE
#' @param customIDs Custom list of select database IDs (HMDB, KEGG, etc.) to be
#' used for annotation. This should be a data frame.  Run data(customIDs) to
#' see an example.
#' @param missing.value How are missing values represented in the input peak
#' intensity matrix? e.g.: NA
#' @param deepsplit How finely should the clusters be split? e.g.: 2 Please see
#' WGCNA for reference.
#' @param networktype Please see WGCNA for reference: e.g: "unsigned" or
#' "signed".
#' @param minclustsize Minimum cluster size. e.g: 10
#' @param module.merge.dissimilarity Maximum dissimilarity measure (i.e.,
#' 1-correlation) to be used for merging modules in WGCNA.  e.g.:0.2
#' @param filter.by Require the presence of certain adducts for a high
#' confidence match. e.g.: c("M+H")
#' @param redundancy_check Should stage 5 that involves redundancy based
#' filtering be performed? e.g.: TRUE or FALSE
#' @param min_ions_perchem Minimum number of adducts/isotopes to be present for
#' a match to be considered high confidence. e.g.:2
#' @param biofluid.location Used only for HMDB. e.g.: NA, "Blood" ,"Urine",
#' "Saliva" Set to NA to turn off this option. Please see
#' http://www.hmdb.ca/metabolites or run data(hmdbAllinf); head(hmdbAllinf) for
#' more details.
#' @param origin Used only for HMDB. e.g.: NA, "Endogenous", "Exogenous", etc.
#' Set to NA to turn off this option. Please see http://www.hmdb.ca/metabolites
#' or run data(hmdbAllinf); head(hmdbAllinf) for more details.
#' @param status Used for HMDB. e.g.: NA, "Detected", "Detected and
#' Quantified", "Detected and Not Quantified", "Expected and Not Quantified".
#' Set to NA to turn off this option. Please see http://www.hmdb.ca/metabolites
#' or run data(hmdbAllinf) for more details.
#' @param boostIDs Databased IDs of previously validated metabolites. e.g.:
#' c("HMDB00696"). Set to NA to turn off this option.
#' @param max_isp Maximum number of expected isotopes. e.g.: 5
#' @param MplusH.abundance.ratio.check Should MplusH be the most abundant
#' adduct? e.g. TRUE or FALSE
#' @param customDB Custom database. Run: data(custom_db); head(custom_db) to
#' see more details on formatting. Set to NA to turn off this option
#' @param HMDBselect How to select metabolites based on HMDB origin, biolfuid,
#' and status filters? e.g.: "all" to take union, "intersect" to take
#' intersection
#' @param mass_defect_window Mass defect window in daltons for finding
#' isotopes. e.g.: 0.01
#' @param mass_defect_mode "pos" for finding positive isotopes; "neg" for
#' finding unexpected losses/fragments; "both" for finding isotopes and
#' unexpected losses/fragments
#' @param pathwaycheckmode How to perform pathway based evaluation?  "pm":
#' boosts the scores if there are other metabolites from the same pathway that
#' are also assigned to the same module. "p": boosts the scores if there are
#' other metabolites from the same pathway without accounting for module
#' membership.
#' @return The function generates output at each stage: Stage 1 includes
#' modules and retention time based clustering of features without any
#' annotation Stage 2 includes modules and retention time based clustering of
#' features along with simple m/z based database matching Stage 3 includes
#' scores for annotations assigned in stage 2 Stages 4 and 5 include the
#' confidence levels before and after redundancy (multiple matches) filtering,
#' respectively
#' @author Karan Uppal <kuppal2@@emory.edu>
#' @references Langfelder P, Horvath S. WGCNA: an R package for weighted
#' correlation network analysis. BMC Bioinformatics. 2008 Dec 29;9:559. Wishart
#' DS et al. HMDB 3.0--The Human Metabolome Database in 2013.  Nucleic Acids
#' Res. 2013 Jan;41(Database issue):D801-7.  Kanehisa M. The KEGG database.
#' Novartis Found Symp. 2002; 247:91-101;discussion 101-3, 119-28, 244-52.
#' Review. Lim E, et al.  T3DB: a comprehensively annotated database of common
#' toxins and their targets. Nucleic Acids Res. 2010 Jan;38(Database
#' issue):D781-6. Sleno L. The use of mass defect in modern mass spectrometry.
#' J Mass Spectrom.  2012 Feb;47(2):226-36. Zhang H, Zhang D, Ray K, Zhu M.
#' Mass defect filter technique and its applications to drug metabolite
#' identification by high-resolution mass spectrometry. J Mass Spectrom. 2009
#' Jul;44(7):999-1016.
#' @keywords ~xMSannotator ~multilevelannotation ~annotation
multilevelannotation <- function(dataA, max.mz.diff = 10, 
    max.rt.diff = 10, cormethod = "pearson", num_nodes = 2, 
    queryadductlist = c("all"), gradienttype = "Acetonitrile", 
    mode = "pos", outloc, db_name = "HMDB", adduct_weights = NA, 
    num_sets = 3000, allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE, 
    customIDs = NA, missing.value = NA, missing.permute = 100, deepsplit = 2, networktype = "unsigned", 
    minclustsize = 10, module.merge.dissimilarity = 0.2, 
    filter.by = c("M+H"), redundancy_check = TRUE, min_ions_perchem = 1, 
    biofluid.location = NA, origin = NA, status = NA, boostIDs = NA, 
    max_isp = 5, MplusH.abundance.ratio.check = FALSE, customDB = NA, 
    HMDBselect = "union", mass_defect_window = 0.01, mass_defect_mode = "pos", 
    dbAllinf = NA, pathwaycheckmode = "pm") {
  
  options(warn = -1)
  allowWGCNAThreads(nThreads = num_nodes)
  
  #load the adducts table
  load("adducts_enviPat.rda")
  setwd('~')
  
  ###### MOCK STRUCTURE OF TABLE ####
  # | Adduct     | num_molecules | charge | adductMass |  Mode     | Type | Merge_add | Merge_sub |
  # | M          | 1             | 1      | 0.000000   |  neutral  | S    | <NA>      | <NA>      |
  # | ...        | ...           | ...    | ...        |  ...      | ...  | ...       | ...       |
  # | 2M-2H+NH4  | 2             | 1      | 16.019271  |  negative | S    | N1H4      | H2        |
  
  #define outloc location (output folder path)
  if(outloc == ''){
    setwd('~')
    outloc = file.path(getwd(), 'temp/outputs/', sep = '')
    print(paste('results are saved in: ', outloc, sep = ''))
  }
  
  #Create the output directory and set as the current working directory
  suppressWarnings(dir.create(outloc))
  setwd(outloc)
  
  #check if dataA exists and if not, use the faahko dataset as a model.
  
  if(exists('dataA') == FALSE){
    
    print('Using faahko package data to test workflow')
    library(faahKO)
    library(xcms)
    library(CAMERA)
    
    #group peaks across samples
    xs = group(faahko)
    #convert the grouped peaks matrix in to a xsAnnotate object and then export the peaklist
    dataA = getPeaklist(xsAnnotate(xs))
    
    #refine the peaklist to include the mz, retention time and intensities for each feature in the dataset
    col_indxs = c(which(colnames(dataA) == 'mz'), which(colnames(dataA) == 'rt'), 
                  seq((which(colnames(dataA) == 'WT')+1), (which(colnames(dataA) == 'isotopes')-1), 1))
    
    #subset the dataframe to keep only those columns corresponding to mz, rt and intensities across samples, for each feature
    dataA = dataA[,col_indxs]
    
    #convert column name from 'rt' to 'time'
    colnames(dataA)[which(colnames(dataA) == 'rt')] = "time"
    
    #round mz values to 7 d.p.
    dataA$mz <- round(dataA$mz, 7)
    
    #round retention time values to 1 d.p.
    dataA$time <- round(dataA$time, 1)
  }
  
  #round all intensity values to 1 decimal place
  dataA[, -c(1:2)] <- round(dataA[, -c(1:2)], 1)
    
  #define the dissimiliarity metric for merging dendrogram branches during clustering
  if (is.na(module.merge.dissimilarity) == TRUE) {
    module.merge.dissimilarity = 1 - corthresh
  }
  
  
  if (is.na(customIDs) == FALSE){
    customIDs = as.data.frame(customIDs)
  }
  
  print(paste("Annotating using ", db_name, " database:", sep = ""))
    
  max_diff_rt <- max.rt.diff
  cutheight = 1 - corthresh  #module.merge.dissimilarity
  time_step = 1
  step1log2scale = FALSE
    
  ###### Seemingly not used #######
  # if (FALSE) {
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "2M\\+2H", replacement = "M+H")
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "3M\\+3H", replacement = "M+H")
  #     
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "2M\\+2Na", replacement = "M+Na")
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "3M\\+3Na", replacement = "M+Na")
  #     
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "2M\\+2K", replacement = "M+K")
  #     adduct_table$Adduct <- gsub(adduct_table$Adduct, 
  #         pattern = "3M\\+3K", replacement = "M+K")
  # }
  #
  # adduct_table <- unique(adduct_table)

  #create a subset of valid adducts to consider
  if (queryadductlist == "all" & mode == "pos") {
      adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "positive") | 
                                            (adduct_table$Type == gradienttype & adduct_table$Mode == "positive")]
      adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
  } else {
      if (queryadductlist == "all" & mode == "neg") {
          adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "negative") | 
                                                (adduct_table$Type == gradienttype & adduct_table$Mode == "negative")]
          adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
      } else {
          adduct_names <- adduct_table$Adduct[which(adduct_table$Adduct %in% queryadductlist)]
          adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
      }
  }
    
  #ensure that adduct names are unique after previous filtering step
  adduct_names <- unique(adduct_names)
    
  #set up objects for processing
  res_list <- list()
  db_name_list = db_name
  
  #No value in creating a second object with the same string
  #outloc_allres <- outloc
    
  #check whether results from step1 of the processing pipeline have already been generated
  l1 <- list.files(outloc)
  check_step1 <- which(l1 == "step1_results.Rda")
        
  if(length(check_step1) < 1){
    
    print('Checking if level1A results exist')
    
    check_levelA = which(l1 == 'xMSannotator_levelA_modules.Rda')
    
    if (length(check_levelA) < 1) {
      
      print('Level1A result did not exist, therefore computing modules using WGCNA')
      
      ### PRE-PROCESS THE DATA MATRIX TO ENABLE CLUSTERING OF FEATURES
      
      #replace missing values (valid options for missing value replacement: 'NA' and 'mean')
      #strongly encourage users to use alternative approaches based on robust multivariate imputation
      
      if (is.na(missing.value) == TRUE) {
        
        ## when 'missing.value' is set to NA, replace all missing values with NA
        # causes clustering to fail when too many missing values occur in the matrix
        dataA <- replace(as.matrix(dataA), which(dataA == missing.value), NA)
        dataA <- as.data.frame(dataA)
        
      } else if(missing.value == 'mean'){
        
        #for each feature, set missing value to mean of all non-missing values
        dataA = adply(dataA, 1, function(row){
          
          intens = row[,-c(1:2)]
          
          #only consider rows where > minimum number of values have been recorded
          if(length(which(!is.na(intens)==T)) > 0){
            
            #calculate row mean, excluding missing values
            m = rowMeans(intens, na.rm = TRUE)
            #replace all missing values (indicated by NA) with the mean value 'm'
            intens[is.na(intens)] = m
            #return intens object
            return(intens)    
          }
        })
      }
      
      
      #create a series of unique mzid strings that map directly to the dataA object
      mzid <- paste(dataA$mz, dataA$time, sep = "_")
      
      #remove mzid values if identical mz and retention time pairs exist (extent of rounding might cause inaccuracy here)
      if (length(which(duplicated(mzid) == TRUE)) > 0) {
        dataA <- dataA[-which(duplicated(mzid) == TRUE), ]
      }
      
      #generate string comprising mz_retention time        
      mzid <- paste(dataA$mz, dataA$time, sep = "_")
      
      #perform fast calculation of the correlation coefficients between all features
      system.time(global_cor <- WGCNA::cor(t(dataA[,-c(1:2)]), nThreads = num_nodes, method = cormethod, use = "p"))
      
      #round correlation coefficients to 2 decimal places
      global_cor <- round(global_cor, 2)
      
      #insert feature details (mz_rt) for each row and column of the correlation matrix.
      colnames(global_cor) <- mzid
      rownames(global_cor) <- mzid
      
      #set the working directory to the output directory and save the correlation matrix
      setwd(outloc)
      
      save(global_cor, file = "global_cor.Rda")
      
      mycl_metabs = NA
      
      #THIS FAILS IF MISSING VALUES ARE FOUND, HENCE THE OPTION FOR FILLING WITH 'MEAN' WAS ADDED ABOVE
      
      ######## IMPORTANT: changed distance metrix from (1-global_cor) to (1 - abs(global_cor))
      #the original approach used 1 - global_cor as the distance. Surely for negative correlations we get the wrong distance
      #e.g. correlation coefficient of -0.3 becomes: 1 - (-0.3) = 1.3, implying an impossibly good distance.
      #hr = flashClust(as.dist(1 - global_cor), method = "complete")
      hr = flashClust(as.dist(1 - abs(global_cor)), method = "complete") 
      
      #as above, changed the distance matrix to use 1 minus the absolute correlation coefficient
      dissTOMCormat <- (1 - abs(global_cor)) ######## IMPORTANT: changed distance metrix from (1-global_cor) to (1 - abs(global_cor))
      
      a1 <- apply(global_cor, 1, function(x) {
        length(which(x > corthresh))
      })
      
      fname_c <- paste("NumberOfCorrelationsPerFeature_cor", corthresh, ".csv", sep = "")
      
      write.csv(a1, file = fname_c)
      rm(global_cor)
      
      # dynamically cut the dendrogram generated during clustering, using the distM matrix to define closeness of features
      mycl_metabs <- cutreeDynamic(hr, distM = dissTOMCormat, 
                                   deepSplit = 1, minClusterSize = minclustsize, 
                                   pamRespectsDendro = FALSE, pamStage = TRUE, 
                                   verbose = 0)
      
      if (clustmethod == "WGCNA") {
            
        #####REMOVE DURING CLEAN UP!!!!
        #functions.dir = '\\\\eawag\\userdata\\jonesmar\\My Documents\\1_cyano_peps\\6_Informatics\\Git_xMSannotator\\xMSannotator\\R\\'
        #temp = list.files(functions.dir, pattern = '.R')
        #temp =temp[-which(temp=='multilevelannotation.R')]
        #temp =temp[-which(temp=='multilevelannotation_MJ.R')]
        
        #lapply(temp, function(x, path){
        #  source(paste(path, x, sep =''))
        #}, path = functions.dir )
            
        levelA_res <- get_peak_blocks_modulesvhclust(dataA = dataA, 
          simmat = NA, adjacencyfromsimilarity = FALSE, 
          time_step = time_step, max.rt.diff = max_diff_rt, 
          outloc, column.rm.index = NA, cor.thresh = NA, 
          deepsplit = deepsplit, minclustsize = minclustsize, 
          cutheight = cutheight, cormethod = cormethod, 
          networktype = networktype, num_nodes = num_nodes, 
          step1log2scale = step1log2scale, mycl_metabs = mycl_metabs)
        
        setwd(outloc)
        levelA_res <- levelA_res[, -c(1:4)]
        save(levelA_res, file = "xMSannotator_levelA_modules.Rda")
  
      } else if (clustmethod == "graph") {
        
        levelA_res <- get_peak_blocks_graph(dataA, 
                                            simmat = global_cor, adjacencyfromsimilarity = TRUE, 
                                            time_step = 3, max.rt.diff = max_diff_rt, 
                                            outloc, column.rm.index = NA, cor.thresh = NA, 
                                            cormethod = cormethod, networktype = networktype, 
                                            num_nodes = num_nodes)
                    
      } else if (clustmethod == "hc") {
                
        levelA_res <- get_peak_blocks_hclust(dataA, 
                                             time_step = time_step, max.rt.diff = max_diff_rt, 
                                             outloc, column.rm.index = NA, cor.thresh = NA, 
                                             deepsplit = deepsplit, minclustsize = minclustsize, 
                                             cutheight = cutheight, cormethod = cormethod)
      } else {
        
        print('Clustmethod must be set to one of: "HMDB", "KEGG", "Custom"')
        
      }
      
      #clean up some bulky items
      rm(dissTOMCormat)
      rm(hr)
      rm(a1)
      
      #here the row name labels get screwed up by re-ordering during processing
      #hence, re-assign correct labels!
      rownames(levelA_res) = paste(levelA_res$mz, levelA_res$time, sep = '_')
      write.csv(levelA_res, file = "Stage1.csv", row.names = FALSE)
      
      print(paste('Saving level1A results to: ', file.path(outloc, 'xMSannotator_levelA_modules.Rda', fsep = ''), sep =''))
      print(paste('Saving level1A results to: ', file.path(outloc, 'Stage1.csv', fsep = ''), sep =''))

      
    } else {
    
      print('Level1A results existed - loading from xMSannotator_levelA_modules.Rda')
      #import existing results file from outloc
      setwd(outloc)
      load("xMSannotator_levelA_modules.Rda")
  
    }
  
    print(paste('Mapping level1A results to database', db_name, sep=' '))
    
    setwd(outloc)
    
    #reorder levelA_res data frame based on mass-to-charge and then on retention time of features
    levelA_res <- levelA_res[order(levelA_res$mz,levelA_res$time), ]
    levelA_res <- levelA_res[, c(1:3)]
    
    #calculate average intensity of each feature in the dataA data frame
    dataA <- dataA[order(dataA$mz, dataA$time), ]
    
    #original storage of mean_int_vec (only way of mapping to mz and retention time was to order by mz then by tr)
    mean_int_vec <- apply(dataA[, -c(1:2)], 1, function(x) {mean(x, na.rm = TRUE)})
    
    #new mean_int_vec object that maintained relationship between intensity, mz and time
    mean_int_vec = data.frame('mean_int_vec' = mean_int_vec, dataA$mz, dataA$time)
    
    #keep only the mz and retention of each feature
    dataA <- dataA[, c(1:2)]
  
    #begin suspect screening against the database defined by db_name    
    if (db_name == "HMDB") {
        
      data(hmdbAllinf)
      
      hmdbAllinf$Name <- gsub(hmdbAllinf$Name, pattern = "[\\\"']", replacement = "")
      hmdbAllinfv3.6 = hmdbAllinf
      
      #clean up
      rm(hmdbAllinf)
      rm(hmdbAllinf, envir = .GlobalEnv)
      
      suppressWarnings(if (is.na(customIDs) == TRUE) {
        
        customIDs <- hmdbAllinfv3.6[, c(1, 20)]
        
        #extract information for location or set to ''
        if (is.na(biofluid.location) == FALSE){
          gres_location <- gregexpr(hmdbAllinfv3.6$BioFluidLocation, pattern = biofluid.location, ignore.case = TRUE)
          gres_location <- which(gres_location == 1)
        } else{
          gres_location <- ''
        }
        
        #extract information for origin or set to ''
        if (is.na(origin) == FALSE){
          gres_origin <- gregexpr(hmdbAllinfv3.6$Origin, pattern = origin, ignore.case = TRUE)
          gres_origin <- which(gres_origin == 1)
        }else{
          gres_origin <- ''
        }
        
        #extract HMDB status
        if (is.na(status) == FALSE) {
          gres_status <- gregexpr(hmdbAllinfv3.6$HMDBStatus, pattern = status, ignore.case = TRUE)
          gres_status <- which(gres_status == 1)
        }else{
          gres_status = ''
        }
        
        #generate a list for HMDB location, origin and status
        gres = list('gres_loc' = gres_location, #e.g. gres_location = c(1,20,10,23)
                    'gres_orig' = gres_origin, #e.g. gres_orig = c(20,40,60,80,100)
                    'gres_status' = gres_status) #e.g. gres_status = ''
        
        #empty list for capturing only the HMDB status, location and origin lists that != ''
        gres_merge = list()
        #remove lists == ''
        for(i in 1:length(gres)){
          if(gres[[i]] != ""){
            nm = names(gres)[i]
            gres_merge[[nm]] = gres[[i]]
          }
        }
        
        #if intersect is true, get the common indices for HMDB compounds from the origin, location and status list
        if(HMDBselect == 'intersect'){
          gres = Reduce(intersect, gres_merge)
        }else{
          gres = unique(unlist(gres_merge)) #otherwise, get all unique HMDB compounds
        }
        
        #keep found compounds
        if(!is.null(gres)){
          customIDs <- hmdbAllinfv3.6[gres, c(1, 20)]
        }else{
          customIDs = ''
        }
      
      })
                
      rm(hmdbAllinfv3.6, envir = .GlobalEnv)
      data(hmdbCompMZ)
      
      hmdbCompMZ$mz <- round(as.numeric(as.character(hmdbCompMZ$mz)), 5)
      
      hmdbCompMZ$Name <- gsub(hmdbCompMZ$Name, pattern = "[\\\"']", replacement = "")
      
      
      suppressWarnings(if (is.na(customIDs) == FALSE) {
  
        customIDs <- unique(customIDs)
        hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$HMDBID %in% customIDs[, 1]), ]
        hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names), ]
        
      })
    
      hmdbCompMZ <- hmdbCompMZ[which(hmdbCompMZ$Adduct %in% adduct_names), ]
                  
      chemCompMZ <- hmdbCompMZ
      
      print("Dimension of precomputed HMDB m/z database")
      print(dim(chemCompMZ))
                  
      hmdbCompMZ <- 1
      hmdbAllinf <- 1
      
      try(rm(hmdbCompMZ), silent = TRUE)
      try(rm(hmdbCompMZ, envir = .GlobalEnv), silent = TRUE)
      try(rm(hmdbAllinf, envir = .GlobalEnv), silent = TRUE)
    
    } else if (db_name == "KEGG") {
  
      data(keggCompMZ)
      
      #round all KEGG compound mz values to 5 d.p.            
      keggCompMZ$mz <- round(as.numeric(as.character(keggCompMZ$mz)), 5)
      
      #fix the naming of KEGG compounds
      keggCompMZ$Name <- gsub(keggCompMZ$Name, pattern = "[\\\"']", replacement = "")
      
      suppressWarnings(if (is.na(customIDs) == FALSE) {
        
        #map KEGG compound mz values to those provided in the user-specified customIDs object (column#1 = mz values in custom list)
        customIDs <- unique(customIDs)
        keggCompMZ <- keggCompMZ[which(keggCompMZ$KEGGID %in% customIDs[, 1]), ]
        keggCompMZ <- keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names), ]
      
      })
      
      #confirm which KEGG compound adducts are found in the adduct names list
      keggCompMZ <- keggCompMZ[which(keggCompMZ$Adduct %in% adduct_names), ]
      chemCompMZ <- keggCompMZ
      
      print("Dimension of precomputed KEGG m/z database")
      print(dim(chemCompMZ))
                  
      rm(keggCompMZ)
      rm(keggCompMZ, envir = .GlobalEnv)
        
    } else if (db_name == "LipidMaps") {
      
      #import the lipidmaps compounds list
      data(lipidmapsCompMZ)
        
      #map lipidmaps adduct forms to those in the adduct_names list (i.e. confirm which lipidmaps adducts are being considered)
      lipidmapsCompMZ <- lipidmapsCompMZ[which(lipidmapsCompMZ$Adduct %in% adduct_names), ]
        
      #ensure lipidmap compound mz values are rounded to 5 d.p.
      lipidmapsCompMZ$mz <- round(as.numeric(as.character(lipidmapsCompMZ$mz)), 5)
      
      #fix naming of lipidmap compounds
      lipidmapsCompMZ$Name <- gsub(lipidmapsCompMZ$Name, pattern = "[\\\"']", replacement = "")
    
      chemCompMZ <- lipidmapsCompMZ
    
      print("Dimension of precomputed LipidMaps m/z database")
      print(dim(chemCompMZ))
    
      #clean up
      rm(lipidmapsCompMZ)
      rm(lipidmapsCompMZ, envir = .GlobalEnv)
    
    } else if (db_name == "T3DB") {
                      
      #import T3db compounds
      data(t3dbCompMZ)
      
      #ensure only those adducts in adduct_names object are being considered from the T3db object
      t3dbCompMZ <- t3dbCompMZ[which(t3dbCompMZ$Adduct %in% adduct_names), ]
      
      #round t3db compound mz values to 5 d.p.      
      t3dbCompMZ$mz <- round(as.numeric(as.character(t3dbCompMZ$mz)), 5)
      
      #fix names for t3db compounds
      t3dbCompMZ$Name <- gsub(t3dbCompMZ$Name,pattern = "[\\\"']", replacement = "")
      chemCompMZ <- t3dbCompMZ
      
      print("Dimension of precomputed T3DB m/z database")
      print(dim(chemCompMZ))
      
      rm(t3dbCompMZ)
      rm(t3dbCompMZ, envir = .GlobalEnv)
      
      
    } else if (db_name == "Custom") {
      
      #read in suspect screening database
      customDB = customDB
      
      colnames(customDB)[which(colnames(customDB)=='ID')] = 'CompoundID'
      colnames(customDB)[which(colnames(customDB)=='Name')] = 'Name'
      colnames(customDB)[which(colnames(customDB)=='Formula')] = 'Formula'
      colnames(customDB)[which(colnames(customDB)=='MonoisotopicMass')] = 'MonoisotopicMass'
      
      #outline of custom database format
      # CompoundID    |   Name            |   Monoisotopic mass   | Formula     | 
      # CID1001       |   Phenethylamine  |   121.089             | C8H11N      |
      # CID6305       |   Tryptophan      |   204.09              | C11H12N2O2  |
      # ...           |   ...             |   ...                 | ...         |
      # CID6057       |   tyrosine        |   181.074             | C9H11NO3    |
      
      #remove any duplicate values in the database (based on compound ID value)
      if (length(which(duplicated(customDB$CompoundID) == TRUE)) > 0) {
        inputmassmat <- customDB[-which(duplicated(customDB$CompoundID) == TRUE), ]
      }
      
      #for each compound in the custom database, extract the mz, database ID, name and formula
      if(length(which(colnames(customDB)=='CompoundID')) == 0){
        stop("column name 'CompoundID' is missing in the custom database")
      }
      
      if(length(which(colnames(customDB)=='Formula')) == 0){
        stop("column name 'Formula' is missing in the custom database")
      }
      
      if(length(which(colnames(customDB)=='MonoisotopicMass')) == 0){
        stop("column name 'MonoisotopicMass' is missing in the custom database")
      }
      
      if(length(which(colnames(customDB)=='Name')) == 0){
        stop("column name 'Name' is missing in the custom database")
      }
      
      #create suspect screening database (update overcomes annoyance of having all columns as factors)
      mz_search_list <- ddply(customDB, ~CompoundID, function(row.in){
        
        #extract compound information from the input database and define the skeleton for the dataframe containing this information
        mz_search_list <- get_mz_by_monoisotopicmass(monoisotopicmass = row.in$MonoisotopicMass, 
                                                     dbid = row.in$CompoundID, 
                                                     name = row.in$Name, 
                                                     formula = row.in$Formula,
                                                     queryadductlist = adduct_names, 
                                                     adduct_table = adduct_table)
        
        return(mz_search_list)
      
      })
      
      #write the suspect screening dataframe to "mz_search_list.rda" object
      save(mz_search_list, file = "mz_search_list.Rda")
      print(paste('Database successfully transformed in to data frame (obj. mz_search_list). Saving to: ', 
                  file.path(outloc, 'mz_search_list.Rda', fsep =''), sep = ''))
        
      #define chemCompMZ as the input suspect screening list
      chemCompMZ <- mz_search_list
        
      #clean up
      rm(mz_search_list)
      rm(customDB)
      rm(customDB, envir = .GlobalEnv)
  
    } else {
        
        stop("db_name should be: KEGG, HMDB, T3DB, LipidMaps, or Custom")
    }

    #if adduct_weights was not defined when executing function, set [M+/-H] as the only adducts with defined weights
    if (is.na(adduct_weights) == TRUE) {
        data(adduct_weights)
        adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
        adduct_weights1[1, ] <- c("M+H", 1)
        adduct_weights1[2, ] <- c("M-H", 1)
        adduct_weights <- as.data.frame(adduct_weights1)
        colnames(adduct_weights) <- c("Adduct", "Weight")
        adduct_weights <- as.data.frame(adduct_weights)
        rm(adduct_weights1)
    }
    
    #fix naming of some compounds
    chemCompMZ$Name <- gsub("([-:();])|[[:punct:]]", "\\1", chemCompMZ$Name)
    chemCompMZ$Formula <- gsub("([-:();])|[[:punct:]]","\\1", chemCompMZ$Formula)
    
    #ensure naming of columns is consistent with earlier function names
    if(grep('CompoundID', colnames(chemCompMZ)) > 0){
      chemCompMZ$chemical_ID = NULL
      colnames(chemCompMZ)[grep('CompoundID', colnames(chemCompMZ))] = 'chemical_ID'
    }else{
      #original command to rename column from ID to chemical_ID
      colnames(chemCompMZ)[which(colnames(chemCompMZ) == 'ID')] = 'chemical_ID'
    }
    
    #if less-than 1000 peaks in the matrix, set the randsetsize to that value, otherwise leave as 1000
    if (dim(dataA)[1] < 1000) {
      randsetsize <- dim(dataA)[1]
    }else{
      randsetsize <- 1000
    }
  
    setwd(outloc)
  
    #target list and associated adducts are saved
    save(chemCompMZ, file = "chemCompMZ.Rda")
    
    #confirm whether levelB checks have already been performed
    l1 <- list.files(outloc)
    
    print('Checking for existing level1B results')
    
    check_levelB <- which(l1 == "xMSannotator_levelB.Rda")
  
    if (length(check_levelB) < 1) {
      
      print(paste('Level1B results NOT found in', outloc, '- performing level1B procedures', sep = ''))
      
      #set up cluster for parallel processing purposes
      cl <- makeSOCKcluster(num_nodes)
      
      #load all required packages
      clusterApply(cl, list('XML', 'R2HTML', 'RCurl', 'SSOAP', 'limma', 'plyr', 'png'), require, character.only = T)
      
      clusterEvalQ(cl, library('plyr'))
  
      #setwd(file.path(getwd(), '/', fsep=''))
      #funcs = list.files(getwd(), full.names = F, pattern = '.R', recursive = F)
      #funcs = funcs[-which(funcs=='multilevelannotation.R')]
      #funcs = funcs[-which(funcs=='xMSannotator_multilevelannotation_MJ.R')]
      #lapply(funcs, function(x) {source(x)})
      
      #parse required functions to cluster
      clusterEvalQ(cl, "processWSDL")
      clusterExport(cl, "Annotationbychemical_IDschild_multilevel")
      clusterExport(cl, "Annotationbychemical_IDschild")
      clusterExport(cl, "find.Overlapping.mzs")
      clusterExport(cl, "find.Overlapping.mzsvparallel")
      clusterExport(cl, "overlapmzchild")
      clusterExport(cl, "getVenn")
      
      if (length(which(duplicated(chemCompMZ$Formula) == TRUE)) > 0) {
          
        #with new approach to making database, no need for this step        
        #if (db_name == "Custom") {
        #  chemCompMZ$mz <- as.numeric(as.character(chemCompMZ$mz))
        #}
        
        chemCompMZ_unique_formulas <- chemCompMZ[-which(duplicated(chemCompMZ$Formula) == TRUE), ]
        
        chemCompMZ_unique_formulas <- rbind(chemCompMZ_unique_formulas, chemCompMZ[which(chemCompMZ$chemical_ID %in% 
            chemCompMZ_unique_formulas$chemical_ID),])
        
        chemCompMZ_unique_formulas <- unique(chemCompMZ_unique_formulas)
        chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[order(chemCompMZ_unique_formulas$chemical_ID),]
      } else {
        chemCompMZ_unique_formulas <- chemCompMZ
      }
                  
                  
      # save(chemCompMZ_unique_formulas,file='chemCompMZ_unique_formulas.Rda')
      save(chemCompMZ, file = "chemCompMZ.Rda")
      rm(chemCompMZ)
      
      #set up a unique string for each formula in the suspect screening database
      formula_table <- table(chemCompMZ_unique_formulas$Formula)
      uniq_formulas <- names(formula_table)
      formula_ID <- paste("Formula", seq(1:length(uniq_formulas)), sep = "_")
      
      #generate a dataframe linking unique suspect screening formulae with unique formulae ID values
      formula_id_mat <- data.frame("Formula_ID" = formula_ID, "Formula" = uniq_formulas)
      
      #associate each chemical ID with a unique formula ID (IMPORTANT TABLE FOR TRACEBACK)
      #THIS APPROACH ENSURES THAT COMPOUNDS WITH IDENTICAL FORMULAE BUT DIFFERENT NAMES, ARE RETAINED IN THE DATA FRAME
      chemCompMZ_unique_formulas <- merge(chemCompMZ_unique_formulas, formula_id_mat, by = "Formula")
      chemCompMZ_unique_formulas$chemical_ID <- chemCompMZ_unique_formulas$Formula_ID
      chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[, c("mz", "chemical_ID", "Name", "Formula", 
          "MonoisotopicMass", "Adduct", "AdductMass")]
      
      #fix name of mz column here in order to prevent two 'mz' column names being generated during parLapply below
      #if two columns with same name are generated, then when one merges the lists to a data frame, the second mz column is dropped
      colnames(chemCompMZ_unique_formulas)[which(colnames(chemCompMZ_unique_formulas)=='mz')] = 'theoretical.mz'
      
      #generate factor object containing all unique chemical ID values
      chemIDs <- unique(chemCompMZ_unique_formulas$chemical_ID)
      
      s1 <- seq(1, length(adduct_names))
      
      ######################## STAGE 2 OF SCORING PROCESS BEGINS HERE ######################
      
      print("Status 2: Mapping m/z to metabolites:")
                  
      # annotation_multilevel_ save(list=ls(),file='test.Rda')
      # system.time(l2<-parLapply(cl,s1,Annotationbychemical_IDschild_multilevel,dataA=dataA,
      # queryadductlist=c(adduct_names),adduct_type=c('S',gradienttype),
      # adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,chemIDs=chemIDs,num_nodes=num_nodes))
      
      # a1<-Annotationbychemical_IDschild_multilevel(1,dataA=dataA,
      # queryadductlist=c(adduct_names),adduct_type=c('S',gradienttype),
      # adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,chemIDs=chemIDs,num_nodes=num_nodes)
      
      # a1<-Annotationbychemical_IDschild(1,dataA=dataA,
      # queryadductlist=c(adduct_names),adduct_type=c('S',gradienttype),
      # adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,keggCompMZ=chemCompMZ_unique_formulas,otherdbs=FALSE,otherinfo=FALSE,num_nodes=num_nodes)
      
      #run stage2 of the processing            
      l2 <- parLapply(cl, s1, Annotationbychemical_IDschild, 
        dataA = dataA, queryadductlist = adduct_names, 
        adduct_type = c("S", gradienttype), adduct_table = adduct_table, 
        max.mz.diff = max.mz.diff, outloc = outloc, 
        keggCompMZ = chemCompMZ_unique_formulas, 
        otherdbs = FALSE, otherinfo = FALSE, num_nodes = num_nodes)
      
      stopCluster(cl)
      
      levelB_res = ldply(l2, rbind)
      
      #clean up
      rm(chemCompMZ)
      rm(l2)
            
      #inform user if no matches found
      if (nrow(levelB_res) < 1) {
        stop("No matches found.")
      }
    
      #define the MatchCategory column (all are assigned multiple at this point 
      #this is updated to unique later, if just one DB match found)
      levelB_res$MatchCategory = 'Multiple'
    
      #levelB_res$mz <- as.numeric(as.character(levelB_res$mz))
      #levelB_res$time <- as.numeric(as.character(levelB_res$time))
      #levelB_res <- as.data.frame(levelB_res)
      #levelB_res <- cbind(levelB_res, MatchCategory)
      
      print("DB matches")
      levelB_res <- unique(levelB_res)
      print(dim(levelB_res))
              
      uniq_formula <- as.character(unique(levelB_res$Formula))
                  
      bad_formula <- which(is.na(uniq_formula) == TRUE)
      
      if (length(bad_formula) > 0) {
        uniq_formula <- uniq_formula[-c(bad_formula)]
      }
                  
      #set up for checking formula validity
      cl <- makeSOCKcluster(num_nodes)
      clusterExport(cl, "check_golden_rules")
      clusterExport(cl, "check_element")
                  
      levelB_res_check <- parLapply(cl, 1:length(uniq_formula), function(j, uniq_formula, NOPS_check) {
        
        curformula <- as.character(uniq_formula[j])
        return(check_golden_rules(curformula, NOPS_check = NOPS_check))
      
      }, uniq_formula = uniq_formula, NOPS_check = NOPS_check)
                
  
      stopCluster(cl)
    
      #convert list to data frame
      levelB_res_check2 <- ldply(levelB_res_check, rbind)
                  
      #keep only those hits that had valid molecular formulae based on golden rules and element check
      levelB_res_check3 <- levelB_res_check2[which(levelB_res_check2$bool_check == 1), ]
                  
      #retaim those hits with valid formulae
      levelB_res <- levelB_res[which(levelB_res$Formula %in% as.character(levelB_res_check3$curformula)),]
      
      #
      water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
      water_adduct_ind <- which(levelB_res$Adduct %in% water_adducts)
                  
      cl <- makeSOCKcluster(num_nodes)
      clusterExport(cl, "check_element")
      
      if (length(water_adduct_ind) > 0) {
        
        levelB_res2 <- levelB_res[c(water_adduct_ind),]
        levelB_res <- levelB_res[-c(water_adduct_ind),]
        sind1 <- seq(1:dim(levelB_res2)[1])
        
        levelB_res_check3 <- parLapply(cl, sind1, function(j, levelB_res2) {
            
          #confirm that a valid elemental formula has been assigned
          adduct <- as.character(levelB_res2$Adduct[j])
          curformula <- as.character(levelB_res2$Formula[j])
          numoxygens <- check_element(curformula, "O")
            
          #this check should be toggled / accessible in the function variables (some people work with non-oxygen containing species)
          if (numoxygens > 0) {
            bool_check <- 1
          } else {
            bool_check <- 0
          }
          
          #assigned colnames here to avoid possible referencing errors later in the process 
          #return data frame indicating which are valid and invalid formulae, with validity indicated by bool(ean) col.
          res <- data.frame('Formula' = curformula, 'bool' = bool_check)
          res <- as.data.frame(res)
          return(res)
        }, levelB_res2)
        
        #merge lists in to data frame
        levelB_res_check4 <- ldply(levelB_res_check3,rbind)
    
        if (length(which(levelB_res_check4$bool == 1)) > 0) {
          #keep only valid formulae
          levelB_res_check4 <- levelB_res_check4[which(levelB_res_check4$bool == 1), ]
          #subset the earlier levelB_res2 object to keep those assignments associated with a valid formula
          valid_form <- which(levelB_res2$Formula %in% as.character(levelB_res_check4$Formula))
        }else{
          valid_form <- {}
        }
    
        #added option for instances where invalid formulae were detected (might be reasonable)
        if (length(valid_form) > 0) {
          levelB_res2 <- levelB_res2[valid_form,]
          levelB_res <- rbind(levelB_res, levelB_res2)
        }else{
          print('no valid formulae according to Golden Rules, however, returning all formulae assigned')
          levelB_res #check if this object is in the same format as rbind(levelB_res, levelB_res2)
        }
      }
        
      #unnecessary re-naming
      #colnames(levelB_res) <- c("theoretical.mz","chemical_ID", "Name", 
      #                          "Formula", "MonoisotopicMass", "Adduct", 
      #                          "AdductMass", 'mz', "time", "MatchCategory")
         
      #save object to outloc     
      save(levelB_res, file = "xMSannotator_levelB.Rda")
                
    } else {
      
      #if Stage 2 processing results already existed, load these in to globalenv.
      print("Status 2: Using existing m/z mapping results:")
      load("xMSannotator_levelB.Rda")
    }
            
    #general clean up
    try(rm(hmdbAllinf, envir = .GlobalEnv), silent = TRUE)
    try(rm(hmdbAllinfv3.6), silent = TRUE)
    try(rm(hmdbCompMZ), silent = TRUE)
    
    #round mz values in levelA and levelB results
    levelA_res$mz <- round(levelA_res$mz, 5)
    levelB_res$mz <- round(levelB_res$mz, 5)
            
    #determine number of results (rows) in levelB
    last_col_ind <- dim(levelB_res)[2]
            
    levels_A <- names(table(levels(levelA_res$Module_RTclust)))
    
    levels_A <- names(levels_A)
    levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ] #not required if we re-define mzid here (less risky)
    
    #re-define the mzid values here - this seems safer than relying on correct ordering from earlier generated mzid
    mzid = paste(levelA_res$mz, levelA_res$time, sep = '_')
    
    #called from earlier-generated mzid (generated during clustering stage)        
    if (length(which(duplicated(mzid) == TRUE)) > 0) {
      #remove duplicated mz-retention time pairs
      levelA_res <- levelA_res[-which(duplicated(mzid) == TRUE), ]
      dataA <- dataA[-which(duplicated(mzid) == TRUE), ]
    }
    
    #order levelA_res by Module_RTclust
    levelA_res <- levelA_res[order(levelA_res$Module_RTclust),]
    levels_A <- levels(unique(levelA_res$Module_RT))
    
    #extract the module number from the Module_RTcluster
    module_num <- gsub(levelA_res$Module_RTclust, pattern = "_[0-9]{1,}", replacement = "")
     
    #store data frame linking Module_RTcluster with the associated mz values in the module       
    #changed numeric references to named references
    levelA_res_all <- levelA_res[, c('Module_RTclust', 'mz')] 
    
    #store information relating to original modules from level1A processing
    orig_module_labels <- levelA_res$Module_RTclust
    
    #overwrite Module_RTclust information with just the number number
    levelA_res_all$Module_RTclust <- module_num
    
    #calculate the mass defect, which is the difference between the measured mz and the nominal mz of the measured mass
    mzdefect <- 1 * ((levelA_res$mz - floor(levelA_res$mz)))
    
    #determine the density of mass defect values across the level1A_res features
    d1 <- density(mzdefect, bw = "nrd", from = min(mzdefect), 
        to = (0.01 + max(mzdefect)))

    #distribution of mass defect values (may be very informative if )
    plot(d1$x, d1$y, xlab = 'mass defect (accurate mz - nominal mz)', ylab = 'density', main = 'distribution of mass defect values in input dataset')
    
    t1 <- d1$x
    
    #group identification number
    gid <- paste("ISgroup", dim(dataA)[1], sep = "")
    
    levelA_res2 <- data.frame('mzdefect' = mzdefect, levelA_res)
    
    #updated to ensure the user-defined mass_defect_window width is applied for splitting mass defects in to groups
    #divide mass defects in to groups based on bins that span the range mass_defect_window
    #e.g. [0-0.009, 0.01-0.0199, 0.02-0.0299 ... 0.09-1]
    massdefect_cor_groups <- sapply(list(myData1 = levelA_res2), function(x) {
      split(x, cut(levelA_res2$mzdefect, breaks = seq(0, 1, mass_defect_window))) })
          
    #determine which mzdefect bins have an observation stored
    if (length(which(levelA_res2$mzdefect == 0) > 0)) {
        massdefect_cor_groups[[length(massdefect_cor_groups) + 1]] <- levelA_res2[which(levelA_res2$mzdefect == 0), ]
    }
    
    #set up empty object for capturing outputs
    diffmatB = {}
    
    #for each mass defect group, 
    diffmatB = lapply(1:length(massdefect_cor_groups), function(gnum) {
      
      cur_group <- {}
          
      if (dim(massdefect_cor_groups[[gnum]])[1] > 0) {
        
        #create an identifier for the mass defect processing group
        #structure is: 'ISgroup'_(module_number_pt1)_(module_number_pt2)_(massdefect_group)
        ISgroup <- paste("ISgroup", massdefect_cor_groups[[gnum]]$Module_RTclust, gnum, sep = "_")
            
        # print(massdefect_cor_groups[[gnum]])
        cur_group <- as.data.frame(massdefect_cor_groups[[gnum]])
        cur_group <- data.frame('ISGroup' = ISgroup,cur_group)
            
        if (length(cur_group) > 0) {
          cur_group <- cur_group[order(cur_group$mz, cur_group$time), ]

          if (length(diffmatB) < 1) {
            # diffmatB<-cur_group
          } else if (dim(cur_group)[1] > 0) {
            # diffmatB<-rbind(diffmatB,cur_group)
          }
        }else {
          print(gnum)
        }
      }
    
      return(cur_group)
    
    })
        
    #generate dataframe from results (by default, drops mass defect groups where no peaks found)
    diffmatB <- ldply(diffmatB, rbind)
    
    if (dim(diffmatB)[1] < dim(dataA)[1]) {
      diffmatC <- levelA_res2[-which(levelA_res$mz %in% diffmatB$mz), ]
                
      if (nrow(diffmatC) > 0) {
        t1 <- table(diffmatB[, 1])
        isop_last <- paste("ISgroup_", diffmatC$Module_RTclust, "_", 1, sep = "")
                  
        diffmatC <- cbind(isop_last, diffmatC)
        colnames(diffmatC) <- colnames(diffmatB)
                  
        diffmatD <- rbind(diffmatB[, c(1:10)], diffmatC[, c(1:10)])
                  
        rm(diffmatC)
        rm(diffmatB)
      } else {
        diffmatD <- diffmatB  #[,c(1:10)]
        rm(diffmatB)
      }
    } else {
      diffmatD <- diffmatB  #[,c(1:10)]
      rm(diffmatB)
    }

    rm(levelA_res2)
    diffmatD <- as.data.frame(diffmatD)
    diffmatD$mz <- as.numeric(diffmatD$mz)
    
    #re-order the diffmatD object by mz, then by retention time
    diffmatD <- diffmatD[order(diffmatD$mz, diffmatD$time),]
    
    #re-order the levelA_res object by mz, then by retention time
    levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
    
    #combine the levelA_res and diffmatD matrices
    levelA_res1 <- cbind(diffmatD$ISGroup, levelA_res)
    
    ################################
    # when replicated mzid values are found in the dataset, it might be necessary to remove corresponding rows from the 
    # mean_int_vec object (else it could have more rows than in the filtered dataA object)
    
    
    #cleaner combination of columns through use of the mzid variable
    diffmatD$mzid = paste(round(diffmatD$mz,5), round(diffmatD$time,2), sep = '_')
    levelA_res$mzid = paste(round(levelA_res$mz,5), round(levelA_res$time,2), sep = '_')
    mean_int_vec$mzid = paste(round(mean_int_vec$dataA.mz, 5), round(mean_int_vec$dataA.time, 2), sep='_')
    
    colnames(mean_int_vec) = gsub('mean_int_vec', 'AvgIntensity', colnames(mean_int_vec))
    
    #merge data frames on mzid
    levelA_res1 = merge(x = diffmatD, y = levelA_res, by = 'mzid')
    
    #reset mzid label (it gets rounded during merge)
    levelA_res1$mzid = paste(round(levelA_res1$mz.x,5), round(levelA_res1$time.x,2), sep = '_')
    
    #merge levelA_res1 and mean_int_vec data frames (VERY IMPORTANT mzid INPUTS WERE ROUNDED TO ENSURE CORRECT MAPPING)
    levelA_res1 = merge(x = levelA_res1, y = mean_int_vec, by = 'mzid')
    
    #clean up levelA_res1 object
    colnames(levelA_res1) = gsub('.x', '', colnames(levelA_res1))
    levelA_res1 = levelA_res1[,c('ISGroup', 'Module_RTclust', 'mz', 'time', 'AvgIntensity', 'mzdefect')]
    colnames(levelA_res1) = gsub('mzdefect', 'MD', colnames(levelA_res1))
    colnames(levelA_res1) = gsub('ISGroup', 'ISgroup', colnames(levelA_res1))
   
    #combine intensity data frame with diffmatD and define column names
    isop_res_md = merge(x = diffmatD, y = mean_int_vec, by ='mzid')
    colnames(isop_res_md)[grep('mzdefect', colnames(isop_res_md))] = 'MD'
    colnames(isop_res_md)[grep('ISGroup', colnames(isop_res_md))] = 'ISgroup'
    
    #keep relevant columns from the merged object
    isop_res_md = isop_res_md[,c('mz', 'time', 'ISgroup', 'Module_RTclust', 'AvgIntensity', 'MD')]
    
    
    #original approach (not using mzid for merging)
    #original approach was based on ordering data frames by mz and then by time, then cbind together
    #isop_res_md <- cbind(diffmatD[, c(4, 5, 1, 3)], mean_int_vec$mean_int_vec, diffmatD[, c(2)])
    #colnames(isop_res_md) <- c("mz", "time", "ISgroup", 
    #    "Module_RTclust", "AvgIntensity", "MD")
    # MD <- diffmatD[, c(2)]
    # levelA_res1 <- cbind(levelA_res1[, c(1:4)], mean_int_vec, 
    #     MD)
    # rm(MD)
    #cnames <- colnames(levelA_res1)
    #cnames[1] <- "ISgroup"
    
    #IMPORTANT matrix defining which compounds have been mapped to which mz
    multiresmat <- merge(levelB_res, levelA_res1, by = "mz")
    
    #reorder the multiresmat columns
    multiresmat <- multiresmat[, c("mz", "time.x", 
        "MatchCategory", "theoretical.mz", "chemical_ID", 
        "Name", "Formula", "MonoisotopicMass", "Adduct", 
        "ISgroup", "Module_RTclust", "time.y", "AvgIntensity", 
        "MD")]
    
    #fix time.x colname
    colnames(multiresmat) = gsub('.x', '', colnames(multiresmat))
    
    #order multiresmat matrix by the Module_RT parameter defined in level1_A
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust), ]
    
    #cleaning up
    rm(levelB_res)
    rm(m1)
    
    #determine whether a mz value has one (Unique) or more (Multiple) matches to something in the user-defined database
    multiresmat$MatchCategory = "Multiple"
    
    t2 <- table(multiresmat$mz)
    same1 <- which(t2 == 1)
    uniquemz <- names(same1)
    multiresmat$MatchCategory[which(multiresmat$mz %in% uniquemz)] <- "Unique"

    #reorder the multiresmat by Module_RTclust and chemical_ID
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust, multiresmat$chemical_ID), ]
            
    #this looks to be a repeat of the steps given above and can therefore be ignored
    #dupmz <- multiresmat$mz[which(duplicated(multiresmat$mz) == TRUE)]
    #tablemz <- table(multiresmat$Module_RTclust, multiresmat$mz)
    #multiresmat$MatchCategory <- rep("Multiple", dim(multiresmat)[1])
    #multiresmat$MatchCategory[-which(multiresmat$mz %in% dupmz)] <- "Unique"
    
    #remove matches where no chemical ID existed in the original database (i.e. chemical_ID column was blank)
    if (length(which(multiresmat$chemical_ID == "-")) > 0) {
        multiresmat <- multiresmat[-which(multiresmat$chemical_ID =="-"), ]
    }
            
    #table indicating the link between the chemical_ID (formula) and mz values (from dataA)
    tall <- table(multiresmat$chemical_ID, multiresmat$mz)
    tall_mznames <- colnames(tall)

    #determine whether mz associated with more-than one compound in the user-defined database
    tall_checkmultichems <- apply(tall, 2, sum)
    tall_multichems <- tall[, which(tall_checkmultichems > 1)]
    names_tall_multichems <- colnames(tall_multichems)
    
    #determine whether compound in database was associated with more-than one mz
    tall_checkmultimz <- apply(tall, 1, sum)
    tall_multimzperchem <- tall[which(tall_checkmultimz > 1), ]
    names_tall_multimzperchem <- rownames(tall_multimzperchem)
    chemids <- names_tall_multimzperchem  #chemids[which(t1chemids>0)]
    
    #compounds associated with unique mz value
    tall_unimzperchem <- tall[which(tall_checkmultimz == 1), ]
    single_mz_chemids <- rownames(tall_unimzperchem)
    
    #set up empty objects for use below
    chem_score <- list()
    i <- 1
    del_mz <- {}
    multiresmat_filt <- {}
     
    #already labelled as time       
    #cnames <- colnames(dataA)
    #cnames[2] <- "time"
    dataA <- unique(dataA)
    
    #cnames <- colnames(multiresmat)
    #cnames[10] <- "ISgroup"
    #colnames(multiresmat) <- cnames
       
    #VERY IMPORTANT OBJECT FOR ISOTOPE ASSIGNMENT (all unique features from dataA, after clustering)     
    level_module_isop_annot <- levelA_res1
    
    ### this whole block is already achieved below using unique(mchemdata$chemical_ID)
    ## determines which chemical IDs were mapped to a specific compound from the user-defined database
    #chem_ids <- table(multiresmat$chemical_ID)
    #chem_ids_1 <- chem_ids[which(chem_ids >= 1)]
    #chem_ids_1 <- names(chem_ids_1)
    #chem_ids_2 <- chem_ids_1

    #unnecessary steps - column 'time' is already called time
    #cnames <- colnames(multiresmat)
    #cnames[2] <- "time"
    #colnames(multiresmat) <- cnames
    
    #clean up
    rm(levelA_res)
    rm(m2)
    
    #save multiresmat as mchemdata (a data frame detailing which compounds and mz were linked together)
    mchemdata <- multiresmat
    rm(multiresmat)
    
    #remove rows where times from different arms of the annotation workflow don't map correctly?
    bad_rows <- which(abs(mchemdata$time - mchemdata$time.y) > 0)
    if (length(bad_rows) > 0) {
        mchemdata <- mchemdata[-bad_rows, ]
    }
    
    ## chemids corresponds to 'Formula_X' mapped to each mz value
    chemids <- unique(mchemdata$chemical_ID)
    
    chemscoremat <- {}
    
    #define the groups in which isotopes will be searched
    chemids_split = as.list(as.character(chemids))
    names(chemids_split) = as.character(chemids)
    
    #tidy-up mchemdata object
    mchemdata$mz <- round(as.numeric(as.character(mchemdata$mz)), 5)
    mchemdata$time <- round(as.numeric(as.character(mchemdata$time)), 1)
    mchemdata$MD <- round(as.numeric(as.character(mchemdata$MD)), 3)
    mchemdata$AvgIntensity <- round(as.numeric(as.character(mchemdata$AvgIntensity)), 1)
    mchemdata$time.y <- round(as.numeric(as.character(mchemdata$time.y)), 1) #why keep this?
   
    if (max_diff_rt >= 9999) {
        module_num <- gsub(mchemdata$Module_RTclust, 
          pattern = "_[0-9]{1,}", replacement = "")
        module_num <- paste(module_num, "_0", sep = "")
        mchemdata$Module_RTclust <- module_num
    }
    
    write.csv(mchemdata, file = "Stage2.csv", row.names = FALSE)
    
    save(list = c("outloc", "adduct_weights", "boostIDs", 
        "pathwaycheckmode", "adduct_table", "max_diff_rt", 
        "corthresh", "filter.by", "max.rt.diff", 
        "max_isp", "min_ions_perchem", "max.mz.diff", 
        "db_name", "allsteps", "redundancy_check", 
        "num_sets"), file = "tempobjects.Rda")
    
    save(list = c("mchemdata", "chemids", "adduct_table", 
        "mzid", "max_diff_rt", "isop_res_md", "corthresh", 
        "level_module_isop_annot", "chemids_split", 
        "corthresh", "max.mz.diff", "outloc", "num_sets", 
        "db_name", "num_nodes", "num_sets", "adduct_weights", 
        "filter.by", "max.rt.diff", "max_isp", "MplusH.abundance.ratio.check", 
        "mass_defect_window", "mass_defect_mode", 
        "allsteps"), file = "step1_results.Rda")
    
    rm(mchemdata)
    rm(chemids)
    rm(mzid)
    try(rm(global_cor), silent = TRUE)
    rm(isop_res_md)
    rm(level_module_isop_annot)
    rm(dataA)
    
    rm(tall_unimzperchem)
    rm(tablemz)
    rm(levelB_res2)
    rm(levelA_res1)
    
    #rm(list = ls())
    load("tempobjects.Rda")
    load("step1_results.Rda")
    
  } else {
    print("Status 1: Skipping step 1.")
    print("Status 2: Using existing step1_results.Rda file.")
    
    allsteps_temp <- allsteps
    load("tempobjects.Rda")
    allsteps <- allsteps_temp
  }

  #############################################################################
  ######################## BEGIN STAGE 3 OF PROCESSING ########################
  ## ISOTOPES ARE ASSIGNED IN THIS STEP
  
  #############################################################################
  
  if (allsteps == TRUE) {
    
    print("Status 3: Calculating scores for individual chemicals/metabolites")
    
    if (length(chemids_split) > num_nodes) {
      
      cl <- makeSOCKcluster(num_nodes)
          
      clusterEvalQ(cl, "multilevelannotationstep2")
      clusterExport(cl, "multilevelannotationstep2")
      clusterEvalQ(cl, "library(Rdisop)")
      clusterEvalQ(cl, "library(plyr)")
      clusterEvalQ(cl, "library(enviPat)")
      clusterExport(cl, "get_chemscorev1.6.71_custom")
      clusterExport(cl, "getMolecule")
      clusterExport(cl, "ldply")
      clusterExport(cl, "get_confidence_stage2")
      clusterExport(cl, "check_element")
      clusterExport(cl, "group_by_rt_histv2")
      clusterExport(cl, "adduct_table")
      clusterExport(cl, "adduct_weights")
      clusterExport(cl, "outloc")
      
      parLapply(cl, 1:length(chemids_split), function(arg1) {
        
        cur_fname <- paste(outloc, "/stage2/chem_score", arg1, ".Rda", sep = "")
        check_if_exists <- suppressWarnings(try(load(cur_fname)))
        
        if (is(check_if_exists, "try-error")) {
          multilevelannotationstep2(outloc1 = outloc, list_number = arg1)
        } else {
          print(paste("List ", arg1, " already exists.", sep = ""))
        }
        
      })
      
      # max.time.diff=max.rt.diff,filter.by=filter.by,max_isp=max_isp,numnodes=1,MplusH.abundance.ratio.check=MplusH.abundance.ratio.check,mass_defect_window=mass_defect_window,mass_defect_mode=mass_defect_mode
      
      stopCluster(cl)
      
  } else {
                
    for (arg1 in 1:num_sets) {
      multilevelannotationstep2(outloc1 = outloc, list_number = arg1)
    }
    
  }
        
  setwd(outloc)
  
  #clean environment
  rm(list = ls())
  
  try(rm(hmdbCompMZ), silent = TRUE)
  
  #load in earlier results
  load("tempobjects.Rda")

  print("Status 4: Pathway evaluation")
  
  annotres <- multilevelannotationstep3(outloc = outloc, adduct_weights = adduct_weights, 
                                        boostIDs = boostIDs, pathwaycheckmode = pathwaycheckmode)
  
  setwd(outloc)
  
  rm(list = ls())
  try(rm(hmdbCompMZ), silent = TRUE)
  try(rm(hmdbCompMZ, env = .GlobalEnv), silent = TRUE)
  
  load("tempobjects.Rda")

  print("Status 5: Assigning confidence levels")
  
  annotresstage4 <- multilevelannotationstep4(outloc = outloc, max.rt.diff = max_diff_rt, filter.by = filter.by,
                                              adduct_weights = adduct_table, min_ions_perchem = min_ions_perchem,
                                              boostIDs = boostIDs, max_isp = max_isp, max.mz.diff = max.mz.diff)
  
  rm(list = ls())
  load("tempobjects.Rda")
 
  try(rm(hmdbAllinf, env = .GlobalEnv), silent = TRUE)
  
  
  
  if (redundancy_check == TRUE) {
    
    print("Status 6:Redundancy filtering")
    
    rm(list = ls())
    
    load("tempobjects.Rda")
    
    suppressWarnings(annotresstage5 <- multilevelannotationstep5(outloc = outloc, 
      max.rt.diff = max_diff_rt, filter.by = filter.by, 
      adduct_weights = adduct_weights, min_ions_perchem = min_ions_perchem, 
      boostIDs = boostIDs, max_isp = max_isp, 
      db_name = db_name, max.mz.diff = max.mz.diff))
    
    # print('Stage 5 confidence level distribution for unique
    # chemical/metabolite IDs')
    # print(table(annotresstage5$Confidence[-which(duplicated(annotresstage5$chemical_ID)==TRUE)]))
    
  } else {
    
    annotresstage5 <- {
    }
  }
}
        
    } else {
        
        annotres <- {
        }  #mchemdata
    }
    
    
    
}
    
    
    print("################")
    print("Final: Processing complete")
    
    
    print("Output files description:")
    print("Stage 1 includes clustering of features based on intensity and retention time without annotations")
    print("Stage 2 includes clustering results along with simple m/z based database matching")
    
    if (allsteps == TRUE) {
        print("Stage 3 includes scores for annotations assigned in stage 2 based on multiple criteria")
        print("Stages 4 and 5 include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
    }
    
    suppressWarnings(sink(file = NULL))
    
    setwd(outloc)
    sink(file = "Readme.txt")
    print("Output files description:")
    print("Stage1.csv includes clustering of features based on intensity and retention time without annotations")
    print("DBresults/Stage2.csv includes clustering results along with simple m/z based database matching")
    
    if (allsteps == TRUE) {
        print("DBresults/Stage3.csv includes scores for annotations assigned in stage 2 based on multiple criteria")
        print("DBresults/Stage4.csv and DBresults/Stage5.csv include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
    }
    suppressWarnings(sink(file = NULL))
    
    if (allsteps == TRUE) {
        # return(list('stage4'=annotresstage4,'stage5'=annotresstage5))
    }
    
}
