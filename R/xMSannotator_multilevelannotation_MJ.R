dataA = temp
max.mz.diff = 10
max.rt.diff = 10
cormethod = "pearson"
num_nodes = 2
queryadductlist = c("all") 
gradienttype = "Acetonitrile" 
mode = "pos"
outloc
db_name = "Custom"
adduct_weights = NA
num_sets = 3000
allsteps = TRUE
corthresh = 0.7
NOPS_check = TRUE
customIDs = NA
missing.value = NA
deepsplit = 2
networktype = "unsigned"
minclustsize = 10
module.merge.dissimilarity = 0.2
filter.by = c("M+H")
redundancy_check = TRUE
min_ions_perchem = 1
biofluid.location = NA
origin = NA
status = NA
boostIDs = NA
max_isp = 5
MplusH.abundance.ratio.check = FALSE
customDB = NA
HMDBselect = "union"
mass_defect_window = 0.01
mass_defect_mode = "pos"
dbAllinf = NA
pathwaycheckmode = "pm"
outloc = 'C:/Users/jonesmar/Documents/Temp_data/xMSannotator_suspectScreening'

###### IMPORTANT FILE!!!!! ######## read in updated table of adducts from .csv file
adduct_table = read.csv(file.path(getwd(), '/xMSannotator/adducts_xMSannotator_enviPat.csv', fsep=''), header = T, stringsAsFactors = FALSE)
adduct_table$Weight = 1


options(warn = -1)
allowWGCNAThreads(nThreads = num_nodes)
dataA <- as.data.frame(dataA)
dataA$mz <- round(as.numeric(dataA$mz), 5)
dataA$time <- round(as.numeric(dataA$time), 1)

dataA[, -c(1:2)] <- round(dataA[, -c(1:2)], 1)

if (is.na(module.merge.dissimilarity) == TRUE) {
  
  module.merge.dissimilarity = 1 - corthresh
}

if (is.na(customIDs) == FALSE) {
  
  customIDs = as.data.frame(customIDs)
}

data(adduct_table)
max_diff_rt <- max.rt.diff
cutheight = 1 - corthresh
time_step = 1
step1log2scale = FALSE

#does absolutely nothing
if (FALSE) {
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "2M\\+2H", replacement = "M+H")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "3M\\+3H", replacement = "M+H")
  
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "2M\\+2Na", replacement = "M+Na")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "3M\\+3Na", replacement = "M+Na")
  
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "2M\\+2K", replacement = "M+K")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, 
                              pattern = "3M\\+3K", replacement = "M+K")
}

adduct_table <- unique(adduct_table)
suppressWarnings(dir.create(outloc))
setwd(outloc)

#refined adduct table generated
if (queryadductlist == "all" & mode == "pos") {
  adduct_names <- adduct_table$Adduct[(adduct_table$Type == "S" & adduct_table$Mode == "positive") | (adduct_table$Type == gradienttype & adduct_table$Mode == "positive")]
  adduct_table <- adduct_table[which(adduct_table$Adduct %in% adduct_names), ]
}

adduct_names <- unique(adduct_names)
res_list <- new("list")
db_name_list = db_name

outloc_allres <- outloc

{
  l1 <- list.files(outloc_allres) #list all files in output directory
  check_step1 <- which(l1 == "step1_results.Rda") #look for file in out directory
  
  if (length(check_step1) < 1) {
  
    print("Status 1: Computing modules using WGCNA")
    check_levelA <- which(l1 == "xMSannotator_levelA_modules.Rda")
    
    if (is.na(missing.value) == FALSE) {
      dataA <- replace(as.matrix(dataA), which(dataA == 
                                                 missing.value), NA)
      dataA <- as.data.frame(dataA)
    }
    
    dataA <- unique(dataA)
    
    mzid <- paste(dataA$mz, dataA$time, sep = "_") #generate vector containing mz_id pairs
    
    if (length(which(duplicated(mzid) == TRUE)) > 
        0) {
      dataA <- dataA[-which(duplicated(mzid) == 
                              TRUE), ]
    }
    
    mzid <- paste(dataA$mz, dataA$time, sep = "_")
    
    system.time(global_cor <- WGCNA::cor(t(dataA[,-c(1:2)]), nThreads = num_nodes, 
                                         method = cormethod, use = "p"))
    
    global_cor <- round(global_cor, 2)
    
    dataA <- unique(dataA)
    setwd(outloc_allres)
    mzid <- paste(dataA$mz, dataA$time, sep = "_")
    
    colnames(global_cor) <- mzid #assign column and row names to the global correction matrix
    rownames(global_cor) <- mzid
    
    save(global_cor, file = "global_cor.Rda") #save the output
    
    mycl_metabs = NA

    #m = m[rowSums(is.na(m)) != ncol(m), ]
    #m = m[,colSums(is.na(m)) != nrow(m)]
    
    
    
    
    
    glc = global_cor
    
    ####### FILLIING MISSING VALUES ADDED BY ME #######
    
    #permute random values between -1 and 1 for NA values and store output matrix
    capture = list()
    for(i in 1:12){
      while(i < 12){
        glc_temp = glc
        glc_temp[is.na(glc_temp)] = runif(n = length(which(is.na(glc_temp) == TRUE)), -1, 1) 
        capture[[i]] = glc_temp
        i = i + 1
      }
    }
    
    #export matrices from list and convert to array of dimensions (#rows, #cols, #matrices)
    arr <- array(unlist(capture), dim = c(dim(capture[[1]]),length(capture)))
    
    #calculate 'average' correlation coefficient (for values that were not NA earlier, this will remain unchanged)
    result = rowMeans(arr, dim = 2) 
    
    colnames(result) <- mzid #assign column and row names to the global correction matrix
    rownames(result) <- mzid
    
    global_cor = result
    rm(result)

    hr = flashClust(as.dist(1 - global_cor), method = "complete")
    dissTOMCormat <- (1 - global_cor)
    
    a1 <- apply(global_cor, 1, function(x) {
      length(which(x > corthresh))
    }) #determine how many other features a given feature is correlated with
    
    fname_c <- paste("NumberOfCorrelationsPerFeature_cor", 
                     corthresh, ".csv", sep = "")
    
    write.csv(a1, file = fname_c)
    
    rm(global_cor)
    
    mycl_metabs <- cutreeDynamic(hr, distM = dissTOMCormat, 
                                 deepSplit = 1, minClusterSize = minclustsize, 
                                 pamRespectsDendro = FALSE, pamStage = TRUE, 
                                 verbose = 0)
    
    clustmethod = "WGCNA"
    
    if (length(check_levelA) < 1) {
      
      if (clustmethod == "WGCNA") {
        
        
        levelA_res <- get_peak_blocks_modulesvhclust(dataA = dataA, 
                                                     simmat = NA, adjacencyfromsimilarity = FALSE, 
                                                     time_step = time_step, max.rt.diff = max_diff_rt, 
                                                     outloc, column.rm.index = NA, cor.thresh = NA, 
                                                     deepsplit = deepsplit, minclustsize = minclustsize, 
                                                     cutheight = cutheight, cormethod = cormethod, 
                                                     networktype = networktype, num_nodes = num_nodes, 
                                                     step1log2scale = step1log2scale, mycl_metabs = mycl_metabs)
        
        
        setwd(outloc_allres)
        levelA_res <- levelA_res[, -c(1:4)]
        save(levelA_res, file = "xMSannotator_levelA_modules.Rda")
        
        
        
      } else {
        
        if (clustmethod == "graph") {
          
          levelA_res <- get_peak_blocks_graph(dataA, 
                                              simmat = global_cor, adjacencyfromsimilarity = TRUE, 
                                              time_step = 3, max.rt.diff = max_diff_rt, 
                                              outloc, column.rm.index = NA, cor.thresh = NA, 
                                              cormethod = cormethod, networktype = networktype, 
                                              num_nodes = num_nodes)
          
        } else {
          
          
          if (clustmethod == "hclust") {
            levelA_res <- get_peak_blocks_hclust(dataA, 
                                                 time_step = time_step, max.rt.diff = max_diff_rt, 
                                                 outloc, column.rm.index = NA, cor.thresh = NA, 
                                                 deepsplit = deepsplit, minclustsize = minclustsize, 
                                                 cutheight = cutheight, cormethod = cormethod)
          }
          
        }
        
      }
      
      setwd(outloc_allres)
      
      
      write.csv(levelA_res, file = "Stage1.csv", 
                row.names = FALSE)
      # write.table(levelA_res,file='Stage1.txt',sep='\t',row.names=FALSE)
      
    } else {
      setwd(outloc_allres)
      load("xMSannotator_levelA_modules.Rda")
      
      # alldegrees<-levelA_res[,c(1:4)]
      # levelA_res<-levelA_res[,-c(1:4)]
      
      
      
    }
    
    setwd(outloc)
    
    levelA_res <- levelA_res[order(levelA_res$mz, 
                                   levelA_res$time), ]
    
    levelA_res <- levelA_res[, c(1:3)]
    
    dataA <- dataA[order(dataA$mz, dataA$time), ]
    
    mean_int_vec <- apply(dataA[, -c(1:2)], 1, function(x) {
      mean(x, na.rm = TRUE)
    })
    
    dataA <- dataA[, c(1:2)]
    
    if (queryadductlist == "all" & mode == "pos") {
      
      adduct_names <- adduct_table$Adduct[(adduct_table$Type == 
                                             "S" & adduct_table$Mode == "positive") | 
                                            (adduct_table$Type == gradienttype & adduct_table$Mode == 
                                               "positive")]
      
      adduct_table <- adduct_table[which(adduct_table$Adduct %in% 
                                           adduct_names), ]
      
    } else {
      if (queryadductlist == "all" & mode == "neg") {
        
        adduct_names <- adduct_table$Adduct[(adduct_table$Type == 
                                               "S" & adduct_table$Mode == "negative") | 
                                              (adduct_table$Type == gradienttype & 
                                                 adduct_table$Mode == "negative")]
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% 
                                             adduct_names), ]
      } else {
        
        adduct_names <- adduct_table$Adduct[which(adduct_table$Adduct %in% 
                                                    queryadductlist)]
        
        adduct_table <- adduct_table[which(adduct_table$Adduct %in% 
                                             adduct_names), ]
        
      }
    }
    
    
    adduct_names <- unique(adduct_names)
    
    if (db_name == "Custom") {
     
      data(adduct_table)
      
      inputmassmat <- custDB
      
      if (length(which(duplicated(inputmassmat[,1]) == TRUE)) > 0) {
        inputmassmat <- inputmassmat[-which(duplicated(inputmassmat[,1]) == TRUE), ]
      }
      
      mz_search_list <- {
      }
      
      mz_search_list <- lapply(1:dim(inputmassmat)[1], function(m) {
                                 
         # print(head(adduct_table))
         adduct_names <- as.character(adduct_names)
         # print(adduct_names)
         
         mz_search_list <- get_mz_by_monoisotopicmass(monoisotopicmass = as.numeric(as.character(inputmassmat[m,4])), 
                                                      dbid = inputmassmat[m,1], name = as.character(inputmassmat[m,2]), 
                                                      formula = as.character(inputmassmat[m,3]), 
                                                      queryadductlist = adduct_names, 
                                                      adduct_table = adduct_table)
         
         return(mz_search_list)
       })
      
      mz_search_list <- ldply(mz_search_list, 
                              rbind)
      save(mz_search_list, file = "mz_search_list.Rda")
      chemCompMZ <- mz_search_list
      rm(mz_search_list)
      rm(customDB)
      rm(customDB, envir = .GlobalEnv)
      
    }
      
    data(adduct_table)
    
    if (is.na(adduct_weights) == TRUE) {
      data(adduct_weights)
      adduct_weights1 <- matrix(nrow = 2, ncol = 2, 
                                0)
      adduct_weights1[1, ] <- c("M+H", 1)
      adduct_weights1[2, ] <- c("M-H", 1)
      adduct_weights <- as.data.frame(adduct_weights1)
      colnames(adduct_weights) <- c("Adduct", "Weight")
      adduct_weights <- as.data.frame(adduct_weights)
      
    }
    
    chemCompMZ$Name <- gsub("([-:();])|[[:punct:]]", 
                            "\\1", chemCompMZ$Name)
    chemCompMZ$Formula <- gsub("([-:();])|[[:punct:]]", 
                               "\\1", chemCompMZ$Formula)
    
    
    cnames <- colnames(chemCompMZ)
    
    cnames[2] <- "chemical_ID"
    colnames(chemCompMZ) <- cnames
    
    
    randsetsize <- 1000
    if (dim(dataA)[1] < 1000) {
      
      randsetsize <- dim(dataA)[1]
    }
    
    chemCompMZ <- as.data.frame(chemCompMZ)
    
    setwd(outloc)
    l1 <- list.files(outloc)
    check_levelB <- which(l1 == "xMSannotator_levelB.Rda")
    
    save(chemCompMZ, file = "chemCompMZ.Rda")
    
    if (length(check_levelB) < 1) {
      
      cl <- makeSOCKcluster(num_nodes)
      
      clusterEvalQ(cl, library(XML))
      clusterEvalQ(cl, library(R2HTML))
      clusterEvalQ(cl, library(RCurl))
      clusterEvalQ(cl, library(SSOAP))
      clusterEvalQ(cl, library(limma))
      
      clusterEvalQ(cl, library(plyr))
      
      clusterEvalQ(cl, "processWSDL")
      clusterEvalQ(cl, library(png))
      clusterExport(cl, "Annotationbychemical_IDschild_multilevel")
      
      clusterExport(cl, "Annotationbychemical_IDschild")
      
      clusterExport(cl, "find.Overlapping.mzs")
      clusterExport(cl, "find.Overlapping.mzsvparallel")
      
      clusterExport(cl, "overlapmzchild")
      clusterExport(cl, "getVenn")
      
      
      if (length(which(duplicated(chemCompMZ$Formula) ==TRUE)) > 0) {
      
        
        if (db_name == "Custom") {
          # chemCompMZ_unique_formulas<-chemCompMZ
          
          # chemCompMZ_unique_formulas<-as.data.frame(chemCompMZ_unique_formulas)
          
          # chemCompMZ_unique_formulas$mz<-as.numeric(as.character(chemCompMZ_unique_formulas$mz))
          
          
          chemCompMZ$mz <- as.numeric(as.character(chemCompMZ$mz))
          
          
        }
        
        
        #THIS PROCEDURE REMOVES COMPOUNDS WITH OVERLAPPING FORMULAE!!!
        
        #list the unique compounds in the chemCompMZ object
        chemCompMZ_unique_formulas <- chemCompMZ[-which(duplicated(chemCompMZ$Formula) == TRUE), ]

        
        #extract all associated adduct forms for compounds that are deemed unique (based on Formula)
        chemCompMZ_unique_formulas <- rbind(chemCompMZ_unique_formulas, 
                                            chemCompMZ[which(chemCompMZ$chemical_ID %in% 
                                            chemCompMZ_unique_formulas$chemical_ID), ])
        
        #ensure that all entries are unique
        chemCompMZ_unique_formulas <- unique(chemCompMZ_unique_formulas)
        #re-order based on chemical ID
        chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[order(chemCompMZ_unique_formulas$chemical_ID),]
      
      } else {
        
        chemCompMZ_unique_formulas <- chemCompMZ
      }
      
      
      # save(chemCompMZ_unique_formulas,file='chemCompMZ_unique_formulas.Rda')
      save(chemCompMZ, file = "chemCompMZ.Rda")
      rm(chemCompMZ)
      
      #generate table
      formula_table <- table(chemCompMZ_unique_formulas$Formula) #summarises number of instances of each formula
      uniq_formulas <- names(formula_table) #gets formula from table colname
      formula_ID <- paste("Formula", seq(1:length(uniq_formulas)), 
                          sep = "_") #unique identifier for each formula
      
      formula_id_mat <- cbind(formula_ID, uniq_formulas)
      formula_id_mat <- as.data.frame(formula_id_mat) #link each formula_ID to each formula in dataframe
      colnames(formula_id_mat) <- c("Formula_ID", 
                                    "Formula")
      
      chemCompMZ_unique_formulas <- merge(chemCompMZ_unique_formulas, 
                                          formula_id_mat, by = "Formula")
      
      chemCompMZ_unique_formulas$chemical_ID <- chemCompMZ_unique_formulas$Formula_ID #overwrite chemicalID with formulaID
      
      chemCompMZ_unique_formulas <- chemCompMZ_unique_formulas[,c("mz", "chemical_ID", "Name", "Formula", 
                                                                 "MonoisotopicMass", "Adduct", "AdductMass")]
      
      chemCompMZ_unique_formulas <- as.data.frame(chemCompMZ_unique_formulas)
      
      chemIDs <- unique(chemCompMZ_unique_formulas$chemical_ID) #unique(chemCompMZ$Formula) #
      
      s1 <- seq(1, length(adduct_names))
      print("Status 2: Mapping m/z to metabolites:")
      
      adduct_names <- as.character(adduct_names)
      
      l2 <- parLapply(cl, s1, Annotationbychemical_IDschild, 
                      dataA = dataA, queryadductlist = c(adduct_names), 
                      adduct_type = c("S", gradienttype), adduct_table = adduct_table, 
                      max.mz.diff = max.mz.diff, outloc = outloc, 
                      keggCompMZ = chemCompMZ_unique_formulas, 
                      otherdbs = FALSE, otherinfo = FALSE, num_nodes = num_nodes)
      
      
      
      stopCluster(cl)
      
      # print(format(object.size(l2),unit='auto'))
      # save(l2,file='l2.Rda')
      
      rm(chemCompMZ)
      levelB_res <- {
      }
      for (j in 1:length(l2)) {
        if (length(l2[[j]]) > 1) {
          levelB_res <- rbind(levelB_res, l2[[j]])
        }
      }
      
      MatchCategory = rep("Multiple", dim(levelB_res)[1])
      
      
      levelB_res$mz <- as.numeric(as.character(levelB_res$mz))
      
      levelB_res$time <- as.numeric(as.character(levelB_res$time))
      
      levelB_res <- as.data.frame(levelB_res)
      
      levelB_res <- cbind(levelB_res, MatchCategory)
      
      print("DB matches")
      levelB_res <- unique(levelB_res)
      print(dim(levelB_res))
      
      uniq_formula <- as.character(unique(levelB_res$Formula))
      
      bad_formula <- which(is.na(uniq_formula) == 
                             TRUE)
      if (length(bad_formula) > 0) {
        uniq_formula <- uniq_formula[-c(bad_formula)]
      }
      
      cl <- makeSOCKcluster(num_nodes)
      
      
      clusterExport(cl, "check_golden_rules")
      clusterExport(cl, "check_element")
      
      levelB_res_check <- parLapply(cl, 1:length(uniq_formula), 
                                    function(j, uniq_formula, NOPS_check) {
                                      
                                      curformula <- as.character(uniq_formula[j])
                                      return(check_golden_rules(curformula, 
                                                                NOPS_check = NOPS_check))
                                      
                                    }, uniq_formula = uniq_formula, NOPS_check = NOPS_check)
      stopCluster(cl)
      
      #check the elemental composition is appropriate
      levelB_res_check2 <- ldply(levelB_res_check, 
                                 rbind)
      
      #keep only those molecular formula that passed the element check
      levelB_res_check3 <- levelB_res_check2[which(levelB_res_check2[, 
                                                                     2] == 1), ]
      
      #filter levelB_res to keep all rows where the element check was passed
      levelB_res <- levelB_res[which(levelB_res$Formula %in% 
                                       as.character(levelB_res_check3[, 1])), 
                               ]
      
      water_adducts <- c("M+H-H2O", "M+H-2H2O", 
                         "M-H2O-H")
      
      water_adduct_ind <- which(levelB_res$Adduct %in% 
                                  water_adducts)
      
      cl <- makeSOCKcluster(num_nodes)
      
      
      clusterExport(cl, "check_element")
      
      levelB_res2_bkup = levelB_res2
      
      
      if (length(water_adduct_ind) > 0) {
        
        #subset dataframe to keep only those containing a water adduct
        levelB_res2 <- levelB_res[c(water_adduct_ind),]
        
        #subset dataframe to keep only those not containing water adduct
        levelB_res <- levelB_res[-c(water_adduct_ind),]
        
        #set up index from 1 to number of features associated with water adduct
        sind1 <- seq(1:dim(levelB_res2)[1])
        
        #in parallel, perform a validity check for element 'O'
        levelB_res_check3 <- parLapply(cl, sind1,function(j, levelB_res2){
          adduct <- as.character(levelB_res2$Adduct[j])
          curformula <- as.character(levelB_res2$Formula[j])
          numoxygens <- check_element(curformula,"O")
          
          if (numoxygens > 0) {
            bool_check <- 1
          } else {
            bool_check <- 0
          }
                                         
          res <- cbind(curformula, bool_check)
          res <- as.data.frame(res)
          return(res)}, levelB_res2)
        
        #turn element check output in to a dataframe
        levelB_res_check4 <- ldply(levelB_res_check3, 
                                   rbind)
        
        #set up empty object for capturing valid formulae
        valid_form <- {
        }
        
        #filter to retain only those with a valid number of oxygen atoms
        if (length(which(levelB_res_check4[, 2] == 
                         1)) > 0) {
          levelB_res_check4 <- levelB_res_check4[which(levelB_res_check4[, 
                                                                         2] == 1), ]
          
          
          valid_form <- which(levelB_res2$Formula %in% 
                                as.character(levelB_res_check4[, 1]))
        }
        
        #combine 'valid' water-adduct based features, with non-water adduct-based features
        if (length(valid_form) > 0) {
          levelB_res2 <- levelB_res2[valid_form, 
                                     ]
          levelB_res <- rbind(levelB_res, levelB_res2)
        }
        
      }
      # levelB_res<-levelB_res[,-c(1)]
      
      
      colnames(levelB_res) <- c("theoretical.mz", 
                                "chemical_ID", "Name", "Formula", "MonoisotopicMass", 
                                "Adduct", "AdductMass", "mz", "time", "MatchCategory")
      
      save(levelB_res, file = "xMSannotator_levelB.Rda")
      
      
      
    } else {
      
      print("Status 2: Using existing m/z mapping results:")
      load("xMSannotator_levelB.Rda")
    }
      
    levelA_res$mz <- round(levelA_res$mz, 5)
    levelB_res$mz <- round(levelB_res$mz, 5)
      
    last_col_ind <- dim(levelB_res)[2] 
    
    levels_A <- levels(levelA_res$Module_RTclust)  
    levels_A <- table(levelA_res$Module_RTclust) 
    levels_A <- names(levels_A)
    
    levelA_res <- levelA_res[order(levelA_res$mz, 
                                   levelA_res$time), ]
    
    
    if (length(which(duplicated(mzid) == TRUE)) > 
        0) {
      levelA_res <- levelA_res[-which(duplicated(mzid) == 
                                        TRUE), ]
      
      dataA <- dataA[-which(duplicated(mzid) == 
                              TRUE), ]
    }
    
    levelA_res <- levelA_res[order(levelA_res$Module_RTclust),]
      
    module_num <- gsub(levelA_res$Module_RTclust, 
                       pattern = "_[0-9]{1,}", replacement = "")
    
    levelA_res_all <- levelA_res[, c(1:2)]
    
    orig_module_labels <- levelA_res$Module_RTclust
    
    levelA_res_all$Module_RTclust <- module_num
    
    mzdefect <- 1 * ((levelA_res$mz - floor(levelA_res$mz)))
    
    d1 <- density(mzdefect, bw = "nrd", from = min(mzdefect), 
                  to = (0.01 + max(mzdefect)))
    
    t1 <- d1$x
    
    
    gid <- {
    }
    
    gid <- paste("ISgroup", dim(dataA)[1], sep = "")
    
    if (length(which(levelA_res2$mzdefect == 0) > 
               0)) {
      massdefect_cor_groups[[length(massdefect_cor_groups) + 
                               1]] <- levelA_res2[which(levelA_res2$mzdefect == 
                                                          0), ]
    }
    
    diffmatB <- {
    }
    
    # temp = ldply(c(1:length(massdefect_cor_groups)), function(x, massdefect_cor_groups){
    #   
    #   out = cbind(x, massdefect_cor_groups[[x]])
    #   
    # }, massdefect_cor_groups)
    
    diffmatB <- lapply(1:length(massdefect_cor_groups), function(gnum) {
                    cur_group <- {}
                    
                    if (dim(massdefect_cor_groups[[gnum]])[1] > 0) {
                      ISgroup <- paste("ISgroup", 
                                       massdefect_cor_groups[[gnum]]$Module_RTclust, 
                                       gnum, sep = "_")
                           
                       # print(massdefect_cor_groups[[gnum]])
                       cur_group <- as.data.frame(massdefect_cor_groups[[gnum]])
                       cur_group <- cbind(ISgroup, cur_group)
                           
                       if (length(cur_group) > 0) {
                         
                         cur_group <- cur_group[order(cur_group$mz, 
                                                      cur_group$time), ]
                         
                         
                         if (length(diffmatB) < 1) {
                           # diffmatB<-cur_group
                         } else {
                           
                           
                           if (dim(cur_group)[1] > 0) {
                             # diffmatB<-rbind(diffmatB,cur_group)
                           }
                         }
                         
                       }
                     } else {
                       print(gnum)
                     }
                     return(cur_group)
                   })

    diffmatB <- ldply(diffmatB, rbind)
    diffmatB <- as.data.frame(diffmatB)
    
    if (dim(diffmatB)[1] > 0) {
      cnames <- colnames(diffmatB)
      cnames[1] <- "ISGroup"
      
      colnames(diffmatB) <- cnames
    }
    
    if (dim(diffmatB)[1] < dim(dataA)[1]) {
      diffmatC <- levelA_res2[-which(levelA_res$mz %in% 
                                       diffmatB$mz), ]
      
      if (nrow(diffmatC) > 0) {
        t1 <- table(diffmatB[, 1])
        isop_last <- paste("ISgroup_", diffmatC$Module_RTclust, 
                           "_", 1, sep = "")
        
        diffmatC <- cbind(isop_last, diffmatC)
        colnames(diffmatC) <- colnames(diffmatB)
        
        
        diffmatD <- rbind(diffmatB[, c(1:10)], 
                          diffmatC[, c(1:10)])
        
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
    
    diffmatD <- diffmatD[order(diffmatD$mz, diffmatD$time),]
    
    levelA_res <- levelA_res[order(levelA_res$mz, 
                                   levelA_res$time), ]
    
    levelA_res1 <- cbind(diffmatD[, 1], levelA_res)
    
    isop_res_md <- cbind(diffmatD[, c(4, 5, 1, 3)], 
                         mean_int_vec, diffmatD[, c(2)])
    
    colnames(isop_res_md) <- c("mz", "time", "ISgroup", 
                               "Module_RTclust", "AvgIntensity", "MD")
    
    MD <- diffmatD[, c(2)]
    
    cnames <- colnames(levelA_res1)
    cnames[1] <- "ISgroup"
    
    colnames(levelA_res1) <- cnames
    
    multiresmat <- merge(levelB_res, levelA_res1, 
                         by = "mz")
    
    multiresmat <- multiresmat[, c("mz", "time.x", 
                                   "MatchCategory", "theoretical.mz", "chemical_ID", 
                                   "Name", "Formula", "MonoisotopicMass", "Adduct", 
                                   "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", 
                                   "MD")]
    
    colnames(multiresmat) <- c("mz", "time", "MatchCategory", 
                               "theoretical.mz", "chemical_ID", "Name", 
                               "Formula", "MonoisotopicMass", "Adduct", 
                               "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", 
                               "MD")
    
    rm(levelB_res)
    
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust),]
    
    rm(m1) 
    
    t2 <- table(multiresmat$mz)
    
    same1 <- which(t2 == 1)
    
    uniquemz <- names(same1)
    
    multiresmat$MatchCategory = rep("Multiple", dim(multiresmat)[1])
    
    multiresmat$MatchCategory[which(multiresmat$mz %in% 
                                      uniquemz)] <- "Unique"
    
    setwd(outloc)
    
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust, 
                                     multiresmat$chemical_ID), ]
    
    dupmz <- multiresmat$mz[which(duplicated(multiresmat$mz) == 
                                    TRUE)]
    
    #generate a matrix indicating the number of occurrences of each retention time cluster by mz pair.
    tablemz <- table(multiresmat$Module_RTclust, 
                     multiresmat$mz)
    
    multiresmat$MatchCategory <- rep("Multiple",dim(multiresmat)[1])
    
    multiresmat$MatchCategory[-which(multiresmat$mz %in% dupmz)] <- "Unique"
    
    #remove instances where formula is indicated by a hyphen
    if (length(which(multiresmat$chemical_ID == "-")) > 0) {
      multiresmat <- multiresmat[-which(multiresmat$chemical_ID == "-"), ]
    }
    
    #matrix indicating occurrence of each mz versus each formula
    tall <- table(multiresmat$chemical_ID, multiresmat$mz)
    
    #list of all mz values
    tall_mznames <- colnames(tall)
    
    #get instances where multiple formula per mz (sum > 1 means this is TRUE)
    tall_checkmultichems <- apply(tall, 2, sum)
    
    tall_multichems <- tall[, which(tall_checkmultichems > 1)]
    
    #get mz values associated with multiple formulae
    names_tall_multichems <- colnames(tall_multichems)
    
    #check which formulae are associated with multiple mz values
    tall_checkmultimz <- apply(tall, 1, sum)
    
    #get the instances where multiple mz per formula
    tall_multimzperchem <- tall[which(tall_checkmultimz > 1), ]
    
    #get instances where one mz matches to one formula
    tall_unimzperchem <- tall[which(tall_checkmultimz == 1), ]
    
    #get the identifier for formulae associated with multiple mz
    names_tall_multimzperchem <- rownames(tall_multimzperchem)
    
    #store names of formulae associated with multiple mz in object chemids
    chemids <- names_tall_multimzperchem
    
    #store names of formulae associated with just one mz
    single_mz_chemids <- rownames(tall_unimzperchem)
    
    chem_score <- new("list")
    i <- 1
    
    del_mz <- {
    }
    
    multiresmat_filt <- {
    }
    
    cnames <- colnames(dataA)
    
    cnames[2] <- "time"
    
    colnames(dataA) <- as.character(cnames)
    
    dataA <- unique(dataA)
    
    cnames <- colnames(multiresmat)
    cnames[10] <- "ISgroup"
    colnames(multiresmat) <- cnames
    
    #create object for isotope annotation from levelA_res1 object
    level_module_isop_annot <- levelA_res1
    
    #determine number of times each compound is found from the database
    chem_ids <- table(multiresmat$chemical_ID)
    #for those compounds found more than once, store their formula ID
    chem_ids_1 <- chem_ids[which(chem_ids >= 1)]
    chem_ids_1 <- names(chem_ids_1)
    
    chem_ids_2 <- chem_ids_1
    
    chemids <- chem_ids_2
    
    mchemdata <- multiresmat
    
    mchemdata <- as.data.frame(mchemdata)
    
    mchemdata$mz <- as.numeric(as.character(mchemdata$mz))
    mchemdata$time <- as.numeric(as.character(mchemdata$time))
    
    bad_rows <- which(abs(mchemdata$time - mchemdata$time.y) > 0)
    
    if (length(bad_rows) > 0) {
      
      mchemdata <- mchemdata[-bad_rows, ]
      
    }
    
    chemids <- unique(mchemdata$chemical_ID)
    chemids <- unique(chemids)
    
    chemscoremat <- {
    }
    
    if (length(chemids) > 10) {
      num_sets <- length(chemids)/2
      
    } else {
      
      num_sets <- 1
    }
    
    list_winsize <- num_sets
    list_size <- round(length(chemids)/list_winsize)
    
    #added by me
    list_size <- list_winsize
    
    list_size = 1
    
    if (length(chemids) > list_winsize) {
      
      g <- seq(1, length(chemids), list_size)
      g <- factor(g)
      chemids_split <- split(1:length(chemids), 
                             f = g)
      split_size <- 1:list_winsize
      
    } else {
      chemids_split <- split(1:length(chemids), 
                             f = length(chemids))
      split_size <- c(1)
    }
    
    num_sets = length(chemids_split)
    
    {
      mchemdata$mz <- round(as.numeric(as.character(mchemdata$mz)), 
                            5)
      mchemdata$time <- round(as.numeric(as.character(mchemdata$time)), 
                              1)
      mchemdata$MD <- round(as.numeric(as.character(mchemdata$MD)), 
                            3)
      mchemdata$mean_int_vec <- round(as.numeric(as.character(mchemdata$mean_int_vec)), 
                                      1)
      mchemdata$time.y <- round(as.numeric(as.character(mchemdata$time.y)), 
                                1)
    }
    
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
                  "num_sets", "chemids_split"), file = "tempobjects.Rda")
    
    save(list = c("mchemdata", "chemids", "adduct_table", 
                  "mzid", "max_diff_rt", "isop_res_md", "corthresh", 
                  "level_module_isop_annot", "chemids_split", 
                  "corthresh", "max.mz.diff", "outloc", "num_sets", 
                  "db_name", "num_nodes", "num_sets", "adduct_weights", 
                  "filter.by", "max.rt.diff", "max_isp", "MplusH.abundance.ratio.check", 
                  "mass_defect_window", "mass_defect_mode", 
                  "allsteps", "chemids_split"), file = "step1_results.Rda")
    
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
    
    rm(list = ls())
    load("tempobjects.Rda")
    
    if (allsteps == TRUE) {
      
      print("Status 3: Calculating scores for individual chemicals/metabolites")
      
      
      
      if (num_sets > num_nodes) {
    
        cl <- makeSOCKcluster(num_nodes)
        
        clusterEvalQ(cl, "multilevelannotationstep2")
        
        clusterExport(cl, "multilevelannotationstep2")
        
        clusterEvalQ(cl, "library(Rdisop)")
        clusterEvalQ(cl, "library(plyr)")
        # clusterEvalQ(cl, 'library(pryr)') clusterEvalQ(cl,
        # 'library(profmem)') clusterEvalQ(cl, 'library(gdata)')
        clusterExport(cl, "get_chemscorev1.6.71")
        clusterExport(cl, "getMolecule")
        clusterExport(cl, "ldply")
        # clusterExport(cl, 'mem_used')
        clusterExport(cl, "get_confidence_stage2")
        
        # clusterExport(cl, 'll')
        clusterExport(cl, "check_element")
        clusterExport(cl, "group_by_rt_histv2")
        clusterExport(cl, "adduct_table")
        clusterExport(cl, "adduct_weights")
        clusterExport(cl, "outloc")
        
        parLapply(cl, 1:num_sets, function(arg1) {
          
          cur_fname <- paste(outloc, "/stage2/chem_score", 
                             arg1, ".Rda", sep = "")
          check_if_exists <- suppressWarnings(try(load(cur_fname)))
          
          if (is(check_if_exists, "try-error")) {
            
            
            # suppressWarnings(multilevelannotationstep2(outloc1=outloc,list_number=arg1))
            multilevelannotationstep2(outloc1 = outloc, 
                                      list_number = arg1)
          } else {
            
            print(paste("List ", arg1, " already exists.", 
                        sep = ""))
          }
          
        })
        
        
        outloc1 = outloc    
        setwd(outloc1)
        list_number = 1
        
        load("step1_results.Rda")
        load("global_cor.Rda")
        unlink("allmatches_with_isotopes.txt")
        
        if (is.na(max.rt.diff) == FALSE) {
          
          max_diff_rt <- max.rt.diff
        }
        
        if (list_number > length(chemids_split)) {
          list_number <- length(chemids_split)
          # stop(paste('Invalid list number. Must be less than or
          # equal to ',length(chemids_split),sep=''))
          return(0)
        }
        
        if (list_number > num_sets) {
          list_number <- length(chemids_split)
          # stop(paste('Invalid list number. Must be less than or
          # equal to ',num_sets,sep=''))
          return(0)
        }
        
        outloc1 <- paste(outloc, "/stage2/", sep = "")
        suppressWarnings(dir.create(outloc1))
        setwd(outloc1)
        
        if (is.na(adduct_weights) == TRUE) {
          data(adduct_weights)
          
          adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
          adduct_weights1[1, ] <- c("M+H", 1)
          adduct_weights1[2, ] <- c("M-H", 1)
          adduct_weights <- as.data.frame(adduct_weights1)
          colnames(adduct_weights) <- c("Adduct", "Weight")
        }
        
        cnames <- colnames(mchemdata)
        cnames[2] <- "time"
        colnames(mchemdata) <- as.character(cnames)
        
        mchemdata$mz <- as.numeric(as.character(mchemdata$mz))
        mchemdata$time <- as.numeric(as.character(mchemdata$time))
        
        library(enviPat)
        data(isotopes)
        
        #for each set of formulae (as defined by chemids_split and list_number) get the list of associated 
        #    features based on mass defect and retention time
        chem_score <- lapply(chemids_split, function(j) {
          
          print(j)
          
          chemid <- chemids[j[1]]
          chemscoremat <- {
          }
          curmchemdata <- mchemdata[which(mchemdata$chemical_ID == 
                                            chemid), ]
          curmchemdata$mz <- as.numeric(as.character(curmchemdata$mz))
          curmchemdata$time <- as.numeric(as.character(curmchemdata$time))
          curmchemdata <- as.data.frame(curmchemdata)
          curmchemdata$Module_RTclust <- gsub(curmchemdata$Module_RTclust, 
                                              pattern = "_[0-9]*", replacement = "")
          isop_res_md$Module_RTclust <- gsub(isop_res_md$Module_RTclust, 
                                             pattern = "_[0-9]*", replacement = "")
          isp_masses_mz_data <- {
          }
          isp_masses_mz_data <- lapply(1:length(curmchemdata$mz),function(m) {
            isp_group <- as.character(curmchemdata$ISgroup[m])
            module_rt_group <- as.character(curmchemdata$Module_RTclust[m])
            module_rt_group <- gsub(module_rt_group,
                                    pattern = "_[0-9]*", replacement = "")
            query_md <- curmchemdata$mz[m] - round(curmchemdata$mz[m])
            put_isp_masses_curmz_data <- isop_res_md[which(abs((isop_res_md$MD) - (query_md)) < mass_defect_window & isop_res_md$Module_RTclust == module_rt_group), ]
                                         put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
            return(put_isp_masses_curmz_data)})
          
          isp_masses_mz_data <- ldply(isp_masses_mz_data, rbind)
          cnames <- colnames(isp_masses_mz_data)
          cnames[5] <- "AvgIntensity"
          colnames(isp_masses_mz_data) <- cnames
          isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)
          
          if (is.na(mass_defect_mode) == TRUE) {
            mass_defect_mode = "pos"
          }
          
          isp_masses_mz_data$mz <- as.numeric(as.character(isp_masses_mz_data$mz))
          isp_masses_mz_data$time <- as.numeric(as.character(isp_masses_mz_data$time))
          isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)
            
          out = get_chemscorev1.6.71_custom(chemicalid = chemid,
                                                   mchemicaldata = curmchemdata,
                                                   corthresh = corthresh,
                                                   global_cor = global_cor,
                                                   mzid = mzid,
                                                   max_diff_rt = max_diff_rt,
                                                   level_module_isop_annot = isp_masses_mz_data,
                                                   adduct_table = adduct_table,
                                                   adduct_weights = adduct_table$Weights,
                                                   filter.by=c("M+H"),
                                                   max_isp = max_isp,
                                                   MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
                                                   mass_defect_window = mass_defect_window,
                                                   mass_defect_mode = mass_defect_mode,
                                                   outlocorig = outloc,
                                                   iso_ppm_tol = 5,
                                                   iso_int_tol = 0.05)
          
          return(out)
          
          
          
        })
        
        # stopCluster(cl)
        cur_fname <- paste("chem_score", list_number, "_a.Rda", 
                           sep = "")
        # save(chem_score,file=cur_fname)
        
        # mc.cores=num_nodes,mc.preschedule=FALSE)
        # stopCluster(cl)
        
        chem_score2 <- chem_score[which(chem_score != "NULL")]
        
        rm(chem_score)
        
        curchemscoremat <- ldply(chem_score2, rbind)
        rm(chem_score2)
        
        # save(chem_score2,file='chem_scoreA.Rda')
        cur_fname <- paste("chem_score", list_number, ".Rda", 
                           sep = "")
        save(curchemscoremat, file = cur_fname)
        # chemscoremat<-rbind(chemscoremat,curchemscoremat)
        
        Sys.sleep(1)
        
        rm("curchemscoremat", "mchemdata", "chemids", "adduct_table", 
           "global_cor", "mzid", "max_diff_rt", "isop_res_md", 
           "corthresh", "level_module_isop_annot", "chemids_split", 
           "corthresh", "max.mz.diff", "outloc", "num_sets", 
           "db_name", "num_nodes", "num_sets", "adduct_weights", 
           "filter.by")
        
        
        rm(list = ls())
        
          
  }
        
        if (length(chem_score) > 0) {
          if (chem_score$chemical_score >= (-100)) {
            chem_score$filtdata <- chem_score$filtdata[order(chem_score$filtdata$mz),]
          
          if (nrow(chem_score$filtdata) > 0) {
            cur_chem_score <- chem_score$chemical_score
            chemscoremat <- cbind(cur_chem_score, chem_score$filtdata)
          } else {
            # print(chem_score)
          }
          
          rm(chem_score)
          chemscoremat <- na.omit(chemscoremat)
          chemscoremat <- as.data.frame(chemscoremat)
          }
        } else {
          rm(chem_score)
        }
        
        rm("curmchemdata", "isp_masses_mz_data", "mzid_cur", 
           "chemid")
      
        suppressWarnings(rm(hmdbCompMZ))
        suppressWarnings(rm(hmdbAllinf))
          
        # chemscoremat<-chemscoremat[,c(1:7,9,11,14)]
        # chemscoremat$MatchCategory[which(chemscoremat$MatchCategory=='Multiple')]<-'M'
        # chemscoremat$MatchCategory[which(chemscoremat$MatchCategory=='Unique')]<-'U'
        # chemscoremat<-as.matrix(chemscoremat)
        
        return(chemscoremat)
        #,mc.cores=num_nodes,mc.preschedule=FALSE)
      }
    }
