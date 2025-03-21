###COVERAGE MODEL HYPERPARAMETER TUNING - STAGE 1

set.seed(123)

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Generate combinations of antibiotics in a list
ab_combiner <- function(abx_list,newname,combname) {
  
  #get list of all combinations
  comblist <- combn(abx_list, 2, function(x) paste(x, collapse = "_"))
  
  #combine singles with combinations
  all_list <- c(abx_list,comblist)
  
  #assign to global environment
  assign(newname,comblist,envir = .GlobalEnv)
  assign(combname,all_list,envir = .GlobalEnv)
  
}

###Generate combinations of antibiotics in a mapped list
abmap_combiner <- function(ab_map) {
  
  #get all 2-ab combinations in the map
  abcombos_map <- combn(names(ab_map), 2, simplify = FALSE)
  
  #bind combos on to singles
  c(ab_map,
    setNames(
      
      #bind on combo elements themselves
      lapply(abcombos_map, function(x) paste(ab_map[x], collapse = "_")),
      
      #bind on combo element names
      sapply(abcombos_map, function(x) paste(x, collapse = "_"))
    )
  )
  
}

###Make filtered, full and short antimicrobial maps
abcombo_variants <- function(abmap,fullname,mainname,shortname,abreflist) {
  
  #get vector of all antibiotics as unlabelled list
  fullmap <- abmap %>% unlist()
  names(fullmap) <- NULL
  
  #get combined ab map just for abs in prescription df
  abmap2 <- abmap[names(abmap) %in% abreflist]
  
  #get unlabelled vector of the shortened combo list
  shortmap <- abmap2 %>% unlist()
  names(shortmap) <- NULL
  
  #assign new objects to glob env
  assign(fullname,fullmap,envir = .GlobalEnv)
  assign(mainname,abmap2,envir = .GlobalEnv)
  assign(shortname,shortmap,envir = .GlobalEnv)
  
}

###Replace antimicrobial short name with long name including for combinations
abcombo_replace <- function(abtarget, map) {
  
  #switch names of combined ab list with the list elements themselves
  flip_abmap <- setNames(names(map), map)
  
  abtarget %>%
    
    #convert column to character format
    as.character() %>%
    
    #if target value is in names of flipped map, return the element
    sapply(function(x) if (x %in% names(flip_abmap)) flip_abmap[[x]] else x)
  
}

###Replace antimicrobial short name with long name including for combinations
abcombo_reverse <- function(abtarget, map) {
  
  abtarget %>%
    
    #convert to character
    as.character() %>%
    
    #if target ab is in the map elements, return that name
    sapply(function(x) if (x %in% map) names(map)[map == x] else x)
  
}

###Preprocessing of model and microsimulation dataframes
ur_datpreproc <- function(ur_df,ur_outcomename,ur_predictorname,ur_comboname,
                          microsim_df,microsim_outcomename,microsim_predictorname,
                          microsim_comboname) {
  
  mutmar <- function(dfx) {
    
    dfx %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                            TRUE~marital_status))
    
  }
  binariseast <- function(dfx) {
    
    dfx %>% select(all_of(shortmap)) %>%
      mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                     ifelse(. == "S" | . == "I", 1, NA))))
    
  }
  dummyer <- function(dfx) {
    
    urfeatdummies <- dummyVars(" ~ .", data = dfx)
    predict(urfeatdummies, newdata = dfx)
    
  }
  
  
  ur_df_outcomes <- ur_df %>% mutmar() %>% binariseast()
  microsim_df_outcomes <- microsim_df %>% mutmar() %>% binariseast()
  
  ur_df_predictors <- ur_df %>% select(!all_of(fullmap)) %>% dummyer()
  microsim_df_predictors <- microsim_df %>%
    select(any_of(colnames(ur_df_predictors))) %>% dummyer()
  
  ur_df_combined <- as.data.frame(cbind(ur_df_outcomes, ur_df_predictors))
  microsim_df_combined <- as.data.frame(cbind(microsim_df_outcomes, microsim_df_predictors))
  
  assign(ur_outcomename,ur_df_outcomes,.GlobalEnv)
  assign(ur_predictorname,ur_df_predictors,.GlobalEnv)
  assign(ur_comboname,ur_df_combined,.GlobalEnv)
  assign(microsim_outcomename,microsim_df_outcomes,.GlobalEnv)
  assign(microsim_predictorname,microsim_df_predictors,.GlobalEnv)
  assign(microsim_comboname,microsim_df_combined,.GlobalEnv)
  
  
}

###Generate latin hypercube parameter grid for 2 hyperparameters
paramgrid_lhs <- function(n_samples,hypname_1,ranglow_1,rangup_1,
                          hypname_2,ranglow_2,rangup_2) {
  
  rang_1 <- c(ranglow_1,rangup_1)
  rang_2 <- c(ranglow_2, rangup_2)
  lhs_sample <- randomLHS(n_samples, 2)
  hypset1 <- round(lhs_sample[, 1] * (rang_1[2] - rang_1[1]) + rang_1[1])
  hypset2 <- round(lhs_sample[, 2] * (rang_2[2] - rang_2[1]) + rang_2[1])
  hyppairs <- data.frame(col1 = hypset1, col2 = hypset2)
  colnames(hyppairs) <- c(hypname_1,hypname_2)
  hyppairs
  
}

###Generate XGBoost matrix from model train or test dataset
model_matrixmaker <- function(dataset,predictors,outc) {
  
  predictorcols <- colnames(predictors)
  selected_columns <- intersect(predictorcols, colnames(dataset))
  
  xgb.DMatrix(data = as.matrix(dataset %>% select(all_of(selected_columns))), 
              label = dataset[[outc]])
  
  
}

###Final preprocessing

##Read-in

train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
hadm <- read_csv("admissions.csv")
pats <- read_csv("patients.csv")

##Model tuning

###Antimicrobial mapping lists
antimicrobial_map <- c("Ampicillin" = "AMP","Ampicillin-sulbactam" = "SAM",
                       "Piperacillin-tazobactam" = "TZP","Cefazolin" = "CZO","Ceftriaxone" = "CRO",
                       "Ceftazidime" = "CAZ","Cefepime" = "FEP","Meropenem" = "MEM","Ciprofloxacin" = "CIP",
                       "Gentamicin" = "GEN","Trimethoprim-sulfamethoxazole" = "SXT","Nitrofurantoin" = "NIT",
                       "Vancomycin" = "VAN")
all_singles <- antimicrobial_map %>% unname()
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM","GEN","SXT")
oral_singles <- c("AMP","SAM","CIP","SXT","NIT")
access_singles <- c("AMP","SAM","GEN","SXT","NIT","CZO")
watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
singlelists <- c("all_singles","iv_singles","oral_singles","access_singles","watch_singles")
comblists <- str_replace(singlelists,"singles","combos")
alllists <- c("all_abs","all_ivs","all_orals","all_access","all_watch")
for (i in seq_along(singlelists)) {ab_combiner(get(singlelists[i]),comblists[i],alllists[i])}
combined_antimicrobial_map <- abmap_combiner(antimicrobial_map)
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
abcombo_variants(combined_antimicrobial_map,
                 "fullmap","combined_antimicrobial_map","shortmap",abx_in_train)

###Final preprocessing and dataset train/test splits
ur_datpreproc(urines5,"urines5_outcomes","urines5_predictors","urines5_combined",
              ur_xg,"ur_xg_outcomes","ur_xg_predictors","ur_xg_combined")

###Latin hypercube hyperparameter grid for max depth and min child weight
lhs_hypgrid <- paramgrid_lhs(10,"max_depth",2,9,"min_child_weight",1,10)


max_child_bestparams <- c()

for (outcome in colnames(urines5_outcomes)) {
  
  best_auc <- 0
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    set.seed(123)
    
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = 0.8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
    
    for (i in 1:nrow(lhs_hypgrid)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = lhs_hypgrid %>% select(max_depth) %>% dplyr::slice(i) %>% unlist(),
        min_child_weight = lhs_hypgrid %>% select(min_child_weight) %>% dplyr::slice(i) %>% unlist(),
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = urtrain_matrix,
        nrounds = 50,
        nfold = 5,
        early_stopping_rounds = 50,
        verbose = 1,
      )
      
      best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
      best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
      cv_model$evaluation_log$test_logloss_mean
      if (best_iteration_auc > best_auc) {
        best_auc <- best_iteration_auc
        best_params <- params
        best_nrounds <- best_iteration_index
      }
      
    }
    
    max_child_bestparams[[outcome]] <- best_params
    
  }
}
for (i in 1:length(max_child_bestparams)) {
  
  maxy <- data.frame(max_child_bestparams[i])
  
  write_csv(maxy,glue("max_child_{combined_antimicrobial_map[i]}.csv"))
  
}
for (outcome in colnames(urines5_outcomes)){
  print(max_child_bestparams[[outcome]]$min_child_weight)
}



