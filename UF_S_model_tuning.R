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
  
  #sub-function for filling in marital status NAs
  mutmar <- function(dfx) {
    
    dfx %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                            TRUE~marital_status))
    
  }
  
  #sub-function to select and binarising ast results to 1s and 0s
  binariseast <- function(dfx) {
    
    dfx %>% select(all_of(shortmap)) %>%
      mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                     ifelse(. == "S" | . == "I", 1, NA))))
    
  }
  
  #sub-function to make dummy variables
  dummyer <- function(dfx) {
    
    urfeatdummies <- dummyVars(" ~ .", data = dfx)
    predict(urfeatdummies, newdata = dfx)
    
  }
  
  #fill marital status and select/binarise ast results
  ur_df_outcomes <- ur_df %>% mutmar() %>% binariseast()
  microsim_df_outcomes <- microsim_df %>% mutmar() %>% binariseast()
  
  #make predictor dataframes
  ur_df_predictors <- ur_df %>% select(!all_of(fullmap)) %>% dummyer()
  microsim_df_predictors <- microsim_df %>%
    select(any_of(colnames(ur_df_predictors))) %>% dummyer()
  
  #make combined dataframes
  ur_df_combined <- as.data.frame(cbind(ur_df_outcomes, ur_df_predictors))
  microsim_df_combined <- as.data.frame(cbind(microsim_df_outcomes, microsim_df_predictors))
  
  #assign resulting dfs to global environment
  assign(ur_outcomename,ur_df_outcomes,.GlobalEnv)
  assign(ur_predictorname,ur_df_predictors,.GlobalEnv)
  assign(ur_comboname,ur_df_combined,.GlobalEnv)
  assign(microsim_outcomename,microsim_df_outcomes,.GlobalEnv)
  assign(microsim_predictorname,microsim_df_predictors,.GlobalEnv)
  assign(microsim_comboname,microsim_df_combined,.GlobalEnv)
  
  
}

###Generate latin hypercube parameter grid for 2 hyperparameters
paramgrid_lhs <- function(n_samples,hypname_1,ranglow_1,rangup_1,
                          hypname_2,ranglow_2,rangup_2,roundval) {
  
  #get latin hypercube samples between zero and 1
  lhs_sample <- randomLHS(n_samples, 2)
  
  #multiply samples by range then add minimum to get into range
  hypset1 <- round(lhs_sample[, 1] * (rangup_1 - ranglow_1) + ranglow_1,roundval)
  hypset2 <- round(lhs_sample[, 2] * (rangup_2 - ranglow_2) + ranglow_2,roundval)
  
  #bind two sets of hyperparameters together
  hyppairs <- data.frame(col1 = hypset1, col2 = hypset2)
  
  #set names for df
  colnames(hyppairs) <- c(hypname_1,hypname_2)
  hyppairs
  
}

###Train-test splitting
TTsplitter <- function(dataset,outc,trainprop){
  
  #get partition index based on specified proportion
  trainindex <- createDataPartition(dataset[[outcome]], p = trainprop, list = FALSE, times = 1)
  
  #index training dataset
  urdftrain <- dataset[trainindex, ]
  
  #index testing dataset
  urdftest <- dataset[-trainindex, ]
  
  #assign to objects with 'Train' and 'Test' suffixes replacing '_combined' suffix
  assign(glue("{deparse(substitute(urines5_combined)) %>% str_remove('_combined')}Train"),
         urdftrain,.GlobalEnv)
  assign(glue("{deparse(substitute(urines5_combined)) %>% str_remove('_combined')}Test"),
         urdftest,.GlobalEnv)
  
}

###Generate XGBoost matrix from model train or test dataset
model_matrixmaker <- function(dataset,predictors,outc) {
  
  #get predictor feature names
  predictorcols <- colnames(predictors)
  
  #ensure that all predictor column names are in specified data
  selected_columns <- intersect(predictorcols, colnames(dataset))
  
  #make xgboost matrix
  xgb.DMatrix(data = as.matrix(dataset %>% select(all_of(selected_columns))), 
              label = dataset[[outc]])
  
  
}

###Hyperparameter pair tuner
hypparam_tuner <- function(traindat,et,md,mcw,ss,csb,hypparam_1,hypparam_2,
                           lhsgrid,itera,outc,aucname,paramname,nroundname,lhcpairs,
                           roundstorun){
  
  #update message
  print(glue("Running CV {itera} for {outc}..."))
  
  #make default hyperparameter list
  hypparamlist <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = et,
    max_depth = md,
    min_child_weight = mcw,
    subsample = ss,
    colsample_bytree = csb
  )
  
  if (lhcpairs==TRUE) {
  
    #replace specified hyperparameters with latin hypercube selections for row i
    hypparamlist[[hypparam_1]] <- lhsgrid %>% select(!!sym(hypparam_1)) %>% dplyr::slice(itera) %>% unlist()
    hypparamlist[[hypparam_2]] <- lhsgrid %>% select(!!sym(hypparam_2)) %>% dplyr::slice(itera) %>% unlist()
  
  }
  
  #run 5-fold cross-validation with 50 rounds
  cv_model <- xgb.cv(
    params = params,
    data = traindat,
    nrounds = roundstorun,
    nfold = 5,
    early_stopping_rounds = 50,
    verbose = 1,
  )
  
  #get boosting round number of best auc
  best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
  
  #get best auc
  best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
  
  #if auc better than previous round, update values
  if (best_iteration_auc > best_auc) {
    best_auc <- best_iteration_auc
    best_params <- params
    best_nrounds <- best_iteration_index
  }
  
  #get number of rounds run (for learning rate tuning)
  n_roundrn <- cv_model$evaluation_log %>% nrow()
  
  #write values to global environment under specified names
  assign(aucname,best_auc,.GlobalEnv)
  assign(paramname,best_params,.GlobalEnv)
  assign(nroundname,best_nrounds,.GlobalEnv)
  assign("n_rounds_run",n_roundrn,.GlobalEnv)
  
  
}

###Save hyperparameter list to csvs
save_hypparams <- function(hyplist,namestring,namesmap){
  
  #iterate over list of models
  for (i in 1:length(hyplist)) {
    
    #make dataframe from parameters
    maxy <- data.frame(hyplist[i])
    
    #write to csv
    write_csv(maxy,glue("{namestring}{namesmap[i]}.csv"))
    
  }
  
}

###Reading in hyperparameters back to list from csvs
hypparamreader <- function(namestring,namesmap) {
  
  #write file nanmes
  file_names <- glue("{namestring}{namesmap}.csv")
  
  #apply read_csv across name list
  hypparlist <- lapply(file_names, function(file) {
    read_csv(file)
  })
  
  #add names
  namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
                "best_nrounds")
  
  #iteratively add names to sublists
  for (i in 1:length(hypparlist)) {
    names(hypparlist[[i]])[3:8] <- namelist
  }
  
  #add model names to list top level
  names(hypparlist) <- namesmap
  
  hypparlist
  
}

###Final preprocessing

##Read-in

train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
hadm <- read_csv("admissions.csv")
pats <- read_csv("patients.csv")

##Antimicrobial mapping lists

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

##Final preprocessing and dataset train/test splits

ur_datpreproc(urines5,"urines5_outcomes","urines5_predictors","urines5_combined",
              ur_xg,"ur_xg_outcomes","ur_xg_predictors","ur_xg_combined")

##Hyperparameter tuning 1/3 (max depth and min child weight)

###Latin hypercube hyperparameter grid for max depth and min child weight
depchild_lhsgrid <- paramgrid_lhs(10,"max_depth",2,9,"min_child_weight",1,10,0)

###Set first empty hyperparameter list
max_child_bestparams <- c()

###Run first hyperparameter tune (max depth and min child weight)
for (outcome in colnames(urines5_outcomes)) {
  
  #set starting point for auc
  best_auc <- 0
  
  #check there are no na values
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #split df to train and test
    urines5_combined %>% TTsplitter(outcome,0.8)
    
    #make xgboost training matrix
    urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
    
    #iterate over latin hypercube grid
    for (i in 1:nrow(depchild_lhsgrid)) {
      
      #run cross-validation and extract best hyperparameters
      cv_model <- urtrain_matrix %>%
        hypparam_tuner(et=0.05,md=1,mcw=1,ss=0.8,csb=0.8,
                       "max_depth","min_child_weight",depchild_lhsgrid,
                       i,outcome,"best_auc","best_params","best_nrounds",TRUE,
                       50)
      
    }
    
    #add hyperparameters to model list
    max_child_bestparams[[outcome]] <- best_params
    
  }
}

###Save first set of tuned hyperparameters to csvs
save_hypparams(max_child_bestparams,"max_child_",
               combined_antimicrobial_map)

###Read first hyperparameter set back in to list to check save
max_child_bestparams <- hypparamreader("max_child_",combined_antimicrobial_map)

##Hyperparameter tuning 2/3 (subsample and col sample by tree)

###Latin hypercube hyperparameter grid for subsample and col sample by tree
colsub_lhsgrid <- paramgrid_lhs(10,"subsample",0.5,1,"colsample_bytree",0.5,1,2)

###Set second empty hyperparameter list
col_sub_bestparams <- c()

###Run second hyperparameter tune (subsample and colsample by tree)
for (outcome in colnames(urines5_outcomes)) {
  
  #set starting point for best auc
  best_auc <- 0
  
  #check for nas
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #train-test split
    urines5_combined %>% TTsplitter(outcome,0.8)
    
    #xgboost training matrix
    urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
    
    #iterate across hyperparameter values
    for (i in 1:nrow(colsub_lhsgrid)) {
      
      #run cross_validation and extract best hyperparameter values
      cv_model <- urtrain_matrix %>%
        hypparam_tuner(et=0.05,md=max_child_bestparams[[outcome]]$max_depth,
                       mcw=max_child_bestparams[[outcome]]$min_child_weight,
                       ss=0.8,csb=0.8,"subsample","colsample_bytree",
                       colsub_lhsgrid,i,outcome,
                       "best_auc","best_params","best_nrounds",TRUE,50)
      
    }
    
    #add hyperparameters to model list
    col_sub_bestparams[[outcome]] <- best_params
    
  }
}

###Save updated set of tuned hyperparameters to csvs
save_hypparams(col_sub_bestparams,"col_sub_",
               combined_antimicrobial_map)

###Read second hyperparameter set back in to list to check save
col_sub_bestparams <- hypparamreader("col_sub_",combined_antimicrobial_map)

##Hyperparameter tuning 3/3 (subsample and col sample by tree)

###Set final empty hyperparameter list
final_bestparams <- c()

###Final hyperparameter tuning (learning rate)
for (outcome in colnames(urines5_outcomes)) {
  
  #set starting point for auc
  best_auc <- 0
  
  #check for nas
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #train-test split
    urines5_combined %>% TTsplitter(outcome,0.8)
    
    #xgboost training matrix
    urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
    
    #starting point for n rounds run
    n_rounds_run <- 0
    
    #starting opoint for learning rate
    parameter_val <- 0.1
    
    #starting point for iteration number
    i <- 1
    
    #set max tree depth outside of while loop so it can be reassigned if required
    md_val <- max_child_bestparams[[outcome]]$max_depth
    
    #run until nrounds is between 200 and 1,000
    while(n_rounds_run<200|n_rounds_run==1000) {
      
      #run cross-validation and extract best values
      urtrain_matrix %>%
        hypparam_tuner(et=parameter_val,md=md_val,
                       mcw=max_child_bestparams[[outcome]]$min_child_weight,
                       ss=col_sub_bestparams[[outcome]]$subsample,
                       csb=col_sub_bestparams[[outcome]]$colsample_bytree,
                       hypparam_1="",hypparam_2="",lhsgrid=NULL,
                       i,outcome,"best_auc","best_params","best_nrounds",FALSE,
                       1000)
      
      #if <200 rounds, halve learning rate
      if(n_rounds_run<200) {
        
        parameter_val <- parameter_val/2
        
        #if > 1,000 rounds, add 0.1 to learning rate
      } else if (n_rounds_run==1000) {
        
        parameter_val <- parameter_val+0.1
        
      }
      
      #next iteration
      i <- i+1
      
      #if learning rate goes beyond 0.3, fix max tree depth at 0.6 and try again
      if(parameter_val==0.4 & md_val!=0.6) {
        
        md_val <- 0.6
        
        parameter_val <- 0.3
        
        #if that doesn't work, abort
      } else if (parameter_val==0.4 & md_val==0.6) {
        
        break
        
        print("Aborting. Review model parameters")
        
      }
    
      
    }
    
    #add best hyperparameters to model list
    final_bestparams[[outcome]] <- best_params
    final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
  }
  
  
}

###Save final set of tuned hyperparameters to csvs
save_hypparams(final_bestparams,"final_params_",
               combined_antimicrobial_map)

###Read final hyperparameter set back in to list to check save
final_bestparams <- hypparamreader("final_params_",combined_antimicrobial_map)

##Saving interim CSV
write_csv(urines5_combined,"urines5_combined.csv")
