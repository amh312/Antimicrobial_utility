##PRESCRIPTION MODEL HYPERPARAMETER TUNING

set.seed(123)

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Binding missing demographic information
binddems <- function(df) {
  
  hadm_key <- hadm %>% select(subject_id,race,language,marital_status) %>% 
    distinct(subject_id,.keep_all=T)
  age_key <- pats %>% select(subject_id,anchor_age) %>% distinct(subject_id,.keep_all = T)
  dems_key <- left_join(hadm_key,age_key,by="subject_id")
  
  assign("age_key",age_key,.GlobalEnv)
  
  df %>% left_join(dems_key,by="subject_id") %>% mutate(anchor_age=as.numeric(anchor_age))
  
}

###Preprocessing of model and microsimulation dataframes
abx_ur_datpreproc <- function(ab_df,ab_predname,ab_outcname,ab_combname,
                              abur_df,abur_predname,abur_outcname,abur_combname,
                              probs_df) {
  
  set.seed(123)
  
  #convert prescription df outcomes to 1s and 0s
  ab_outcomes <- ab_df %>%
    select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                       overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
  
  #select predictor variables of choice
  ab_predictors <- ab_df %>% select(pHADM:age65,prAKI:pDIAB,pCARD:curr_service,pICU:pSEPSIS,ob_freq,highCRP,
                                    pc_dyspnea:SIRS,
                                    ab_name_Ampicillin_Ceftriaxone:ab_name_Ampicillin,race:language,anchor_age)
  
  #make boolean dummy features from predictors
  ab_dummy <- dummyVars(" ~ .", data = ab_predictors)
  ab_predictors <- predict(ab_dummy, newdata = ab_predictors)
  
  #put outcomes and predictors back together
  ab_combined <- as.data.frame(cbind(ab_outcomes, ab_predictors))
  
  #convert microsim df outcomes to 1s and 0s
  ur_outcomes <- abur_df %>%
    select(micro_specimen_id,CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                                         overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
  
  #join microsim outcomes to antimicrobials in probs df
  ur_outcomes <- probs_df %>%
    left_join(ur_outcomes, relationship = 'many-to-one') %>% 
    select(-c(id_no,I:subject_id))
  
  #filter microsim to predictor columns based on those in prescriptions df
  common_columns <- intersect(names(ab_predictors),names(abur_df))
  ur_predictors <- abur_df %>% select(micro_specimen_id,all_of(common_columns))
  
  #join predictor vars to antimicrobials in probs df
  ur_predictors <- probs_df %>% left_join(ur_predictors,
                                          relationship = 'many-to-one') %>% 
    select(-c(id_no,I:subject_id)) %>% rename(ab_name="Antimicrobial")
  
  #make dummy features from microsim predictors
  dummy_ur <- dummyVars(" ~ .", data = ur_predictors)
  ur_predictors <- predict(dummy_ur, newdata = ur_predictors)
  
  #put microsim predictors and outcomes back together
  ur_combined <- as.data.frame(cbind(ur_outcomes, ur_predictors))
  
  #standardise naming of ab name feature variables
  colnames(ur_combined)[grepl("ab_name",colnames(ur_combined))] <- colnames(ur_combined)[grepl("ab_name",colnames(ur_combined))] %>% 
    str_replace("ab_name","ab_name_") %>% 
    str_replace_all("-",".")
  
  assign(ab_predname,ab_predictors,envir = .GlobalEnv)
  assign(ab_outcname,ab_outcomes,envir = .GlobalEnv)
  assign(ab_combname,ab_combined,envir = .GlobalEnv)
  assign(abur_predname,ur_predictors,envir = .GlobalEnv)
  assign(abur_outcname,ur_outcomes,envir = .GlobalEnv)
  assign(abur_combname,ur_combined,envir = .GlobalEnv)
  
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
TTsplitter <- function(dataset,outc,trainprop,chosnametrain,chosnametest){
  
  #get partition index based on specified proportion
  trainindex <- createDataPartition(dataset[[outc]], p = trainprop, list = FALSE, times = 1)
  
  #index training dataset
  urdftrain <- dataset[trainindex, ]
  
  #index testing dataset
  urdftest <- dataset[-trainindex, ]
  
  #assign to objects with 'Train' and 'Test' suffixes replacing '_combined' suffix
  assign(chosnametrain,
         urdftrain,.GlobalEnv)
  assign(chosnametest,
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
                           roundstorun,eta_tune=FALSE,earlystop=50){
  
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
    params = hypparamlist,
    data = traindat,
    nrounds = roundstorun,
    nfold = 5,
    early_stopping_rounds = earlystop,
    verbose = 1,
  )
  
  #get boosting round number of best auc
  best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
  
  #get best auc
  best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
  
  #if auc better than previous round, update values
  if ((best_iteration_auc > best_auc & eta_tune==FALSE)|
      eta_tune==TRUE){
    
    best_auc <- best_iteration_auc
    best_params <- hypparamlist
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

##Read-in

probs_df_overall <- read_csv("probs_df_overall_postur.csv")
pats <- read_csv("patients.csv")
ur_xg <- read_csv("interim_ur_util.csv")
hadm <- read_csv("admissions.csv")
train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")

##Final preprocessing and dataset train/test splits

###Make full abx dataframe, factorise outcomes, and add missing demographics
abx <- data.frame(rbind(train_abx,test_abx)) %>% factorise() %>% binddems()
ur_xg <- ur_xg %>% left_join(age_key,by="subject_id")

###Preprocessing outcomes and predictors
abx_ur_datpreproc(abx,"abx_predictors","abx_outcomes","abx_combined",
                  ur_xg,"ur_abx_predictors","ur_abx_outcomes","ur_abx_combined",
                  probs_df_overall)

#Write CSVs for later training/validation
write_csv(abx_combined,"abx_combined.csv")
write_csv(ur_abx_combined,"ur_abx_combined.csv")

###Latin hypercube hyperparameter grid for max depth and min child weight
depchild_lhsgrid <- paramgrid_lhs(10,"max_depth",2,9,"min_child_weight",1,10,0)

###Set first empty hyperparameter list
cdi_tox_max_child_bestparams <- c()

###Run first hyperparameter tune (max depth and min child weight)
for (outcome in colnames(abx_outcomes)) {
  
  ###Set starting point for AUROC
  best_auc <- 0
  
  #ensure no na values
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #split df to train and test
    abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
    
    #make xgboost training matrix
    abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
    
    #iterate over lhs-selected hyperparameter pairs
    for (i in 1:nrow(depchild_lhsgrid)) {
      
      #run cross-validation and extract best hyperparameters
      cv_model <- abxtrain_matrix %>%
        hypparam_tuner(et=0.05,md=1,mcw=1,ss=0.8,csb=0.8,
                       "max_depth","min_child_weight",depchild_lhsgrid,
                       i,outcome,"best_auc","best_params","best_nrounds",TRUE,
                       50)
      
    }
    
    #add best parameters to list
    cdi_tox_max_child_bestparams[[outcome]] <- best_params
    
  }
}

cdi_tox_max_child_bestparams[['overall_tox']] %>% view()

###Save first set of tuned hyperparameters to csvs
save_hypparams(cdi_tox_max_child_bestparams,"max_child_",
               names(abx_outcomes))

###Read first hyperparameter set back in to list to check save
cdi_tox_max_child_bestparams <- hypparamreader("max_child_",names(abx_outcomes))

##Hyperparameter tuning 2/3 (subsample and col sample by tree)

###Latin hypercube hyperparameter grid for subsample and col sample by tree
colsub_lhsgrid <- paramgrid_lhs(10,"subsample",0.5,1,"colsample_bytree",0.5,1,2)

###Set second empty hyperparameter list
cdi_tox_col_sub_bestparams <- c()

###Run second hyperparameter tune (subsample and colsample by tree)
for (outcome in colnames(abx_outcomes)) {
  
  #set starting point for best auc
  best_auc <- 0
  
  #check for nas
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #train-test split
    abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
    
    #xgboost training matrix
    abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
    
    #iterate across hyperparameter values
    for (i in 1:nrow(colsub_lhsgrid)) {
      
      #run cross_validation and extract best hyperparameter values
      cv_model <- abxtrain_matrix %>%
        hypparam_tuner(et=0.05,md=cdi_tox_max_child_bestparams[[outcome]]$max_depth,
                       mcw=cdi_tox_max_child_bestparams[[outcome]]$min_child_weight,
                       ss=0.8,csb=0.8,"subsample","colsample_bytree",
                       colsub_lhsgrid,i,outcome,
                       "best_auc","best_params","best_nrounds",TRUE)
      
    }
    
    #add hyperparameters to model list
    cdi_tox_col_sub_bestparams[[outcome]] <- best_params
    
  }
}

###Save updated set of tuned hyperparameters to csvs
save_hypparams(cdi_tox_col_sub_bestparams,"cdi_tox_col_sub_",
               colnames(abx_outcomes))

###Read second hyperparameter set back in to list to check save
cdi_tox_col_sub_bestparams <- hypparamreader("cdi_tox_col_sub_",colnames(abx_outcomes))

##Hyperparameter tuning 3/3 (subsample and col sample by tree)

###Set final empty hyperparameter list
cdi_tox_final_bestparams <- c()

###Final hyperparameter tuning (learning rate)
for (outcome in colnames(abx_outcomes)) {
  
  #set starting point for auc
  best_auc <- 0
  
  #check for nas
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    #set seed
    set.seed(123)
    
    #train-test split
    abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
    
    #xgboost training matrix
    abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
    
    #starting point for n rounds run
    n_rounds_run <- 0
    
    #starting opoint for learning rate
    parameter_val <- 0.1
    
    #starting point for iteration number
    i <- 1
    
    #set max tree depth outside of while loop so it can be reassigned if required
    md_val <- cdi_tox_max_child_bestparams[[outcome]]$max_depth
    
    #set round limit
    roundlimit <- 1000
    
    stoprounds <- 50
    
    #run until nrounds is between 300 and 1,000
    while(n_rounds_run<300|n_rounds_run==roundlimit) {
      
      #run cross-validation and extract best values
      abxtrain_matrix %>%
        hypparam_tuner(et=parameter_val,md=md_val,
                       mcw=cdi_tox_max_child_bestparams[[outcome]]$min_child_weight,
                       ss=cdi_tox_col_sub_bestparams[[outcome]]$subsample,
                       csb=cdi_tox_col_sub_bestparams[[outcome]]$colsample_bytree,
                       hypparam_1="",hypparam_2="",lhsgrid=NULL,
                       i,outcome,"best_auc","best_params","best_nrounds",FALSE,
                       roundlimit,TRUE,earlystop = stoprounds)
      
      #if <300 rounds, halve learning rate
      if(n_rounds_run<300) {
        
        parameter_val <- parameter_val/2
        
        #if > 1,000 rounds, add 0.1 to learning rate
      } else if (n_rounds_run==roundlimit) {
        
        parameter_val <- parameter_val+0.1
        
      }
      
      #next iteration
      i <- i+1
      
      #if learning rate goes beyond 0.3, fix max tree depth at 6, increase max rounds and reduce stopping rounds
      if(parameter_val==0.4 & md_val>6) {
        
        md_val <- 6
        
        stoprounds <- 10
        
        roundlimit <- 2000
        
        parameter_val <- 0.3
        
        #if still not converging, lower tree depth further
      } else if (parameter_val==0.4 & md_val>=4 & md_val<=6) {
        
        md_val <- md_val-1
        
        parameter_val <- 0.3
        
        #if that doesn't work when the tree depth hits 3, abort
      } else if (parameter_val==0.4 & md_val==3) {
        
        stop("Model not converging. Review hyperparameters")
        
      }
      
    }
    
    #add best hyperparameters to model list
    cdi_tox_final_bestparams[[outcome]] <- best_params
    cdi_tox_final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
  }
  
  
}

###Save final set of tuned hyperparameters to csvs
save_hypparams(cdi_tox_final_bestparams,"cdi_tox_final_params_",
               colnames(abx_outcomes))

###Read final hyperparameter set back in to list to check save
cdi_tox_final_bestparams <- hypparamreader("cdi_tox_final_params_",colnames(abx_outcomes))
