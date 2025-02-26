##CDI AND TOXICITY HYPERPARAMETER TUNING - STAGE 3

set.seed(123)

###Third round of hyperparameter tuning (learning rate and best iteration)

###Read_in previous parameters
file_names <- glue("cdi_tox_col_sub_{names(abx_outcomes)}.csv")
cdi_tox_col_sub_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(cdi_tox_col_sub_bestparams)) {
  names(cdi_tox_col_sub_bestparams[[i]])[3:8] <- namelist
}
names(cdi_tox_col_sub_bestparams) <- names(abx_outcomes)

file_names <- glue("max_child_{names(abx_outcomes)}.csv")
cdi_tox_max_child_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(cdi_tox_max_child_bestparams)) {
  names(cdi_tox_max_child_bestparams[[i]])[3:8] <- namelist
}
names(cdi_tox_max_child_bestparams) <- names(abx_outcomes)

cdi_tox_final_bestparams <- c()

for (outcome in colnames(abx_outcomes)) {

  best_auc <- 0
    
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    set.seed(123)
    trainIndex <- createDataPartition(abx_combined[[outcome]], p = .8, list = FALSE, times = 1)
    abxTrain <- abx_combined[trainIndex, ]
    abxTest <- abx_combined[-trainIndex, ]
    
    predictor_columns <- colnames(abx_predictors)
    selected_columns <- intersect(predictor_columns, colnames(abxTrain))
    missing_cols <- setdiff(selected_columns, colnames(ur_abx_combined))
    ur_abx_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                                label = abxTrain[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                               label = abxTest[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(all_of(selected_columns))), 
                                label = ur_abx_combined[[outcome]])
    
    n_rounds_run <- 0
    parameter_val <- 0.1
    i <- 1
    
    while(n_rounds_run<200|n_rounds_run==1000) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = parameter_val,
        max_depth = 6,
        min_child_weight = cdi_tox_max_child_bestparams[[outcome]]$min_child_weight,
        subsample = cdi_tox_col_sub_bestparams[[outcome]]$subsample,
        colsample_bytree = cdi_tox_col_sub_bestparams[[outcome]]$colsample_bytree
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = train_matrix,
        nrounds = 1000,
        nfold = 5,
        early_stopping_rounds = 10,
        verbose = 1
      )
      
      best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
      best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
      cv_model$evaluation_log$test_logloss_mean
      
        best_auc <- best_iteration_auc
        best_params <- params
        best_nrounds <- best_iteration_index
      
      
      n_rounds_run <- cv_model$evaluation_log %>% nrow()
      
      if(n_rounds_run<200) {
        
        parameter_val <- parameter_val/2
        
      } else if (n_rounds_run==1000) {
        
        parameter_val <- parameter_val+0.1
        
      }
      
      i <- i+1
      
      if (parameter_val >0.3){
        
        print("Learning rate too high - adjust hyperparameters")
        
        break
        
      }
      
    }
    
    cdi_tox_final_bestparams[[outcome]] <- best_params
    cdi_tox_final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
  }
}

for (outcome in 1:ncol(abx_outcomes)) {
  
  param <- data.frame(cdi_tox_final_bestparams[outcome])
  
  write_csv(param,glue("cdi_tox_final_params_{names(abx_outcomes)[outcome]}.csv"))
  
}

###Saving interim CSVs
write_csv(abx_combined,"abx_combined.csv")
write_csv(ur_abx_combined,"ur_abx_combined.csv")
