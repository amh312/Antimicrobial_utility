##CDI AND TOXICITY HYPERPARAMETER TUNING - STAGE 3

set.seed(123)

###Third round of hyperparameter tuning (learning rate and best iteration)
parameter_list <- c(0.1,0.05,0.01)
best_auc <- 0
cdi_tox_final_bestparams <- c()
for (outcome in colnames(abx_outcomes)) {
  
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
    
    for (i in 1:length(parameter_list)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = parameter_list[i],
        max_depth = cdi_tox_max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = cdi_tox_max_child_bestparams[[outcome]]$min_child_weight,
        subsample = cdi_tox_col_sub_bestparams[[outcome]]$subsample,
        colsample_bytree = cdi_tox_col_sub_bestparams[[outcome]]$colsample_bytree
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = train_matrix,
        nrounds = 1000,
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