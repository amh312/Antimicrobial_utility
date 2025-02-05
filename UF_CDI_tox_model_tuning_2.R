##CDI AND TOXICITY HYPERPARAMETER TUNING - STAGE 2

set.seed(123)

###Second round of hyperparameter tuning (subsample and colsample_bytree)
num_samples <- 10
subsample_range <- c(0.5,1)
colsample_bytree_range <- c(0.5, 1)
lhs_sample <- randomLHS(num_samples, 2)
subsample <- round(lhs_sample[, 1] * (subsample_range[2] - subsample_range[1]) + subsample_range[1],2)
colsample_bytree <- round(lhs_sample[, 2] * (colsample_bytree_range[2] - colsample_bytree_range[1]) + colsample_bytree_range[1],2)
parameter_grid <- data.frame(subsample = subsample, colsample_bytree = colsample_bytree)
print(parameter_grid)
cdi_tox_col_sub_bestparams <- c()
best_auc <- 0
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
    
    for (i in 1:nrow(parameter_grid)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = cdi_tox_max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = cdi_tox_max_child_bestparams[[outcome]]$min_child_weight,
        subsample = parameter_grid %>% select(subsample) %>% dplyr::slice(i) %>% unlist(),
        colsample_bytree = parameter_grid %>% select(colsample_bytree) %>% dplyr::slice(i) %>% unlist()
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = train_matrix,
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
    
    cdi_tox_col_sub_bestparams[[outcome]] <- best_params
    
  }
}
for (outcome in 1:ncol(abx_outcomes)) {
  
  coly <- data.frame(cdi_tox_col_sub_bestparams[outcome])
  
  write_csv(coly,glue("cdi_tox_col_sub_{names(abx_outcomes)[outcome]}.csv"))
  
}