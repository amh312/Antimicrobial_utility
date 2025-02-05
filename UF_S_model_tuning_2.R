###COVERAGE MODEL HYPERPARAMETER TUNING - STAGE 2
  
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
col_sub_bestparams <- c()
best_auc <- 0
for (outcome in colnames(urines5_outcomes)) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    set.seed(123)
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    missing_cols <- setdiff(selected_columns, colnames(ur_xg_combined))
    ur_xg_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(all_of(selected_columns))), 
                                label = ur_xg_combined[[outcome]])
    
    for (i in 1:nrow(parameter_grid)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams[[outcome]]$min_child_weight,
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
    
    col_sub_bestparams[[outcome]] <- best_params
    
  }
}
for (i in 1:length(col_sub_bestparams)) {
  
  coly <- data.frame(col_sub_bestparams[i])
  
  write_csv(coly,glue("col_sub_{combined_antimicrobial_map[i]}.csv"))
  
}

