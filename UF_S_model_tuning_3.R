###COVERAGE MODEL HYPERPARAMETER TUNING - STAGE 3

set.seed(123)

###Third round of hyperparameter tuning (learning rate)

###Read_in previous parameters
file_names <- glue("col_sub_{combined_antimicrobial_map}.csv")
col_sub_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(col_sub_bestparams)) {
  names(col_sub_bestparams[[i]])[3:8] <- namelist
}
names(col_sub_bestparams) <- combined_antimicrobial_map

file_names <- glue("max_child_{combined_antimicrobial_map}.csv")
max_child_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(max_child_bestparams)) {
  names(max_child_bestparams[[i]])[3:8] <- namelist
}
names(max_child_bestparams) <- combined_antimicrobial_map

final_bestparams <- c()
for (outcome in colnames(urines5_outcomes)) {
  
  best_auc <- 0
  
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
    
    n_rounds_run <- 0
    parameter_val <- 0.1
    i <- 1
    
    while(n_rounds_run<200|n_rounds_run==1000) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = parameter_val,
        max_depth = max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams[[outcome]]$min_child_weight,
        subsample = col_sub_bestparams[[outcome]]$subsample,
        colsample_bytree = col_sub_bestparams[[outcome]]$colsample_bytree
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
      
    }
    
    final_bestparams[[outcome]] <- best_params
    final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
  }
}

for (i in 1:length(final_bestparams)) {
  
  param <- data.frame(final_bestparams[i])
  
  write_csv(param,glue("final_params_{combined_antimicrobial_map[i]}.csv"))
  
}

###Saving interim CSVs
write_csv(urines5_combined,"urines5_combined.csv")
