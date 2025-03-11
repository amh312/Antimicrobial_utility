##CDI AND TOXICITY HYPERPARAMETER TUNING - STAGE 1

set.seed(123)

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Read-in
probs_df_overall <- read_csv("probs_df_overall.csv")

###Dataset preprocessing
util_probs_df <- probs_df_overall
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()
abx <- data.frame(rbind(train_abx,test_abx))
hadm_key <- hadm %>% select(subject_id,race,language,marital_status) %>% 
  distinct(subject_id,.keep_all=T)
age_key <- pats %>% select(subject_id,anchor_age) %>% distinct(subject_id,.keep_all = T)
dems_key <- left_join(hadm_key,age_key,by="subject_id")
abx <- abx %>% left_join(dems_key,by="subject_id") %>% mutate(anchor_age=as.numeric(anchor_age))
ur_xg <- ur_xg %>% left_join(age_key,by="subject_id")
abx_outcomes <- abx %>%
  select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                     overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
abx_predictors <- abx %>% select(pHADM:age65,prAKI:pDIAB,pCARD:curr_service,pICU:pSEPSIS,ob_freq,highCRP,
                                 pc_dyspnea:SIRS,
                                 ab_name_Ampicillin_Ceftriaxone:ab_name_Ampicillin,race:language,anchor_age)
ur_abx_outcomes <- ur_xg %>%
  select(micro_specimen_id,CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                                       overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
ur_abx_outcomes <- util_probs_df %>% left_join(ur_abx_outcomes,
                                               relationship = 'many-to-one') %>% 
  select(-c(id_no,I:subject_id))
common_columns <- intersect(names(abx_predictors),names(ur_xg))
ur_abx_predictors <- ur_xg %>% select(micro_specimen_id,all_of(common_columns))
ur_abx_predictors <- util_probs_df %>% left_join(ur_abx_predictors,
                                                 relationship = 'many-to-one') %>% 
  select(-c(id_no,I:subject_id)) %>% rename(ab_name="Antimicrobial")
ur_abx_combined <- data.frame(cbind(ur_abx_outcomes,ur_abx_predictors))
set.seed(123)
dummies <- dummyVars(" ~ .", data = abx_predictors)
abx_predictors <- predict(dummies, newdata = abx_predictors)
abx_combined <- as.data.frame(cbind(abx_outcomes, abx_predictors))
dummies <- dummyVars(" ~ .", data = ur_abx_predictors)
ur_abx_predictors <- predict(dummies, newdata = ur_abx_predictors)
ur_abx_combined <- as.data.frame(cbind(ur_abx_outcomes, ur_abx_predictors))
colnames(ur_abx_combined)[grepl("ab_name",colnames(ur_abx_combined))] <- colnames(ur_abx_combined)[grepl("ab_name",colnames(ur_abx_combined))] %>% 
  str_replace("ab_name","ab_name_") %>% 
  str_replace_all("-",".")

###First round of hyperparameter tuning (max tree depth and min child weight)
num_samples <- 10
max_depth_range <- c(2,9)
min_child_weight_range <- c(1, 10)
lhs_sample <- randomLHS(num_samples, 2)
max_depth <- round(lhs_sample[, 1] * (max_depth_range[2] - max_depth_range[1]) + max_depth_range[1])
min_child_weight <- round(lhs_sample[, 2] * (min_child_weight_range[2] - min_child_weight_range[1]) + min_child_weight_range[1])
parameter_grid <- data.frame(max_depth = max_depth, min_child_weight = min_child_weight)
print(parameter_grid)
cdi_tox_max_child_bestparams <- c()

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
    
    for (i in 1:nrow(parameter_grid)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = parameter_grid %>% select(max_depth) %>% dplyr::slice(i) %>% unlist(),
        min_child_weight = parameter_grid %>% select(min_child_weight) %>% dplyr::slice(i) %>% unlist(),
        subsample = 0.8,
        colsample_bytree = 0.8
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
    
    cdi_tox_max_child_bestparams[[outcome]] <- best_params
    
  }
}
for (outcome in 1:ncol(abx_outcomes)) {
  
  maxy <- data.frame(cdi_tox_max_child_bestparams[outcome])
  
  write_csv(maxy,glue("max_child_{names(abx_outcomes)[outcome]}.csv"))
  
}
