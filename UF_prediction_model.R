#PREDICTION MODEL AND FEATURE TUNING

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

##Read-in

train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
hadm <- read_csv("admissions.csv")
pats <- read_csv("patients.csv")

##Model tuning

###Antimicrobial mapping lists
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                 "MEM","CIP","GEN","SXT","NIT")
ab_singles <- all_singles
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                "GEN","SXT")
iv_ab_singles <- iv_singles
iv_combos <- combn(iv_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_ivs <- c(iv_singles, iv_combos)
oral_singles <- c("AMP","SAM","CIP",
                  "SXT","NIT")
oral_ab_singles <- oral_singles
oral_combos <- combn(oral_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_orals <- c(oral_singles, oral_combos)
access_singles <- c("AMP","SAM","GEN",
                    "SXT","NIT","CZO")
access_combos <- combn(access_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_access <- c(access_singles, access_combos)
watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
watch_combos <- combn(watch_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_watch <- c(watch_singles, watch_combos)
antimicrobial_map <- c(
  "Ampicillin" = "AMP",
  "Ampicillin-sulbactam" = "SAM",
  "Piperacillin-tazobactam" = "TZP",
  "Cefazolin" = "CZO",
  "Ceftriaxone" = "CRO",
  "Ceftazidime" = "CAZ",
  "Cefepime" = "FEP",
  "Meropenem" = "MEM",
  "Ciprofloxacin" = "CIP",
  "Gentamicin" = "GEN",
  "Trimethoprim-sulfamethoxazole" = "SXT",
  "Nitrofurantoin" = "NIT",
  "Vancomycin" = "VAN"
)
map_combinations <- combn(names(antimicrobial_map), 2, simplify = FALSE)
combined_antimicrobial_map <- c(
  antimicrobial_map,
  setNames(
    lapply(map_combinations, function(x) paste(antimicrobial_map[x], collapse = "_")),
    sapply(map_combinations, function(x) paste(x, collapse = "_"))
  )
)
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
fullmap <- combined_antimicrobial_map %>% unlist()
names(fullmap) <- NULL
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
shortmap <- combined_antimicrobial_map %>% unlist()
names(shortmap) <- NULL

###Final preprocessing and dataset train/test splits
urines5 <- urines5 %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                       TRUE~marital_status))
ur_xg <- ur_xg %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                   TRUE~marital_status))
urines5_outcomes <- urines5 %>%
  select(all_of(shortmap))
ur_xg_outcomes <- ur_xg %>%
  select(all_of(shortmap))
urines5_outcomes <- urines5_outcomes %>%
  mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                 ifelse(. == "S" | . == "I", 1, NA))))
ur_xg_outcomes <- ur_xg_outcomes %>%
  mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                 ifelse(. == "S" | . == "I", 1, NA))))
urines5_predictors <- urines5 %>% select(!all_of(fullmap))
ur_xg_predictors <- ur_xg %>%
  select(any_of(names(urines5_predictors)))
set.seed(123)
dummies <- dummyVars(" ~ .", data = urines5_predictors)
urines5_predictors <- predict(dummies, newdata = urines5_predictors)
dummies2 <- dummyVars(" ~ .", data = ur_xg_predictors)
ur_xg_predictors <- predict(dummies2, newdata = ur_xg_predictors)
urines5_combined <- as.data.frame(cbind(urines5_outcomes, urines5_predictors))
ur_xg_combined <- as.data.frame(cbind(ur_xg_outcomes, ur_xg_predictors))

###First round of hyperparameter tuning (max tree depth and min child weight)
num_samples <- 10
max_depth_range <- c(2,9)
min_child_weight_range <- c(1, 10)
lhs_sample <- randomLHS(num_samples, 2)
max_depth <- round(lhs_sample[, 1] * (max_depth_range[2] - max_depth_range[1]) + max_depth_range[1])
min_child_weight <- round(lhs_sample[, 2] * (min_child_weight_range[2] - min_child_weight_range[1]) + min_child_weight_range[1])
parameter_grid <- data.frame(max_depth = max_depth, min_child_weight = min_child_weight)
print(parameter_grid)
max_child_bestparams <- c()
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

###Third round of hyperparameter tuning (learning rate)
parameter_list <- c(0.1,0.05,0.01,0.001)
best_auc <- 0
final_bestparams <- c()
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
    
    for (i in 1:length(parameter_list)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = parameter_list[i],
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
      if (best_iteration_auc > best_auc) {
        best_auc <- best_iteration_auc
        best_params <- params
        best_nrounds <- best_iteration_index
      }
      
    }
    
    final_bestparams[[outcome]] <- best_params
    
  }
}
for (i in 1:length(final_bestparams)) {
  
  param <- data.frame(final_bestparams[i])
  
  write_csv(param,glue("final_params_{combined_antimicrobial_map[i]}.csv"))
  
}

###Best iteration
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
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = final_bestparams[[outcome]]$eta,
      max_depth = final_bestparams[[outcome]]$max_depth,
      min_child_weight = final_bestparams[[outcome]]$min_child_weight,
      subsample = final_bestparams[[outcome]]$subsample,
      colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
    )
    
    print("Running CV...")
    
    cv_model <- xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = 1000,
      nfold = 5,
      early_stopping_rounds = 50,
      verbose = 1,
    )
    
    final_bestparams[[outcome]]$best_nrounds <- cv_model$best_iteration

  }
}
for (i in 1:length(final_bestparams)) {
  
  param <- data.frame(final_bestparams[i])
  
  write_csv(param,glue("final_params_{combined_antimicrobial_map[i]}.csv"))
  
}

###Saving interim CSVs
write_csv(urines5_combined,"urines5_combined.csv")
write_csv(urines_ref,"urines_ref_combined.csv")

###Model training, testing and microsimulation probability predictions
test_probs_df <- data.frame(matrix(nrow=floor(nrow(urines5_combined)*0.2),ncol=0))
micro_probs_df <- data.frame(matrix(nrow=nrow(ur_xg_combined),ncol=0))
aucs <- data.frame(matrix(nrow=1,ncol=0))
shap_summary_tables <- list()
metrics_list <- list() 
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
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = final_bestparams[[outcome]]$eta,
      max_depth = final_bestparams[[outcome]]$max_depth,
      min_child_weight = final_bestparams[[outcome]]$min_child_weight,
      subsample = final_bestparams[[outcome]]$subsample,
      colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
    )
    
    print("Training...")
    
    xgb_model <- xgb.train(
      params = params,
      data = train_matrix,
      nrounds = final_bestparams[[outcome]]$best_nrounds
    )
    
    print("Shapping...")
    
    shap_values <- predict(xgb_model, newdata = train_matrix, predcontrib = TRUE)
    shap_df <- as.data.frame(shap_values)
    shap_df <- shap_df[, -ncol(shap_df)]
    shap_summary <- data.frame(
      Feature = colnames(shap_df),
      MeanAbsSHAP = colMeans(abs(shap_df))
    )
    shap_summary <- shap_summary %>% filter(MeanAbsSHAP!=0)
    
    #Feature selection
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(shap_summary %>% pull(Feature))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(shap_summary %>% pull(Feature))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(shap_summary %>% pull(Feature))), 
                                label = ur_xg_combined[[outcome]])
    
    #Run again with selected features
    xgb_model <- xgb.train(
      params = params,
      data = train_matrix,
      nrounds = final_bestparams[[outcome]]$best_nrounds,
    )
    
    pred_prob_test <- predict(xgb_model, newdata = test_matrix)
    roc_result <- roc(urines5Test[[outcome]], pred_prob_test)
    auc_value <- auc(roc_result)
    print(paste("AUC-ROC:", auc_value))
    aucs[[outcome]] <- auc_value
    pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
    plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
    dev.off()
    pred_prob_micro <- predict(xgb_model, newdata = micro_matrix)
    
    test_probs_df[[outcome]] <- pred_prob_test
    micro_probs_df[[outcome]] <- pred_prob_micro
    
    shap_summary <- shap_summary[order(-shap_summary$MeanAbsSHAP), ]
    
    shap_summary_tables[[outcome]] <- shap_summary
    
    pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
    actual_test_class <- urines5Test[[outcome]]
    
    confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
    accuracy <- confusion$overall['Accuracy']
    precision <- confusion$byClass['Precision']
    recall <- confusion$byClass['Recall']
    f1_score <- 2 * (precision * recall) / (precision + recall)
    
    metrics_list[[outcome]] <- list(
      AUC = auc_value,
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1_Score = f1_score
    )
    
  }
}
for (i in 1:length(shap_summary_tables)) {
  
  shappy <- data.frame(shap_summary_tables[i]) %>% mutate(across(2, ~ round(., 3))) %>% filter(if_any(2, ~ . != 0))
  
  colnames(shappy) <- c("Feature","Variable")
  
  write_csv(shappy,glue("SHAP_{combined_antimicrobial_map[i]}.csv"))
  
}
for (i in 1:length(metrics_list)) {
  
  metricky <- data.frame(metrics_list[i])
  
  write_csv(metricky,glue("metrics_{combined_antimicrobial_map[i]}.csv"))
  
}

###Microsimulation dataframe
replace_values <- function(column, map) {
  flipped_map <- setNames(names(map), map)
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(flipped_map)) flipped_map[[x]] else x)
}
micro_probs_df$micro_specimen_id <- ur_xg$micro_specimen_id
micro_probs_df$subject_id <- ur_xg$subject_id
probs_df_overall <- micro_probs_df %>% melt(id.vars=c("micro_specimen_id","subject_id")) %>% 
  mutate(variable=replace_values(variable, combined_antimicrobial_map)) %>% 
  rename(S="value",Antimicrobial="variable") %>% mutate(I=NA,R=1-S,NT=NA)
probs_df_overall <- probs_df_overall %>% relocate(micro_specimen_id,.after = "NT") %>% 
  relocate(subject_id,.after="micro_specimen_id") %>% 
  relocate(I,.before="S") %>% mutate(id_no=1:nrow(probs_df_overall)) %>% 
  relocate(id_no,.before="Antimicrobial")
write_csv(probs_df_overall,"probs_df_overall.csv")

##CDI and toxicity prediction models

###Dataset preprocessing
util_probs_df <- probs_df_overall
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()
abx <- data.frame(rbind(train_abx,test_abx))
hadm_key <- hadm %>% select(subject_id,race,language,marital_status) %>% 
  distinct(subject_id,.keep_all=T)
age_key <- pats %>% select(subject_id,anchor_age) %>% distinct(subject_id,.keep_all = T)
dems_key <- left_join(hadm_key,age_key,by="subject_id")
abx <- abx %>% left_join(dems_key,by="subject_id")
ur_xg <- ur_xg %>% left_join(age_key,by="subject_id")
abx_outcomes <- abx %>%
  select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                     overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
abx_predictors <- abx %>% select(pHADM:age65,prAKI:pDIAB,pCARD:pSEPSIS,temperature:dbp,pc_dyspnea:pc_fever,
                                 abx_name_Ampicillin_Ceftriaxone:ob_freq,highCRP,race:anchor_age)
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
  select(-c(id_no,I:subject_id)) %>% rename(abx_name="Antimicrobial")
ur_abx_combined <- data.frame(cbind(ur_abx_outcomes,ur_abx_predictors))
set.seed(123)
dummies <- dummyVars(" ~ .", data = abx_predictors)
abx_predictors <- predict(dummies, newdata = abx_predictors)
abx_combined <- as.data.frame(cbind(abx_outcomes, abx_predictors))
dummies <- dummyVars(" ~ .", data = ur_abx_predictors)
ur_abx_predictors <- predict(dummies, newdata = ur_abx_predictors)
ur_abx_combined <- as.data.frame(cbind(ur_abx_outcomes, ur_abx_predictors))
colnames(ur_abx_combined)[grepl("abx_name",colnames(ur_abx_combined))] <- colnames(ur_abx_combined)[grepl("abx_name",colnames(ur_abx_combined))] %>% 
  str_replace("abx_name","abx_name_") %>% 
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

###Third round of hyperparameter tuning (learning rate)
parameter_list <- c(0.1,0.05,0.01,0.001)
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
    
  }
}
for (outcome in 1:ncol(abx_outcomes)) {
  
  param <- data.frame(cdi_tox_final_bestparams[outcome])
  
  write_csv(param,glue("cdi_tox_final_params_{names(abx_outcomes)[outcome]}.csv"))
  
}

###Best iteration
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
    
    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = cdi_tox_final_bestparams[[outcome]]$eta,
      max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
      min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
      subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
      colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
    )
    
    print("Running CV...")
    
    cv_model <- xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = 1000,
      nfold = 5,
      early_stopping_rounds = 50,
      verbose = 1,
    )
    
    cdi_tox_final_bestparams[[outcome]]$best_nrounds <- cv_model$best_iteration
    
  }
}
for (outcome in 1:ncol(abx_outcomes)) {
  
  param <- data.frame(cdi_tox_final_bestparams[outcome])
  
  write_csv(param,glue("cdi_tox_final_params_{names(abx_outcomes)[outcome]}.csv"))
  
}

###Read-in chosen model parameters (CDI/toxicity prediction)
file_names <- glue("cdi_tox_final_params_{colnames(abx_outcomes)}.csv")
cdi_tox_final_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(cdi_tox_final_bestparams)) {
  names(cdi_tox_final_bestparams[[i]])[3:8] <- namelist
}
names(cdi_tox_final_bestparams) <- colnames(abx_outcomes)

###Saving interim CSVs
write_csv(abx_combined,"abx_combined.csv")
write_csv(ur_abx_combined,"ur_abx_combined.csv")

###CDI final model run
set.seed(123)
trainIndex <- createDataPartition(abx_combined[['CDI']], p = 0.8, list = FALSE, times = 1)
abxTrain <- abx_combined[trainIndex, ]
abxTest <- abx_combined[-trainIndex, ]
predictor_columns <- colnames(abx_predictors)
selected_columns <- intersect(predictor_columns, colnames(abxTrain))
missing_cols <- setdiff(selected_columns, colnames(ur_abx_combined))
ur_abx_combined[missing_cols] <- 0
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                            label = abxTrain[['CDI']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                           label = abxTest[['CDI']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(all_of(selected_columns))), 
                            label = ur_abx_combined[['CDI']])

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = cdi_tox_final_bestparams[['CDI']]$eta,
  max_depth = cdi_tox_final_bestparams[['CDI']]$max_depth,
  min_child_weight = cdi_tox_final_bestparams[['CDI']]$min_child_weight,
  subsample = cdi_tox_final_bestparams[['CDI']]$subsample,
  colsample_bytree = cdi_tox_final_bestparams[['CDI']]$colsample_bytree
)

print("Training...")

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = cdi_tox_final_bestparams[['CDI']]$best_nrounds
)

print("Shapping...")

shap_values <- predict(xgb_model, newdata = train_matrix, predcontrib = TRUE)
shap_df <- as.data.frame(shap_values)
shap_df <- shap_df[, -ncol(shap_df)]
shap_summary <- data.frame(
  Feature = colnames(shap_df),
  MeanAbsSHAP = colMeans(abs(shap_df))
)

cdi_shap_summary <- shap_summary %>% filter(MeanAbsSHAP!=0)

###CDI feature selection
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(shap_summary %>% pull(Feature))), 
                            label = abxTrain[['CDI']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(shap_summary %>% pull(Feature))), 
                           label = abxTest[['CDI']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(shap_summary %>% pull(Feature))), 
                            label = ur_abx_combined[['CDI']])

###CDI run again with selected features
xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = cdi_tox_final_bestparams[['CDI']]$best_nrounds,
)

pred_prob_test <- predict(xgb_model, newdata = test_matrix)
roc_result <- roc(abxTest[['CDI']], pred_prob_test)
cdi_auc_value <- auc(roc_result)
print(paste("AUC-ROC:", cdi_auc_value))
pdf(glue("CDI_xg_roc.pdf"), width = 10, height = 10)
plot(roc_result, main = glue("CDI ROC Curve"), col = "blue")
dev.off()
cdi_util_probs <- predict(xgb_model, newdata = micro_matrix)

pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
actual_test_class <- abxTest[['CDI']]

cdi_confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
cdi_accuracy <- cdi_confusion$overall['Accuracy']
cdi_precision <- cdi_confusion$byClass['Precision']
cdi_recall <- cdi_confusion$byClass['Recall']
cdi_f1_score <- 2 * (cdi_precision * cdi_recall) / (cdi_precision + cdi_recall)
probs_df_overall$prob_CDI <- cdi_util_probs
write_csv(probs_df_overall,"probs_df_overall.csv")

write_csv(cdi_shap_summary,glue("SHAP_CDI.csv"))

cdi_metrics_list <- data.frame(
  AUC = cdi_auc_value,
  Accuracy = cdi_accuracy,
  Precision = cdi_precision,
  Recall = cdi_recall,
  F1_Score = cdi_f1_score
)

write_csv(cdi_metrics_list,glue("metrics_CDI.csv"))

###Toxicity data partitioning
set.seed(123)
trainIndex <- createDataPartition(abx_combined[['overall_tox']], p = 0.8, list = FALSE, times = 1)
abxTrain <- abx_combined[trainIndex, ]
abxTest <- abx_combined[-trainIndex, ]
predictor_columns <- colnames(abx_predictors)
selected_columns <- intersect(predictor_columns, colnames(abxTrain))
missing_cols <- setdiff(selected_columns, colnames(ur_abx_combined))
ur_abx_combined[missing_cols] <- 0
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                            label = abxTrain[['overall_tox']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                           label = abxTest[['overall_tox']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(all_of(selected_columns))), 
                            label = ur_abx_combined[['overall_tox']])

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = cdi_tox_final_bestparams[['overall_tox']]$eta,
  max_depth = cdi_tox_final_bestparams[['overall_tox']]$max_depth,
  min_child_weight = cdi_tox_final_bestparams[['overall_tox']]$min_child_weight,
  subsample = cdi_tox_final_bestparams[['overall_tox']]$subsample,
  colsample_bytree = cdi_tox_final_bestparams[['overall_tox']]$colsample_bytree
)

print("Training...")

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = cdi_tox_final_bestparams[['overall_tox']]$best_nrounds
)

print("Shapping...")

shap_values <- predict(xgb_model, newdata = train_matrix, predcontrib = TRUE)
shap_df <- as.data.frame(shap_values)
shap_df <- shap_df[, -ncol(shap_df)]
shap_summary <- data.frame(
  Feature = colnames(shap_df),
  MeanAbsSHAP = colMeans(abs(shap_df))
)
overall_tox_shap_summary <- shap_summary %>% filter(MeanAbsSHAP!=0)

###Toxicity feature selection
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(shap_summary %>% pull(Feature))), 
                            label = abxTrain[['overall_tox']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(shap_summary %>% pull(Feature))), 
                           label = abxTest[['overall_tox']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(shap_summary %>% pull(Feature))), 
                            label = ur_abx_combined[['overall_tox']])

###Toxicity run again with selected features
xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = cdi_tox_final_bestparams[['overall_tox']]$best_nrounds,
)

pred_prob_test <- predict(xgb_model, newdata = test_matrix)
roc_result <- roc(abxTest[['overall_tox']], pred_prob_test)
overall_tox_auc_value <- auc(roc_result)
print(paste("AUC-ROC:", overall_tox_auc_value))
pdf(glue("overall_tox_xg_roc.pdf"), width = 10, height = 10)
plot(roc_result, main = glue("overall_tox ROC Curve"), col = "blue")
dev.off()
overall_tox_util_probs <- predict(xgb_model, newdata = micro_matrix)

pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
actual_test_class <- abxTest[['overall_tox']]

overall_tox_confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
overall_tox_accuracy <- overall_tox_confusion$overall['Accuracy']
overall_tox_precision <- overall_tox_confusion$byClass['Precision']
overall_tox_recall <- overall_tox_confusion$byClass['Recall']
overall_tox_f1_score <- 2 * (overall_tox_precision * overall_tox_recall) / (overall_tox_precision + overall_tox_recall)
probs_df_overall$prob_tox <- overall_tox_util_probs
write_csv(probs_df_overall,"probs_df_overall.csv")
write_csv(overall_tox_shap_summary,glue("SHAP_toxicity.csv"))

tox_metrics_list <- data.frame(
  AUC = overall_tox_auc_value,
  Accuracy = overall_tox_accuracy,
  Precision = overall_tox_precision,
  Recall = overall_tox_recall,
  F1_Score = overall_tox_f1_score
)

write_csv(tox_metrics_list,glue("metrics_toxicity.csv"))

