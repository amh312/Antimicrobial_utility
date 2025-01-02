#PREDICTION MODEL AND FEATURE TUNING

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Value replacement
replace_values <- function(column, map) {
  flipped_map <- setNames(names(map), map)
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(flipped_map)) flipped_map[[x]] else x)
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
                 "MEM","CIP","GEN","SXT","NIT","VAN")
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
    final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
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
roc_plots <- list()
confidence_biglist <- list()
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
    roc_plot <- ggroc(roc_result) + 
      ggtitle(glue("{replace_values(outcome,combined_antimicrobial_map)}\nROC Curve"))+
      theme_minimal()+
      labs(x = "1 - Specificity", y = "Sensitivity")+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )+
      geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
                   color = "grey", linetype = "dashed")
    
    roc_plots[[outcome]] <- roc_plot
    
    ggsave(filename = glue("{outcome}_xg_roc.pdf"), plot = roc_plot, width = 10, height = 10)

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
    
    #Bootstrapping for confidence intervals
    
    n_bootstraps <- 1000
    calculate_metrics <- function(df, indices) {
      y_true_boot <- df$y_true[indices]
      y_scores_boot <- df$y_scores[indices]
      y_pred_boot <- df$y_pred[indices]
      
      confusion <- confusionMatrix(factor(y_pred_boot), factor(y_true_boot))
      accuracy <- confusion$overall['Accuracy']
      precision <- confusion$byClass['Precision']
      recall <- confusion$byClass['Recall']
      f1_score <- 2 * (precision * recall) / (precision + recall)
      auroc <- auc(y_true_boot, y_scores_boot)
      
      return(c(auroc = auroc, precision = precision, recall = recall, accuracy = accuracy, f1 = f1_score))
    }
    
    df_ys <- data.frame(y_true = actual_test_class, y_pred = pred_test_class,y_scores = pred_prob_test)
    
    boot_results <- boot(data = df_ys, statistic = calculate_metrics, R = n_bootstraps)
    
    confidence_intervals <- lapply(1:5, function(i) boot.ci(boot_results, type = "perc", index = i))
    
    names(confidence_intervals) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")
    
    confidence_biglist[[outcome]] <- confidence_intervals
    
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
ci_df <- data.frame(matrix(nrow=length(urines5_outcomes),ncol=5))
colnames(ci_df) <- c("AUROC_CI","Precision_CI","Recall_CI","F1_CI","Accuracy_CI")

for (i in 1: length(confidence_biglist)) {
  
  ci_df$AUROC_CI[i] <- glue(" ({round(confidence_biglist[[i]]$AUC$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})")
  ci_df$Precision_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Precision$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})")
  ci_df$Recall_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Recall$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})")
  ci_df$Accuracy_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Accuracy$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})")
  ci_df$F1_CI[i] <- glue(" ({round(confidence_biglist[[i]]$F1_Score$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})")
  
}

write_csv(ci_df,"interim_ci_df.csv")

###Microsimulation dataframe
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

###Third round of hyperparameter tuning (learning rate and best iteration)
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
    cdi_tox_final_bestparams[[outcome]]$best_nrounds <- best_nrounds
    
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
roc_plot <- ggroc(roc_result) + 
  ggtitle(glue("CDI\nROC Curve"))+
  theme_minimal()+
  labs(x = "1 - Specificity", y = "Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
               color = "grey", linetype = "dashed")

roc_plots[['CDI']] <- roc_plot

ggsave(filename = glue("CDI_xg_roc.pdf"), plot = roc_plot, width = 10, height = 10)
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

#Bootstrapping for confidence intervals

n_bootstraps <- 1000
calculate_metrics <- function(df, indices) {
  y_true_boot <- df$y_true[indices]
  y_scores_boot <- df$y_scores[indices]
  y_pred_boot <- df$y_pred[indices]
  
  confusion <- confusionMatrix(factor(y_pred_boot), factor(y_true_boot))
  accuracy <- confusion$overall['Accuracy']
  precision <- confusion$byClass['Precision']
  recall <- confusion$byClass['Recall']
  f1_score <- 2 * (precision * recall) / (precision + recall)
  auroc <- auc(y_true_boot, y_scores_boot)
  
  return(c(auroc = auroc, precision = precision, recall = recall, accuracy = accuracy, f1 = f1_score))
}

df_ys <- data.frame(y_true = actual_test_class, y_pred = pred_test_class,y_scores = pred_prob_test)

boot_results <- boot(data = df_ys, statistic = calculate_metrics, R = n_bootstraps)

confidence_intervals <- lapply(1:5, function(i) boot.ci(boot_results, type = "perc", index = i))

names(confidence_intervals) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")

confidence_biglist[['CDI']] <- confidence_intervals

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
roc_plot <- ggroc(roc_result) + 
  ggtitle(glue("Toxicity\nROC Curve"))+
  theme_minimal()+
  labs(x = "1 - Specificity", y = "Sensitivity")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
               color = "grey", linetype = "dashed")

roc_plots[['overall_tox']] <- roc_plot

ggsave(filename = glue("toxicity_xg_roc.pdf"), plot = roc_plot, width = 10, height = 10)
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

#Bootstrapping for confidence intervals

n_bootstraps <- 1000
calculate_metrics <- function(df, indices) {
  y_true_boot <- df$y_true[indices]
  y_scores_boot <- df$y_scores[indices]
  y_pred_boot <- df$y_pred[indices]
  
  confusion <- confusionMatrix(factor(y_pred_boot), factor(y_true_boot))
  accuracy <- confusion$overall['Accuracy']
  precision <- confusion$byClass['Precision']
  recall <- confusion$byClass['Recall']
  f1_score <- 2 * (precision * recall) / (precision + recall)
  auroc <- auc(y_true_boot, y_scores_boot)
  
  return(c(auroc = auroc, precision = precision, recall = recall, accuracy = accuracy, f1 = f1_score))
}

df_ys <- data.frame(y_true = actual_test_class, y_pred = pred_test_class,y_scores = pred_prob_test)

boot_results <- boot(data = df_ys, statistic = calculate_metrics, R = n_bootstraps)

confidence_intervals <- lapply(1:5, function(i) boot.ci(boot_results, type = "perc", index = i))

names(confidence_intervals) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")

confidence_biglist[['overall_tox']] <- confidence_intervals

###Save all metric confidence intervals to CSV
ci_df <- data.frame(matrix(nrow=length(confidence_biglist),ncol=5))
colnames(ci_df) <- c("AUROC_CI","Precision_CI","Recall_CI","F1_CI","Accuracy_CI")

for (i in 1: length(confidence_biglist)) {
  
  ci_df$AUROC_CI[i] <- ifelse(!is.na(confidence_biglist[[i]]$AUC),glue(" ({round(confidence_biglist[[i]]$AUC$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})"),"NA")
  ci_df$Precision_CI[i] <- ifelse(!is.na(confidence_biglist[[i]]$Precision),glue(" ({round(confidence_biglist[[i]]$Precision$percent[4],3)}-{round(confidence_biglist[[i]]$Precision$percent[5],3)})"),"NA")
  ci_df$Recall_CI[i] <- ifelse(!is.na(confidence_biglist[[i]]$Recall),glue(" ({round(confidence_biglist[[i]]$Recall$percent[4],3)}-{round(confidence_biglist[[i]]$Recall$percent[5],3)})"),"NA")
  ci_df$Accuracy_CI[i] <- ifelse(!is.na(confidence_biglist[[i]]$Accuracy),glue(" ({round(confidence_biglist[[i]]$Accuracy$percent[4],3)}-{round(confidence_biglist[[i]]$Accuracy$percent[5],3)})"),"NA")
  ci_df$F1_CI[i] <- ifelse(!is.na(confidence_biglist[[i]]$F1_Score),glue(" ({round(confidence_biglist[[i]]$F1_Score$percent[4],3)}-{round(confidence_biglist[[i]]$F1_Score$percent[5],3)})"),"NA")
  
}

ci_df$Model <- names(confidence_biglist)

write_csv(ci_df,"final_ci_df.csv")
ci_filter_list <- c(all_singles,"CDI","overall_tox")
ci_singles_table <- ci_df %>% filter(Model %in% ci_filter_list) %>% 
  mutate(Model=case_when(grepl("tox",Model)~"Toxicity",TRUE~Model))
ci_combos_table <- ci_df %>% filter(!Model %in% ci_filter_list)
ci_singles_table <- ci_singles_table %>% mutate(Model = case_when(
  Model %in% combined_antimicrobial_map ~ replace_values(
    Model, combined_antimicrobial_map
  ),TRUE~Model
))
ci_combos_table <- ci_combos_table %>% mutate(Model = case_when(
  Model %in% combined_antimicrobial_map ~ replace_values(
    Model, combined_antimicrobial_map
  ),TRUE~Model
)) %>% mutate(Model = str_replace_all(Model,"_"," & "))

write_csv(ci_singles_table,"ci_singles_table.csv")
write_csv(ci_combos_table,"ci_combos_table.csv")

##Save all ROCs to grid PDFs
roc_plots_main <- c(roc_plots[1:13],roc_plots[38:39])
roc_plots_other <- roc_plots[14:37]

main_grid_plot <- plot_grid(plotlist = roc_plots_main, ncol = 3)
other_grid_plot <- plot_grid(plotlist = roc_plots_other, ncol = 3)

ggsave(glue("roc_plots_main.pdf"), plot = main_grid_plot, width = 15, height = 30)
ggsave(glue("roc_plots_other.pdf"), plot = other_grid_plot, width = 15, height = 30)

##Feature contributions dataframe
other_outcome_map <- list("CDI","toxicity")
names(other_outcome_map) <- c("CDI","Toxicity")
full_outcome_map <- c(combined_antimicrobial_map,other_outcome_map)
feat_table <- data.frame(matrix(nrow=0,ncol=3))
for (i in 1:length(full_outcome_map)) {
  
  to_bind <- read_csv(glue("SHAP_{full_outcome_map[i]}.csv"))
  colnames(to_bind)[2] <- "Shapley value"
  to_bind$Model <- names(full_outcome_map)[i]
  to_bind <- to_bind %>% relocate(Model,.before = 1)
  feat_table <- tibble(rbind(feat_table,to_bind))
  
}

pc <- "pc_"
true <- "TRUE"
false <- "FALSE"
provider_id <- "provider_id"
pr <- "^pr"
p <- "^p"
r <- "r$"
rx <- "rx$"
d7 <- "^d7"
prev_week <- "Previous week"

feat_table <- feat_table %>% mutate(Feature = case_when(
  grepl("pc",Feature)&grepl("TRUE",Feature)~glue("Admission complaint of {Feature %>% str_remove(pc)}"),
  grepl("pc",Feature)&grepl("FALSE",Feature)~glue("admission complaint of {Feature %>% str_remove(pc)}"),
  grepl("^provider",Feature)~glue("Provider ID {Feature %>% str_remove(provider_id)}"),
  grepl("^pr",Feature)~glue("Previous {Feature %>% str_remove(pr)}"),
  grepl("^p",Feature)~glue("Previous {Feature %>% str_remove(p)}"),
  grepl("^d7",Feature)~glue("Previous week {Feature %>% str_remove(d7)}"),
  TRUE~Feature
)) %>% mutate(Feature=case_when(
  grepl("FALSE$",Feature)~glue("No {Feature %>% str_remove(false)}"),
                                grepl("TRUE$",Feature)~glue("{Feature %>% str_remove(true)}"),
                                TRUE~Feature
)) %>% 
  mutate(Feature = case_when(
  grepl("r$",Feature)&!grepl("(insur|Rest|fever)",Feature)~glue("{Feature %>% str_remove(r)} resistance in the last year"),
  grepl("rx$",Feature)~glue("{Feature %>% str_remove(rx)} treatment in the last year"),
  grepl("ICD_",Feature)~str_replace(Feature,"ICD_","ICD diagnosis category "),
  grepl("PROC_",Feature)~str_replace(Feature,"PROC_","ICD procedure category "),
  grepl("NUTR",Feature)~str_replace(Feature,"NUTR","nutrition input in the last year"),
  grepl("CATH",Feature)~str_replace(Feature,"CATH","urinary catheter in the last 28 days"),
  grepl("Surg",Feature)~str_replace(Feature,"Surg","surgery in the least year"),
  grepl("DISC",Feature)~str_replace(Feature,"DISC","discharge in the last 28 days"),
  grepl("Restr",Feature)~str_replace(Feature,"Restr","need for restraints in the last year"),
  grepl("NH",Feature)~str_replace(Feature,"NH","discharge to nursing care"),
  grepl("HADM",Feature)~str_replace(Feature,"HADM","hospital admission in the last year"),
  grepl("Physio",Feature)~str_replace(Feature,"Physio","physiotherapy input in the last year"),
  grepl("Social",Feature)~str_replace(Feature,"Social","social worker input in the last year"),
  grepl("Obese",Feature)~str_replace(Feature,"Obese","obesity recorded in the last 3 years"),
  grepl("TPN",Feature)~str_replace(Feature,"TPN","parenteral nutrition in the last year"),
  grepl("Overweight",Feature)~str_replace(Feature,"Overweight","overweight BMI in the last 3 years"),
  grepl("DNR",Feature)~str_replace(Feature,"DNR","'do not resuscitate' order in the last year"),
  grepl("OT",Feature)&!grepl("OTHER",Feature)~str_replace(Feature,"OT","occupational therapy input in the last year"),
  grepl("ICU",Feature)~str_replace(Feature,"ICU","intensive care admission in the last 28 days"),
  grepl("NGT",Feature)~str_replace(Feature,"NGT","nasogastric tube in the last 28 days"),
  grepl("Neph",Feature)~str_replace(Feature,"Neph","nephrostomy insertion in the last year"),
  grepl("Underweight",Feature)~str_replace(Feature,"Underweight","underweight BMI in the last 3 years"),
  grepl("DIAB",Feature)~str_replace(Feature,"DIAB","diabetes diagnosis in the last year"),
  grepl("SEPSIS",Feature)~str_replace(Feature,"SEPSIS","sepsis in the last year"),
  grepl("CDI",Feature)~str_replace(Feature,"CDI","CDI in the last year"),
  grepl("AKI",Feature)~str_replace(Feature,"AKI","acute kidney injury"),
  grepl("DIAB",Feature)~str_replace(Feature,"DIAB","diabetes diagnosis"),
  grepl("CARD",Feature)~str_replace(Feature,"CARD","cardiovascular disease"),
  grepl(" CA$",Feature)~str_replace(Feature," CA"," cancer "),
  grepl("CKD",Feature)~str_replace(Feature,"CKD","chronic kidney disease"),
  grepl("CVA",Feature)~str_replace(Feature,"CVA","stroke"),
  grepl("standard_age",Feature)~str_replace(Feature,"standard_age","Age group"),
  grepl("ob_freq",Feature)~str_replace(Feature,"ob_freq","Observation frequency on admission"),
  grepl("language",Feature)~str_replace(Feature,"language","Language: "),
  grepl("admission_location",Feature)~str_replace(Feature,"admission_location","Admitted from "),
  grepl("insurance",Feature)~str_replace(Feature,"insurance","Insurance status: "),
  grepl("marital_status",Feature)~str_replace(Feature,"marital_status","Marital status: "),
  grepl("curr_service_",Feature)~str_replace(Feature,"curr_service_","Current inpatient service provider: "),
  grepl("curr_service",Feature)~str_replace(Feature,"curr_service","Current inpatient service provider: "),
  grepl("race",Feature)~str_replace(Feature,"race","Race: "),
  grepl("age65",Feature)~str_replace(Feature,"age65","Over age 65"),
  grepl("highCRP",Feature)~str_replace(Feature,"highCRP","Elevated CRP on admission"),
  grepl("abx_name_",Feature)~str_replace(Feature,"abx_name_","Antibiotic(s): "),
  grepl("anchor_age",Feature)~str_replace(Feature,"anchor_age","Age"),
  grepl("`Alkaline Phosphatase`",Feature)~str_replace(Feature,"`Alkaline Phosphatase`","Alkaline Phosphatase"),
  grepl("`Red blood Cells`",Feature)~str_replace(Feature,"`Red blood Cells`","Red blood Cell Count"),
  grepl("`White Blood Cells`",Feature)~str_replace(Feature,"`White Blood Cells`","Leukocyte Count"),
  grepl("`Platelet Count`",Feature)~str_replace(Feature,"`Platelet Count`","Platelet Count"),
  grepl("`Bilirubin, Total`",Feature)~str_replace(Feature,"`Bilirubin, Total`","Bilirubin"),
  grepl("`Lactate Dehydrogenase`",Feature)~str_replace(Feature,"`Lactate Dehydrogenase`","Lactate Dehydrogenase"),
  grepl("`Alanine Aminotransferase`",Feature)~str_replace(Feature,"`Alanine Aminotransferase`","Alanine Aminotransferase"),
  grepl("temperat",Feature)~str_replace(Feature,"temperature","Temperature on admission"),
  grepl("heartrate",Feature)~str_replace(Feature,"heartrate","Heart rate on admission"),
  grepl("resprate",Feature)~str_replace(Feature,"resprate","Respiratory rate on admission"),
  grepl("abnormalWCC",Feature)~str_replace(Feature,"abnormalWCC","Abnormal white cell count on admission"),
  grepl("o2sat",Feature)~str_replace(Feature,"o2sat","Oxygen saturation on admission"),
  grepl("temperat",Feature)~str_replace(Feature,"temperature","Temperature on admission"),
  grepl("arrival_transport",Feature)~str_replace(Feature,"arrival_transport","Means of arrival "),
  grepl("MALE",Feature)~str_replace(Feature,"MALE","Male"),
  grepl("^sbp",Feature)~str_replace(Feature,"sbp","Systolic blood pressure on admission"),
  grepl("^dbp",Feature)~str_replace(Feature,"dbp","Diastolic blood pressure on admission"),
  grepl("^acuity",Feature)~str_replace(Feature,"acuity","Mental acuity score on admission"),
  grepl("disposition",Feature)~str_replace(Feature,"disposition","Transfer location: "),
  grepl("^`",Feature)~str_replace(Feature,"^`","Admitted on "),
  grepl("Multivitamins",Feature)~str_replace(Feature,"Multivitamins","Admitted on multivitamins"),
  grepl("lowbp",Feature)~str_replace(Feature,"lowbp","hypotension"),
  grepl("flankpain",Feature)~str_replace(Feature,"flankpain","flank pain"),
  grepl("anemia",Feature)~str_replace(Feature,"anemia","anaemia"),
  grepl("prbleed",Feature)~str_replace(Feature,"prbleed","PR bleeding"),
  grepl("diarvom",Feature)~str_replace(Feature,"diarvom","diarrhoea & vomiting"),
  grepl("dyspnea",Feature)~str_replace(Feature,"dyspnea","dyspnoea"),
  grepl("diarrhea",Feature)~str_replace(Feature,"diarrhea","diarrhoea"),
  grepl("backpain",Feature)~str_replace(Feature,"backpain","back pain"),
  grepl("Hyd",Feature)~str_replace(Feature,"Hyd","hydration therapy"),
  TRUE~Feature
)) %>% 
  mutate(
    Feature=case_when(
      grepl("\\?",Feature)~str_replace(Feature,"\\?","UNKNOWN"),
      grepl("^No `",Feature)~str_replace(Feature,"No `","Not admitted on `"),
      TRUE~Feature
    )
  ) %>% 
  mutate(
    Feature=case_when(
      grepl("`",Feature)~str_replace_all(Feature,"`",""),
      TRUE~Feature
    )
  ) %>% mutate(
    Feature=case_when(
      grepl("AMP",Feature)~str_replace(Feature,"AMP",ab_name("AMP")),
      grepl("SAM",Feature)~str_replace(Feature,"SAM",ab_name("SAM")),
      grepl("TZP",Feature)~str_replace(Feature,"TZP",ab_name("TZP")),
      grepl("CZO",Feature)~str_replace(Feature,"CZO",ab_name("CZO")),
      grepl("CRO",Feature)~str_replace(Feature,"CRO",ab_name("CRO")),
      grepl("CAZ",Feature)~str_replace(Feature,"CAZ",ab_name("CAZ")),
      grepl("FEP",Feature)~str_replace(Feature,"FEP",ab_name("FEP")),
      grepl("MEM",Feature)~str_replace(Feature,"MEM",ab_name("MEM")),
      grepl("CIP",Feature)~str_replace(Feature,"CIP",ab_name("CIP")),
      grepl("GEN",Feature)~str_replace(Feature,"GEN",ab_name("GEN")),
      grepl("SXT",Feature)~str_replace(Feature,"SXT",ab_name("SXT")),
      grepl("NIT",Feature)~str_replace(Feature,"NIT",ab_name("NIT")),
      grepl("VAN",Feature)~str_replace(Feature,"VAN",ab_name("VAN")),
      grepl("MTR",Feature)~str_replace(Feature,"MTR",ab_name("MTR")),
      grepl("AZM",Feature)~str_replace(Feature,"AZM",ab_name("AZM")),
      grepl("ERY",Feature)~str_replace(Feature,"ERY",ab_name("ERY")),
      TRUE~Feature
    )
  ) %>% mutate(
    Feature = str_replace(Feature,"_"," & ")
  ) %>% mutate(
    Feature = str_replace(Feature,"\\.","/")
  ) %>% mutate(Feature=case_when(
    grepl("Previous week",Feature)~glue("{Feature %>% str_remove(prev_week)} in the last week"),
    grepl("Ampicillinc",Feature)~str_replace(Feature,"Ampicillinc","AMPc beta-lactamase"),
    grepl("No Over",Feature)~str_replace(Feature,"No Over","Not over"),
    grepl("RDW",Feature)~str_replace(Feature,"RDW","Most recent RDW"),
    grepl("Albumin",Feature)~str_replace(Feature,"Albumin","Most recent Albumin"),
    grepl("PTT",Feature)~str_replace(Feature,"PTT","Most recent prothrombin time"),
    grepl("Alkaline Phosphatase",Feature)~str_replace(Feature,"Alkaline Phosphatase","Most recent Alkaline Phosphatase"),
    grepl("Red blood Cell Count",Feature)~str_replace(Feature,"Red blood Cell Count","Most recent Red blood Cell Count"),
    grepl("Lymphocytes",Feature)~str_replace(Feature,"Lymphocytes","Most recent Lymphocyte Count"),
    grepl("Leukocyte Count",Feature)~str_replace(Feature,"Leukocyte Count","Most recent Leukocyte Count"),
    grepl("Hemoglobin",Feature)~str_replace(Feature,"Hemoglobin","Most recent Haemoglobin"),
    grepl("Platelet Count",Feature)~str_replace(Feature,"Platelet Count","Most recent Platelet Count"),
    grepl("Potassium",Feature)~str_replace(Feature,"Potassium","Most recent Serum Potassium"),
    grepl("Neutrophils",Feature)~str_replace(Feature,"Neutrophils","Most recent Neutrophil Count"),
    grepl("Hematocrit",Feature)~str_replace(Feature,"Hematocrit","Most recent Haematocrit"),
    grepl("Lactate",Feature)~str_replace(Feature,"Lactate","Most recent Lactate"),
    grepl("Monocytes",Feature)~str_replace(Feature,"Monocytes","Most recent Monocyte Count"),
    grepl("Bilirubin",Feature)~str_replace(Feature,"Bilirubin","Most recent Bilirubin"),
    grepl("Sodium",Feature)~str_replace(Feature,"Sodium","Most recent Serum Sodium"),
    grepl("Bicarbonate",Feature)~str_replace(Feature,"Bicarbonate","Most recent Serum Bicarbonate"),
    grepl("Lactate Dehydrogenase",Feature)~str_replace(Feature,"Lactate Dehydrogenase","Most recent Lactate Dehydrogenase"),
    grepl("Creatinine",Feature)~str_replace(Feature,"Creatinine","Most recent Creatinine"),
    grepl("RDW",Feature)~str_replace(Feature,"RDW","Most recent RDW"),
    grepl("No Male",Feature)~str_replace(Feature,"No","Not"),
    TRUE~Feature
  )) %>% mutate(`Shapley value`=round(`Shapley value`,3)) %>% 
  mutate(Model=str_replace_all(Model,"_"," & "))

write_csv(feat_table,"feat_table.csv")
