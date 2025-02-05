###COVERAGE MODEL FINAL TRAINING AND VALIDATION

set.seed(123)

###Read-in chosen model parameters (coverage prediction)
file_names <- glue("final_params_{combined_antimicrobial_map}.csv")
final_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(final_bestparams)) {
  names(final_bestparams[[i]])[3:8] <- namelist
}
names(final_bestparams) <- combined_antimicrobial_map

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
    
    print(glue("Training for {outcome}"))
    
    xgb_model <- xgb.train(
      params = params,
      data = train_matrix,
      nrounds = final_bestparams[[outcome]]$best_nrounds
    )
    
    print(glue("Shapping for {outcome}"))
    
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
    
    print(glue("Second training for {outcome}"))
    
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
    
    print(glue("Bootstrapping for {outcome}"))
    
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
    print(glue("Bootstrapping for {outcome} 2"))
    set.seed(123)
    boot_results <- boot(data = df_ys, statistic = calculate_metrics, R = n_bootstraps)
    boot_results[[1]][is.na(boot_results[[1]])] <- 0
    boot_results[[2]][is.na(boot_results[[2]])] <- 0
    print(glue("Bootstrapping for {outcome} 3"))
    confidence_intervals <- lapply(1:5, function(i) boot.ci(boot_results, type = "perc", index = i))
    print(glue("Bootstrapping for {outcome} 4"))
    names(confidence_intervals) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")
    print(glue("Bootstrapping for {outcome} 5"))
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
  ci_df$Precision_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Precision$percent[4],3)}-{round(confidence_biglist[[i]]$Precision$percent[5],3)})")
  ci_df$Recall_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Recall$percent[4],3)}-{round(confidence_biglist[[i]]$Recall$percent[5],3)})")
  ci_df$Accuracy_CI[i] <- glue(" ({round(confidence_biglist[[i]]$Accuracy$percent[4],3)}-{round(confidence_biglist[[i]]$Accuracy$percent[5],3)})")
  ci_df$F1_CI[i] <- glue(" ({round(confidence_biglist[[i]]$F1_Score$percent[4],3)}-{round(confidence_biglist[[i]]$F1_Score$percent[5],3)})")
  
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