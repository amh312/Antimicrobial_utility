###CDI MODEL FINAL TRAINING AND VALIDATION

set.seed(123)

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
confidence_biglist <- list()

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
set.seed(123)
boot_results <- boot(data = df_ys, statistic = calculate_metrics, R = n_bootstraps)

confidence_intervals <- lapply(1:5, function(i) boot.ci(boot_results, type = "perc", index = i))

names(confidence_intervals) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")

confidence_biglist[['CDI']] <- confidence_intervals