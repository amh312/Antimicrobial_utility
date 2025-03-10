###TOXICITY MODEL FINAL TRAINING AND VALIDATION

set.seed(123)

###Toxicity data partitioning
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
pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
pred_test_class <- relevel(factor(pred_test_class), ref = "1")
actual_test_class <- abxTest[['overall_tox']]
actual_test_class <- relevel(factor(actual_test_class), ref = "1")
roc_result <- roc(actual_test_class, pred_prob_test)
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

calibration_data <- data.frame(
  predicted_prob = pred_prob_test,
  actual = as.numeric(as.character(actual_test_class))
)

calibration_data$bin <- cut(calibration_data$predicted_prob, 
                            breaks = quantile(calibration_data$predicted_prob, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                            include.lowest = TRUE, labels = FALSE)

calibration_summary <- calibration_data %>%
  group_by(bin) %>%
  summarise(mean_pred = mean(predicted_prob),
            observed_freq = mean(actual))

calibration_plot <- ggplot(calibration_summary, aes(x = mean_pred, y = observed_freq)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  theme_minimal() +
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "Mean Predicted Probability", y = "Observed Frequency",
       title = glue("Toxicity\nCalibration Curve")) +
  theme(plot.title = element_text(hjust = 0.5))

calibration_plots[["overall_tox"]] <- calibration_plot

ggsave(filename = glue("overall_tox_calib_plot.pdf"), plot = calibration_plot, width = 10, height = 10)

overall_tox_util_probs <- predict(xgb_model, newdata = micro_matrix)

pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
pred_test_class <- relevel(factor(pred_test_class), ref = "1")
actual_test_class <- abxTest[['overall_tox']]
actual_test_class <- relevel(factor(actual_test_class), ref = "1")

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
set.seed(123)
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

ci_df2 <- read_csv("interim_ci_df.csv")
ci_df <- data.frame(rbind(ci_df2,ci_df))

ci_df$Model <- c(combined_antimicrobial_map,"CDI","overall_tox")

write_csv(ci_df,"final_ci_df.csv")
ci_filter_list <- c(all_singles,"CDI","overall_tox")
ci_df <- ci_df %>% mutate(Model=as.character(Model))
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

##Save all calibration plots to grid PDFs
calibration_plots_main <- c(calibration_plots[1:13],calibration_plots[38:39])
calibration_plots_other <- calibration_plots[14:37]

main_grid_plot <- plot_grid(plotlist = calibration_plots_main, ncol = 3)
other_grid_plot <- plot_grid(plotlist = calibration_plots_other, ncol = 3)

ggsave(glue("calibration_plots_main.pdf"), plot = main_grid_plot, width = 15, height = 30)
ggsave(glue("calibration_plots_other.pdf"), plot = other_grid_plot, width = 15, height = 30)

##Feature contributions dataframe
other_outcome_map <- list("CDI","toxicity")
names(other_outcome_map) <- c("CDI","Toxicity")
ab_filter_list <- probs_df_overall %>% distinct(Antimicrobial) %>% unlist()
full_outcome_map <-combined_antimicrobial_map[names(combined_antimicrobial_map) %in% ab_filter_list]
full_outcome_map <- c(full_outcome_map,other_outcome_map)
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
nt <- "nt$"
s <- "s$"
rx <- "rx$"
d7 <- "^d7"
prev_week <- "Previous week"
space_start <- "^ "

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
    grepl("DIAG_",Feature)~str_replace(Feature,"DIAG_","ICD diagnosis category "),
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
    grepl("highCRP",Feature)~str_replace(Feature,"highCRP","Elevated CRP in the last 24 hours"),
    grepl("ab_name_",Feature)~str_replace(Feature,"ab_name_","Antibiotic(s): "),
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
    grepl("abnormalWCC",Feature)~str_replace(Feature,"abnormalWCC","Abnormal white cell count in last 24 hours"),
    grepl("o2sat",Feature)~str_replace(Feature,"o2sat","Oxygen saturation on admission"),
    grepl("temperat",Feature)~str_replace(Feature,"temperature","Temperature on admission"),
    grepl("arrival_transport",Feature)~str_replace(Feature,"arrival_transport","Means of arrival "),
    grepl("MALE",Feature)~str_replace(Feature,"MALE","Male"),
    grepl("^sbp",Feature)~str_replace(Feature,"sbp","Systolic blood pressure on admission"),
    grepl("^dbp",Feature)~str_replace(Feature,"dbp","Diastolic blood pressure on admission"),
    grepl("^acuity",Feature)~str_replace(Feature,"acuity","Acuity score on admission"),
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
    Feature=case_when(
      grepl("Ampicillins",Feature)~str_replace(Feature,"Ampicillins", "Ampicillin susceptibility in the last year"),
      grepl("Ampicillin/sulbactams",Feature)~str_replace(Feature,"Ampicillin/sulbactams", "Ampicillin/sulbactam susceptibility in the last year"),
      grepl("Piperacillin/tazobactams",Feature)~str_replace(Feature,"Piperacillin/tazobactams", "Piperacillin/tazobactam susceptibility in the last year"),
      grepl("Cefazolins",Feature)~str_replace(Feature,"Cefazolins", "Cefazolin susceptibility in the last year"),
      grepl("Ceftriaxones",Feature)~str_replace(Feature,"Ceftriaxones", "Ceftriaxone susceptibility in the last year"),
      grepl("Ceftazidimes",Feature)~str_replace(Feature,"Ceftazidimes", "Ceftazidime susceptibility in the last year"),
      grepl("Cefepimes",Feature)~str_replace(Feature,"Cefepimes", "Cefepime susceptibility in the last year"),
      grepl("Meropenems",Feature)~str_replace(Feature,"Meropenems", "Meropenem susceptibility in the last year"),
      grepl("Ciprofloxacins",Feature)~str_replace(Feature,"Ciprofloxacins", "Ciprofloxacin susceptibility in the last year"),
      grepl("Gentamicins",Feature)~str_replace(Feature,"Gentamicins", "Gentamicin susceptibility in the last year"),
      grepl("Trimethoprim/sulfamethoxazoles",Feature)~str_replace(Feature,"Trimethoprim/sulfamethoxazoles", "Trimethoprim/sulfamethoxazole susceptibility in the last year"),
      grepl("Nitrofurantoins",Feature)~str_replace(Feature,"Nitrofurantoins", "Nitrofurantoin susceptibility in the last year"),
      grepl("Vancomycins",Feature)~str_replace(Feature,"Vancomycins", "Vancomycin susceptibility in the last year"),
      grepl("Metronidazoles",Feature)~str_replace(Feature,"Metronidazoles", "Metronidazole susceptibility in the last year"),
      grepl("Aztreonams",Feature)~str_replace(Feature,"Aztreonams", "Aztreonam susceptibility in the last year"),
      grepl("Erythromycins",Feature)~str_replace(Feature,"Erythromycins", "Erythromycin susceptibility in the last year"),
      TRUE~Feature
    )
  )  %>% mutate(
    Feature=case_when(
      grepl("Ampicillinnt",Feature)~str_replace(Feature,"Ampicillinnt", "Ampicillin not tested in the last year"),
      grepl("Ampicillin/sulbactamnt",Feature)~str_replace(Feature,"Ampicillin/sulbactamnt", "Ampicillin/sulbactam not tested in the last year"),
      grepl("Piperacillin/tazobactamnt",Feature)~str_replace(Feature,"Piperacillin/tazobactamnt", "Piperacillin/tazobactam not tested in the last year"),
      grepl("Cefazolinnt",Feature)~str_replace(Feature,"Cefazolinnt", "Cefazolin not tested in the last year"),
      grepl("Ceftriaxonent",Feature)~str_replace(Feature,"Ceftriaxonent", "Ceftriaxone not tested in the last year"),
      grepl("Ceftazidiment",Feature)~str_replace(Feature,"Ceftazidiment", "Ceftazidime not tested in the last year"),
      grepl("Cefepiment",Feature)~str_replace(Feature,"Cefepiment", "Cefepime not tested in the last year"),
      grepl("Meropenemnt",Feature)~str_replace(Feature,"Meropenemnt", "Meropenem not tested in the last year"),
      grepl("Ciprofloxacinnt",Feature)~str_replace(Feature,"Ciprofloxacinnt", "Ciprofloxacin not tested in the last year"),
      grepl("Gentamicinnt",Feature)~str_replace(Feature,"Gentamicinnt", "Gentamicin not tested in the last year"),
      grepl("Trimethoprim/sulfamethoxazolent",Feature)~str_replace(Feature,"Trimethoprim/sulfamethoxazolent", "Trimethoprim/sulfamethoxazole not tested in the last year"),
      grepl("Nitrofurantoinnt",Feature)~str_replace(Feature,"Nitrofurantoinnt", "Nitrofurantoin not tested in the last year"),
      grepl("Vancomycinnt",Feature)~str_replace(Feature,"Vancomycinnt", "Vancomycin not tested in the last year"),
      grepl("Metronidazolent",Feature)~str_replace(Feature,"Metronidazolent", "Metronidazole not tested in the last year"),
      grepl("Aztreonamnt",Feature)~str_replace(Feature,"Aztreonamnt", "Aztreonam not tested in the last year"),
      grepl("Erythromycinnt",Feature)~str_replace(Feature,"Erythromycinnt", "Erythromycin not tested in the last year"),
      TRUE~Feature
    )
  ) %>% mutate(
    Feature = str_replace(Feature,"_"," & ")
  ) %>% mutate(
    Feature = str_replace(Feature,"\\.","/")
  ) %>% mutate(Feature=case_when(
    grepl("Previous week",Feature)~glue("{Feature %>% str_remove(prev_week)} in the last week"),
    grepl("AmpicillinC",Feature)~str_replace(Feature,"AmpicillinC","AMPc beta-lactamase"),
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
    grepl("EMERGentamicinCY",Feature)~str_replace(Feature,"entamicin","EN"),
    grepl("abdopain",Feature)~str_replace(Feature,"abdopain","abdominal pain"),
    grepl("chestpain",Feature)~str_replace(Feature,"chestpain","chest pain"),
    TRUE~Feature
  )) %>% mutate(`Shapley value`=round(`Shapley value`,3)) %>% 
  mutate(Feature=str_replace(
    Feature,"in the last year in the last week","in the last week"),
    Feature = str_replace(Feature,"7d",""),
    Feature = str_replace(Feature,"CLI ","Clindamycin "),
    Feature = str_replace(Feature,"AMX","Amoxicillin"),
    Feature = str_replace(Feature,"SURG","surgery"),
    Feature = str_replace(Feature,"TOB","Tobramycin"),
    Feature = str_replace(Feature,"TCY","Tetracycline"),
    Feature = str_replace(Feature,"Provider ID","Presence of provider ID"),
    Feature = str_replace(Feature, "PEN","Benzylpenicillin"),
    Feature = str_replace(Feature, "LNZ","Linezolid"),
    Feature = str_replace(Feature, "LVX","Levofloxacin"),
    Feature=str_remove(Feature,"^ "),
    Feature = str_replace(Feature, "  "," ")) %>% 
  mutate(
    Feature = str_replace(Feature,"CLIs ","Clindamycin susceptibility in the last year"),
    Feature = str_replace(Feature,"Levofloxacins","Levofloxacin susceptibility in the last year"),
    Feature = str_replace(Feature,"Tetracyclines","Tetracycline susceptibility in the last year"),
    Feature = str_replace(Feature,"Benzylpenicillins","Penicillin susceptibility in the last year"),
    Feature = str_replace(Feature,"AMKs","Amikacin susceptibility in the last year"),
    Feature = str_replace(Feature,"AMK","Amikacin")
  ) %>% 
  mutate(Model=str_replace_all(Model,"_"," & "))

write_csv(feat_table,"feat_table.csv")

feat_abs <- c(c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN") %>% ab_name() %>% str_replace_all("/","-"),
"CDI","Toxicity")
feat_table_singles <- feat_table %>% filter(Model %in% feat_abs) %>% 
  mutate(Model = factor(Model, levels = unique(Model))) %>%
  group_by(Model) %>% arrange(desc(`Shapley value`),.by_group = T) %>% ungroup() %>% 
  mutate(`Shapley value`=case_when(`Shapley value`<0.001 ~ as.character("<0.001"),
                                   TRUE~as.character(`Shapley value`)))
write_csv(feat_table_singles,"feat_table_singles.csv")
