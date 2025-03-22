###COVERAGE MODEL FINAL TRAINING AND VALIDATION

set.seed(123)

##Functions

###Reading in hyperparameters back to list from csvs
hypparamreader <- function(namestring,namesmap) {
  
  #write file nanmes
  file_names <- glue("{namestring}{namesmap}.csv")
  
  #apply read_csv across name list
  hypparlist <- lapply(file_names, function(file) {
    read_csv(file)
  })
  
  #add names
  namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
                "best_nrounds")
  
  #iteratively add names to sublists
  for (i in 1:length(hypparlist)) {
    names(hypparlist[[i]])[3:8] <- namelist
  }
  
  #add model names to list top level
  names(hypparlist) <- namesmap
  
  hypparlist
  
}

###Generate combinations of antibiotics in a list
ab_combiner <- function(abx_list,newname,combname) {
  
  #get list of all combinations
  comblist <- combn(abx_list, 2, function(x) paste(x, collapse = "_"))
  
  #combine singles with combinations
  all_list <- c(abx_list,comblist)
  
  #assign to global environment
  assign(newname,comblist,envir = .GlobalEnv)
  assign(combname,all_list,envir = .GlobalEnv)
  
}

###Generate combinations of antibiotics in a mapped list
abmap_combiner <- function(ab_map) {
  
  #get all 2-ab combinations in the map
  abcombos_map <- combn(names(ab_map), 2, simplify = FALSE)
  
  #bind combos on to singles
  c(ab_map,
    setNames(
      
      #bind on combo elements themselves
      lapply(abcombos_map, function(x) paste(ab_map[x], collapse = "_")),
      
      #bind on combo element names
      sapply(abcombos_map, function(x) paste(x, collapse = "_"))
    )
  )
  
}

###Make filtered, full and short antimicrobial maps
abcombo_variants <- function(abmap,fullname,mainname,shortname,abreflist) {
  
  #get vector of all antibiotics as unlabelled list
  fullmap <- abmap %>% unlist()
  names(fullmap) <- NULL
  
  #get combined ab map just for abs in prescription df
  abmap2 <- abmap[names(abmap) %in% abreflist]
  
  #get unlabelled vector of the shortened combo list
  shortmap <- abmap2 %>% unlist()
  names(shortmap) <- NULL
  
  #assign new objects to glob env
  assign(fullname,fullmap,envir = .GlobalEnv)
  assign(mainname,abmap2,envir = .GlobalEnv)
  assign(shortname,shortmap,envir = .GlobalEnv)
  
}

###Replace antimicrobial short name with long name including for combinations
abcombo_replace <- function(abtarget, map) {
  
  #switch names of combined ab list with the list elements themselves
  flip_abmap <- setNames(names(map), map)
  
  abtarget %>%
    
    #convert column to character format
    as.character() %>%
    
    #if target value is in names of flipped map, return the element
    sapply(function(x) if (x %in% names(flip_abmap)) flip_abmap[[x]] else x)
  
}

###Replace antimicrobial short name with long name including for combinations
abcombo_reverse <- function(abtarget, map) {
  
  abtarget %>%
    
    #convert to character
    as.character() %>%
    
    #if target ab is in the map elements, return that name
    sapply(function(x) if (x %in% map) names(map)[map == x] else x)
  
}

###Train-test splitting
TTsplitter <- function(dataset,outc,trainprop,chosnametrain,chosnametest){
  
  #get partition index based on specified proportion
  trainindex <- createDataPartition(dataset[[outc]], p = trainprop, list = FALSE, times = 1)
  
  #index training dataset
  urdftrain <- dataset[trainindex, ]
  
  #index testing dataset
  urdftest <- dataset[-trainindex, ]
  
  #assign to objects with 'Train' and 'Test' suffixes replacing '_combined' suffix
  assign(chosnametrain,
         urdftrain,.GlobalEnv)
  assign(chosnametest,
         urdftest,.GlobalEnv)
  
}

###Line up features between dataframes
lineup_features <- function(microsim_df,pred_df,train_df){
  
  predictor_columns <- colnames(pred_df)
  selected_columns <- intersect(predictor_columns, colnames(train_df))
  missing_cols <- setdiff(selected_columns, colnames(microsim_df))
  microsim_df[missing_cols] <- 0
  
  microsim_df
  
}

###Generate XGBoost matrix from model train or test dataset
model_matrixmaker <- function(dataset,predictors,outc) {
  
  #get predictor feature names
  predictorcols <- colnames(predictors)
  
  #ensure that all predictor column names are in specified data
  selected_columns <- intersect(predictorcols, colnames(dataset))
  
  #make xgboost matrix
  xgb.DMatrix(data = as.matrix(dataset %>% select(all_of(selected_columns))), 
              label = dataset[[outc]])
  
  
}

###Model training
xg_training <- function(trainmat,lr,mtd,mcw,ss,csb,bnr,outc){
  
  urparams <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    eta = lr,
    max_depth = mtd,
    min_child_weight = mcw,
    subsample = ss,
    colsample_bytree = csb
  )
  
  print(glue("Training for {outc}"))
  
  xgb.train(
    params = urparams,
    data = trainmat,
    nrounds = bnr
  )
  
}

###Feature importance calculation
shapper <- function(trainmat,model,outc) {
  
  print(glue("Checking feature importances for {outc}"))
  
  ur_shapvals <- predict(xgb_urinemodel, newdata = urtrain_matrix, predcontrib = TRUE)
  shap_df <- as.data.frame(ur_shapvals)
  shap_df <- shap_df[, -ncol(shap_df)]
  shap_ursmry <- data.frame(
    Feature = colnames(shap_df),
    meanshap = colMeans(shap_df)
  )
  
  shap_ursmry <- shap_ursmry %>% filter(meanshap!=0)
  shap_ursmry <- shap_ursmry[order(-shap_ursmry$meanshap), ]
  
  shap_ursmry
  
}

###Feature selection using SHAP
mat_feat_selector <- function(dataset,shapsum,outc) {
  
  xgb.DMatrix(data = as.matrix(dataset %>% select(shapsum %>% pull(Feature))), 
              label = dataset[[outc]])
  
}

##Read-in

train_abx <- read_csv("train_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
urines5_combined <- read_csv("urines5_combined.csv")
ur_xg_combined <- read_csv("ur_xg_combined.csv")

##Reference lists

###Antimicrobial mapping lists
antimicrobial_map <- c("Ampicillin" = "AMP","Ampicillin-sulbactam" = "SAM",
                       "Piperacillin-tazobactam" = "TZP","Cefazolin" = "CZO","Ceftriaxone" = "CRO",
                       "Ceftazidime" = "CAZ","Cefepime" = "FEP","Meropenem" = "MEM","Ciprofloxacin" = "CIP",
                       "Gentamicin" = "GEN","Trimethoprim-sulfamethoxazole" = "SXT","Nitrofurantoin" = "NIT",
                       "Vancomycin" = "VAN")
all_singles <- antimicrobial_map %>% unname()
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM","GEN","SXT")
oral_singles <- c("AMP","SAM","CIP","SXT","NIT")
access_singles <- c("AMP","SAM","GEN","SXT","NIT","CZO")
watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
singlelists <- c("all_singles","iv_singles","oral_singles","access_singles","watch_singles")
comblists <- str_replace(singlelists,"singles","combos")
alllists <- c("all_abs","all_ivs","all_orals","all_access","all_watch")
for (i in seq_along(singlelists)) {ab_combiner(get(singlelists[i]),comblists[i],alllists[i])}
combined_antimicrobial_map <- abmap_combiner(antimicrobial_map)
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
abcombo_variants(combined_antimicrobial_map,
                 "fullmap","combined_antimicrobial_map","shortmap",abx_in_train)

###Model testing and microsimulation reference dataframes
test_probs_df <- data.frame(matrix(nrow=floor(nrow(urines5_combined)*0.2),ncol=0))
micro_probs_df <- data.frame(matrix(nrow=nrow(ur_xg_combined),ncol=0))

###Blank dataframes
aucs <- data.frame(matrix(nrow=1,ncol=0))
shap_ur_summary_tables <- list()
metrics_list <- list()
roc_plots <- list()
calibration_plots <- list()
confidence_biglist <- list()

##Model training and validation

###Read in hyperparameters
final_bestparams <- hypparamreader("final_params_",combined_antimicrobial_map)

###Iterate over antimicrobial agents
for (outcome in colnames(urines5_outcomes)) {
  
  ###Check no NAs
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    ###Set seed
    set.seed(123)
    
    ###Split df to train and test
    urines5_combined %>% TTsplitter(outcome,0.8,"urines5Train","urines5Test")
    
    ###Ensure features line up between dataframes
    ur_xg_combined <- ur_xg_combined %>%
      lineup_features(urines5_predictors,urines5Train)
    
    ###Make xgboost training matrices
    urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,'AMP')
    urtest_matrix <- urines5Test %>% model_matrixmaker(urines5_predictors,'AMP')
    urmicro_matrix <- ur_xg_combined %>% model_matrixmaker(urines5_predictors,'AMP')
    
    ###First training
    xgb_urinemodel <- urtrain_matrix %>%
      xg_training(lr=final_bestparams[[outcome]]$eta,
                  mtd=final_bestparams[[outcome]]$max_depth,
                  mcw=final_bestparams[[outcome]]$min_child_weight,
                  ss=final_bestparams[[outcome]]$subsample,
                  csb=final_bestparams[[outcome]]$colsample_bytree,
                  bnr=final_bestparams[[outcome]]$best_nrounds,
                  outc=outcome)
    
    ###Feature importances
    shap_ur_summary <- urtrain_matrix %>% shapper(xgb_urinemodel,outcome)
    
    ###Feature selection
    urtrain_matrix <- urines5Train %>% mat_feat_selector(shap_ur_summary,outcome)
    urtest_matrix <- urines5Test %>% mat_feat_selector(shap_ur_summary,outcome)
    urmicro_matrix <- ur_xg_combined %>% mat_feat_selector(shap_ur_summary,outcome)
    
    ###Second training
    xgb_urinemodel <- urtrain_matrix %>%
      xg_training(lr=final_bestparams[[outcome]]$eta,
                  mtd=final_bestparams[[outcome]]$max_depth,
                  mcw=final_bestparams[[outcome]]$min_child_weight,
                  ss=final_bestparams[[outcome]]$subsample,
                  csb=final_bestparams[[outcome]]$colsample_bytree,
                  bnr=final_bestparams[[outcome]]$best_nrounds,
                  outc=outcome)
    
    
    
    
    pred_prob_test <- predict(xgb_urinemodel, newdata = test_matrix)
    pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
    pred_test_class <- relevel(factor(pred_test_class), ref = "1")
    actual_test_class <- urines5Test[[outcome]]
    actual_test_class <- relevel(factor(actual_test_class), ref = "1")
    roc_result <- roc(actual_test_class, pred_prob_test)
    auc_value <- auc(roc_result)
    print(paste("AUC-ROC:", auc_value))
    aucs[[outcome]] <- auc_value
    roc_plot <- ggroc(roc_result) + 
      ggtitle(glue("{abcombo_replace(outcome,combined_antimicrobial_map)}\nROC Curve"))+
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
           title = glue("{abcombo_replace(outcome,combined_antimicrobial_map)}\nCalibration Curve")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    calibration_plots[[outcome]] <- calibration_plot
    print(calibration_plot)
    
    ggsave(filename = glue("{outcome}_calib_plot.pdf"), plot = calibration_plot, width = 10, height = 10)
    
    pred_prob_micro <- predict(xgb_urinemodel, newdata = micro_matrix)
    
    test_probs_df[[outcome]] <- pred_prob_test
    micro_probs_df[[outcome]] <- pred_prob_micro
    
    
    
    shap_ur_summary_tables[[outcome]] <- shap_ur_summary
    
    pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
    pred_test_class <- relevel(factor(pred_test_class), ref = "1")
    actual_test_class <- urines5Test[[outcome]]
    actual_test_class <- relevel(factor(actual_test_class), ref = "1")
    
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
  print(paste("AUC-ROC:", auc_value))
}


for (i in 1:length(shap_ur_summary_tables)) {
  
  shappy <- data.frame(shap_ur_summary_tables[i]) %>% mutate(across(2, ~ round(., 3))) %>% filter(if_any(2, ~ . != 0))
  
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
  
  ci_df$AUROC_CI[i] <- ifelse(
    !is.null(confidence_biglist[[i]]$AUC$percent[4])&!is.null(confidence_biglist[[i]]$AUC$percent[5]),
    glue(" ({round(confidence_biglist[[i]]$AUC$percent[4],3)}-{round(confidence_biglist[[i]]$AUC$percent[5],3)})"),
    NA)
  ci_df$Precision_CI[i] <- ifelse(
    !is.null(confidence_biglist[[i]]$Precision$percent[4])&!is_null(confidence_biglist[[i]]$Precision$percent[5]),
    glue(" ({round(confidence_biglist[[i]]$Precision$percent[4],3)}-{round(confidence_biglist[[i]]$Precision$percent[5],3)})"),
  NA)
  ci_df$Recall_CI[i] <- ifelse(
    !is.null(confidence_biglist[[i]]$Recall$percent[4])&!is_null(confidence_biglist[[i]]$Recall$percent[5]),
    glue(" ({round(confidence_biglist[[i]]$Recall$percent[4],3)}-{round(confidence_biglist[[i]]$Recall$percent[5],3)})"),
  NA)
  ci_df$Accuracy_CI[i] <- ifelse(
    !is.null(confidence_biglist[[i]]$Accuracy$percent[4])&!is_null(confidence_biglist[[i]]$Accuracy$percent[5]),
    glue(" ({round(confidence_biglist[[i]]$Accuracy$percent[4],3)}-{round(confidence_biglist[[i]]$Accuracy$percent[5],3)})"),
  NA)
  ci_df$F1_CI[i] <- ifelse(
    !is.null(confidence_biglist[[i]]$F1_Score$percent[4])&!is_null(confidence_biglist[[i]]$F1_Score$percent[5]),
    glue(" ({round(confidence_biglist[[i]]$F1_Score$percent[4],3)}-{round(confidence_biglist[[i]]$F1_Score$percent[5],3)})"),
  )
}

write_csv(ci_df,"interim_ci_df.csv")

###Microsimulation dataframe
micro_probs_df$micro_specimen_id <- ur_xg$micro_specimen_id
micro_probs_df$subject_id <- ur_xg$subject_id
probs_df_overall <- micro_probs_df %>% melt(id.vars=c("micro_specimen_id","subject_id")) %>% 
  mutate(variable=abcombo_replace(variable, combined_antimicrobial_map)) %>% 
  rename(S="value",Antimicrobial="variable") %>% mutate(I=NA,R=1-S,NT=NA)
probs_df_overall <- probs_df_overall %>% relocate(micro_specimen_id,.after = "NT") %>% 
  relocate(subject_id,.after="micro_specimen_id") %>% 
  relocate(I,.before="S") %>% mutate(id_no=1:nrow(probs_df_overall)) %>% 
  relocate(id_no,.before="Antimicrobial")
write_csv(probs_df_overall,"probs_df_overall.csv")
