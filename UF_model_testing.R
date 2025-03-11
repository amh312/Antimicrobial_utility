#MODEL ADDITIONAL TESTING

set.seed(123)

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Time sens plot
stability_plot <- function(df,metric,perf_metric) {
  
  metric <- enquo(metric)
  
  df <- df %>% 
    filter(!is.na(!!metric)) %>% 
    mutate(Model=case_when(Model=="overall_tox"~"Toxicity",
                                      Model=="CDI"~"CDI",
                                      TRUE~ab_name(Model)))
  df <- df %>% rename(`Training dataset size`="Training_size")
  model_levels <- ab_name(all_singles) %>% append(c("Vancomycin","CDI","Toxicity")) %>% rev()
  df$Model <- factor(df$Model,levels=model_levels)
  df$`Training dataset size` <- as.character(df$`Training dataset size`)
  dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=`Training dataset size`,color=`Training dataset size`)) +
    geom_point()+
    theme_minimal()+
    ggtitle(glue("{perf_metric} after 50 XGBoost training rounds using different\ntraining dataset proportions"))+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    xlim(0,1)
  
  ggsave(glue("stability_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(dfplot)
  
}

###Fairness plot
protchar_plot <- function(df,prot_char,metric,perf_metric,title_bit) {
  
  metric <- enquo(metric)
  
  df <- df %>% filter(grepl(prot_char,Category)) %>% 
    filter(!is.na(!!metric)) %>%
    mutate(Model=case_when(Model=="overall_tox"~"Toxicity",
                           Model=="CDI"~"CDI",
                           TRUE~ab_name(Model)))
  model_levels <- ab_name(all_singles) %>% append(c("Vancomycin","CDI","Toxicity")) %>% rev()
  df$Model <- factor(df$Model,levels=model_levels)
  df$Characteristic <- factor(df$Characteristic,
                              levels=df %>% filter(grepl(prot_char,Category)) %>% 
                                distinct(Characteristic) %>% unlist())
  dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=Characteristic,color=Characteristic)) +
    geom_point()+
    theme_minimal()+
    ggtitle(glue("{perf_metric} after 50 XGBoost training rounds for different\n{title_bit}"))+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    xlim(0,1)
  
  ggsave(glue("protchar_{prot_char}_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(dfplot)
  
}

###Time plot
timesens_plot <- function(df,tr_yr,metric,perf_metric) {
  
  metric <- enquo(metric)
  
  df <- df %>% filter(grepl(tr_yr,Train_year)) %>% 
    filter(!is.na(!!metric)) %>%
    mutate(Model=case_when(Model=="overall_tox"~"Toxicity",
                           Model=="CDI"~"CDI",
                           TRUE~ab_name(Model)))
  model_levels <- ab_name(all_singles) %>% append(c("Vancomycin","CDI","Toxicity")) %>% rev()
  df$Model <- factor(df$Model,levels=model_levels)
  df$Test_year <- factor(df$Test_year,
                         levels=df %>% filter(grepl(tr_yr,Train_year)) %>% 
                           distinct(Test_year) %>% unlist())
  df <- df %>% rename(`Test year range`="Test_year")
  dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=`Test year range`,color=`Test year range`)) +
    geom_point()+
    theme_minimal()+
    ggtitle(glue("{perf_metric} after 50 XGBoost training rounds for different\ntesting timeframes when trained on {tr_yr}"))+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    xlim(0,1)
  
  ggsave(glue("timesens_{tr_yr}_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(dfplot)
  
}

##Read-in and mapping lists

###Read-in
train_abx <- read_csv("train_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
urines_ref <- read_csv("urines_ref.csv")
pats <- read_csv("patients.csv")
train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")
util_probs_df <- read_csv("probs_df_overall.csv")
hadm <- read_csv("admissions.csv")

###Make AUC list
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
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
metrics_ablist <- c(combined_antimicrobial_map,"CDI","toxicity")

###Performance metrics summary
auc_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
auc_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  auc_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("AUC"))
  
}
colnames(auc_df) <- c("Antimicrobial","AUC")
auc_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
auc_df[38,1] <- "CDI"
auc_df[39,1] <- "toxicity"
auc_df$Antimicrobial <- as.character(auc_df$Antimicrobial)
write_csv(auc_df,"auc_df.csv")
recall_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
recall_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  recall_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Recall"))
  
}
colnames(recall_df) <- c("Antimicrobial","recall")
recall_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
recall_df[38,1] <- "CDI"
recall_df[39,1] <- "toxicity"
recall_df$Antimicrobial <- as.character(recall_df$Antimicrobial)
write_csv(recall_df,"recall_df.csv")
precision_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
precision_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  precision_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Precision"))
  
}
colnames(precision_df) <- c("Antimicrobial","precision")
precision_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
precision_df[38,1] <- "CDI"
precision_df[39,1] <- "toxicity"
precision_df$Antimicrobial <- as.character(precision_df$Antimicrobial)
write_csv(precision_df,"precision_df.csv")
accuracy_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
accuracy_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  accuracy_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Accuracy"))
  
}
colnames(accuracy_df) <- c("Antimicrobial","accuracy")
accuracy_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
accuracy_df[38,1] <- "CDI"
accuracy_df[39,1] <- "toxicity"
accuracy_df$Antimicrobial <- as.character(accuracy_df$Antimicrobial)
write_csv(accuracy_df,"accuracy_df.csv")
f1_score_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
f1_score_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  f1_score_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("F1"))
  
}
colnames(f1_score_df) <- c("Antimicrobial","f1_score")
f1_score_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
f1_score_df[38,1] <- "CDI"
f1_score_df[39,1] <- "toxicity"
f1_score_df$Antimicrobial <- as.character(f1_score_df$Antimicrobial)
write_csv(f1_score_df,"f1_score_df.csv")

metrics_manus_table <- auc_df %>% left_join(precision_df,by="Antimicrobial") %>% 
  left_join(recall_df,by="Antimicrobial") %>%
  left_join(f1_score_df,by="Antimicrobial") %>%
  left_join(accuracy_df,by="Antimicrobial") %>% 
  rename(Model="Antimicrobial",AUROC="AUC",Precision="precision",
         Recall="recall",`F1 score`="f1_score", Accuracy="accuracy") %>% 
  mutate(across(AUROC:`Accuracy`, ~ round(.x,3)),
         Model=str_replace_all(Model,"_"," & "))
metrics_singles_table <- metrics_manus_table %>% dplyr::slice(
  c(1:13,38:39)
) %>% mutate(Model=case_when(Model=="toxicity"~"Toxicity",TRUE~Model))
metrics_combos_table <- metrics_manus_table %>% dplyr::slice(
  -c(1:13,38:39)
) %>% mutate(Model=case_when(Model=="toxicity"~"Toxicity",TRUE~Model))

ci_singles_table <- read_csv("ci_singles_table.csv")
ci_combos_table <- read_csv("ci_combos_table.csv")


metrics_singles_table <- metrics_singles_table %>% left_join(ci_singles_table)
metrics_combos_table <- metrics_combos_table %>% left_join(ci_combos_table)

metrics_singles_table <- metrics_singles_table %>% mutate(
  AUROC = case_when(!is.na(AUROC_CI)&!is.na(AUROC)~glue("{AUROC} {AUROC_CI}"),
                    TRUE~glue("{AUROC} (NA)")),
  Precision = case_when(!is.na(Precision_CI)&!is.na(Precision)~glue("{Precision} {Precision_CI}"),
                    TRUE~glue("{Precision} (NA)")),
  Recall = case_when(!is.na(Recall_CI)&!is.na(Recall)~glue("{Recall} {Recall_CI}"),
                    TRUE~glue("{Recall} (NA)")),
  Accuracy = case_when(!is.na(Accuracy_CI)&!is.na(Accuracy)~glue("{Accuracy} {Accuracy_CI}"),
                    TRUE~glue("{Accuracy} (NA)")),
  `F1 score` = case_when(!is.na(F1_CI)&!is.na(`F1 score`)~glue("{`F1 score`} {F1_CI}"),
                    TRUE~glue("{`F1 score`} (NA)"))) %>% select(Model:Accuracy)

metrics_combos_table <- metrics_combos_table %>% mutate(
  AUROC = case_when(!is.na(AUROC_CI)&AUROC_CI!="NA"~glue("{AUROC} {AUROC_CI}"),
                    TRUE~glue("{AUROC} (NA)")),
  Precision = case_when(!is.na(Precision_CI)&Precision_CI!="NA"~glue("{Precision} {Precision_CI}"),
                        TRUE~glue("{Precision} (NA)")),
  Recall = case_when(!is.na(Recall_CI)&Recall_CI!="NA"~glue("{Recall} {Recall_CI}"),
                     TRUE~glue("{Recall} (NA)")),
  Accuracy = case_when(!is.na(Accuracy_CI)&Accuracy_CI!="NA"~glue("{Accuracy} {Accuracy_CI}"),
                       TRUE~glue("{Accuracy} (NA)")),
  `F1 score` = case_when(!is.na(F1_CI)&F1_CI!="NA"~glue("{`F1 score`} {F1_CI}"),
                         TRUE~glue("{`F1 score`} (NA)"))) %>% select(Model:Accuracy)

write_csv(metrics_singles_table,"metrics_singles_table.csv")
write_csv(metrics_combos_table,"metrics_combos_table.csv")

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
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
fullmap <- combined_antimicrobial_map %>% unlist()
names(fullmap) <- NULL
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
shortmap <- combined_antimicrobial_map %>% unlist()
names(shortmap) <- NULL

###Final preprocessing and dataset train/test splits - urines dataset
urines5 <- urines5 %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                       TRUE~marital_status))
urines5_outcomes <- urines5 %>%
  select(all_of(shortmap))
ur_xg_outcomes <- ur_xg %>%
  select(all_of(shortmap))
urines5_outcomes <- urines5_outcomes %>%
  mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                 ifelse(. == "S" | . == "I", 1, NA))))
urines5_predictors <- urines5 %>% select(!all_of(fullmap))
set.seed(123)
dummies <- dummyVars(" ~ .", data = urines5_predictors)
urines5_predictors <- predict(dummies, newdata = urines5_predictors)
urines5_combined <- as.data.frame(cbind(urines5_outcomes, urines5_predictors))

###Read-in chosen model parameters (resistance prediction)
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

###Model stability analysis (resistance prediction)

###Iterate training and testing over training dataset sizes
tr_size_seq <- c(0.02,0.06,0.10,0.14)
metrics_biglist <- list()
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    metrics_medlist <- list()
    
    for(siz in seq_along(tr_size_seq)) {
      
    metrics_litlist <- list()
    
    for (seedpick in seq(1,6)) {
      
      metrics_list <- list()
      
      iterrun <- glue("{outcome} interation {seedpick} for {tr_size_seq[siz]}")
      print(iterrun)
    
    set.seed(seedpick)
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = tr_size_seq[siz], list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    
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
    
    pred_prob_test <- predict(xgb_model, newdata = test_matrix)
    roc_result <- roc(urines5Test[[outcome]], pred_prob_test)
    auc_value <- auc(roc_result)
    print(paste("AUC-ROC:", auc_value))
    pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
    plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
    dev.off()
    pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
    pred_test_class <- factor(pred_test_class,levels=c(1,0))
    actual_test_class <- urines5Test[[outcome]]
    actual_test_class <- factor(actual_test_class,levels=c(1,0))
    
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
    
    metrics_litlist[[seedpick]] <- metrics_list
    
    }
    
    metrics_medlist[[siz]] <- metrics_litlist
    names(metrics_medlist[[siz]]) <- tr_size_seq[siz]
    
  }
}
  
  metrics_biglist[[outcome]] <- metrics_medlist

}

###Write stability metrics to CSV
metrics_df <- data.frame(matrix(nrow=0,ncol=7))
for (key in 1:length(names(metrics_biglist))) {
  for (i in 1:length(metrics_biglist[[key]])) {
    
    for (j in 1:length(metrics_biglist[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist[[key]][[i]][[j]])) {
   
    print(metrics_biglist[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist)[key],
          Training_size = names(metrics_biglist[[key]][[i]][1]),
          AUC = metrics_biglist[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist[[key]][[i]][[j]][[k]]$F1_Score)
        
        colnames(metrics_df) <- colnames(results)
        
        metrics_df <- data.frame(
          rbind(
            metrics_df,results
          )
        )
     
  }
  }
  }
}
rownames(metrics_df) <- NULL
write_csv(metrics_df,"stability_metrics.csv")

##Model fairness analysis

###Protected characteristics excluding age
ur_prot <- urines5_combined %>% mutate(
  language_NONENG = case_when(languageENGLISH!=1~1,TRUE~0),
  race_ASIAN = rowSums(select(., contains("raceASIAN"))) > 0,
  race_BLACK = rowSums(select(., contains("raceBLACK"))) > 0,
  race_HISPANIC = rowSums(select(., contains("raceHISPANIC"))) > 0,
  race_OTHER = rowSums(select(., matches("(raceOTHER|racePORTUGUESE|raceSOUTH|raceNATIVE|raceAMERICAN|raceMULTIPLE)"))) > 0,
  race_WHITE = rowSums(select(., contains("raceWHITE"))) > 0,
  marital_status_MARRIED = rowSums(select(., matches("(marital_statusMARRIED)"))) > 0,
  marital_status_NONMARRIED = case_when(marital_statusMARRIED!=1~1,TRUE~0)
  )
protchar_index <- which(grepl("(^marital_status_|^language|^MALE|^race_)",colnames(ur_prot))&
                          !grepl("(UNKNOWN|Other|\\?)",colnames(ur_prot)))
protchars <- colnames(ur_prot)[protchar_index]
metrics_biglist2 <- list()
for (outcome in colnames(urines5_outcomes)[1:13]) {
    
    if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
      
      metrics_medlist2 <- list()
        
        for (seedpick in seq(1,6)) {
          
          set.seed(seedpick)
          
          trainIndex <- createDataPartition(urines5_combined[[outcome]], p = 0.8, list = FALSE, times = 1)
          urines5Train <- urines5_combined[trainIndex, ]
          urines5Test <- urines5_combined[-trainIndex, ]
          urines5Test <- urines5Test %>% mutate(
            language_NONENG = case_when(languageENGLISH!=1~1,TRUE~0),
            race_ASIAN = rowSums(select(., contains("raceASIAN"))) > 0,
            race_BLACK = rowSums(select(., contains("raceBLACK"))) > 0,
            race_HISPANIC = rowSums(select(., contains("raceHISPANIC"))) > 0,
            race_OTHER = rowSums(select(., matches("(raceOTHER|racePORTUGUESE|raceSOUTH|raceNATIVE|raceAMERICAN|raceMULTIPLE)"))) > 0,
            race_WHITE = rowSums(select(., contains("raceWHITE"))) > 0,
            marital_status_MARRIED = rowSums(select(., matches("(marital_statusMARRIED)"))) > 0,
            marital_status_NONMARRIED = case_when(marital_statusMARRIED!=1~1,TRUE~0)
          )
          
          predictor_columns <- colnames(urines5_predictors)
          selected_columns <- intersect(predictor_columns, colnames(urines5Train))
          train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                      label = urines5Train[[outcome]])
          test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                                     label = urines5Test[[outcome]])
          
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
          
          metrics_litlist2 <- list()
          
          for(protchar in seq_along(protchar_index)) {
            
            iterrun <- glue("{outcome} iteration {seedpick} for {protchars[protchar]}")
            print(iterrun)
            
            metrics_list2 <- list()
          
          urines5Test2 <- urines5Test[urines5Test[protchar_index[protchar]]==1,]
          
          if (length(urines5Test2[[outcome]] %>% unique()) > 1) {
          
            
            test_matrix2 <- xgb.DMatrix(data = as.matrix(urines5Test2 %>% select(all_of(selected_columns))), 
                                       label = urines5Test2[[outcome]])
            
          pred_prob_test <- predict(xgb_model, newdata = test_matrix2)
          roc_result <- roc(urines5Test2[[outcome]], pred_prob_test)
          auc_value <- auc(roc_result)
          print(paste("AUC-ROC:", auc_value))
          pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
          plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
          dev.off()
          
          pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
          pred_test_class <- factor(pred_test_class,levels=c(1,0))
          actual_test_class <- urines5Test2[[outcome]]
          actual_test_class <- factor(actual_test_class,levels=c(1,0))
          
          confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
          accuracy <- confusion$overall['Accuracy']
          precision <- confusion$byClass['Precision']
          recall <- confusion$byClass['Recall']
          f1_score <- 2 * (precision * recall) / (precision + recall)
          
          metrics_list2[[outcome]] <- list(
            Characteristic = protchars[[protchar]],
            AUC = auc_value,
            Accuracy = accuracy,
            Precision = precision,
            Recall = recall,
            F1_Score = f1_score,
            Test_support = nrow(urines5Test2)
          )
          
          metrics_litlist2[[protchar]] <- metrics_list2
          
          } else {
            
            metrics_litlist2[[protchar]] <- list(
              Characteristic = protchars[[protchar]],
              AUC = NA,
              Accuracy = NA,
              Precision = NA,
              Recall = NA,
              F1_Score = NA,
              Test_support = nrow(urines5Test2)
            )
            
          }
          
        }
        
        metrics_medlist2[[seedpick]] <- metrics_litlist2
      }
    }
      
    metrics_biglist2[[outcome]] <- metrics_medlist2
    
}
metrics_df2 <- data.frame(matrix(nrow=0,ncol=9))
for (key in 1:length(names(metrics_biglist2))) {
  for (i in 1:length(metrics_biglist2[[key]])) {
    
    for (j in 1:length(metrics_biglist2[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist2[[key]][[i]][[j]])) {
        
        print(metrics_biglist2[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist2)[[key]],
          Iteration = i,
          Characteristic = metrics_biglist2[[key]][[i]][[j]][[k]]$Characteristic,
          AUC = metrics_biglist2[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist2[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist2[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist2[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist2[[key]][[i]][[j]][[k]]$F1_Score,
          Test_support = metrics_biglist2[[key]][[i]][[j]][[k]]$Test_support)
        
        colnames(metrics_df2) <- colnames(results)
        
        metrics_df2 <- data.frame(
          rbind(
            metrics_df2,results
          )
        )
        
      }
    }
  }
}
rownames(metrics_df2) <- NULL
write_csv(metrics_df2,"fairness_metrics1.csv")

###Age
ages <- urines5_combined %>% distinct(standard_age) %>% unlist() %>% sort()
metrics_biglist3 <- list()
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    metrics_medlist3 <- list()
      
      for (seedpick in seq(1,6)) {
        
        set.seed(seedpick)
        
          trainIndex <- createDataPartition(urines5_combined[[outcome]], p = 0.8, list = FALSE, times = 1)
          urines5Train <- urines5_combined[trainIndex, ]
          urines5Test <- urines5_combined[-trainIndex, ]
          
          predictor_columns <- colnames(urines5_predictors)
          selected_columns <- intersect(predictor_columns, colnames(urines5Train))
          train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                      label = urines5Train[[outcome]])
          test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                                     label = urines5Test[[outcome]])
          
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
          
          metrics_litlist3 <- list()
          
          for(age in seq_along(ages)) {
            
            iterrun <- glue("{outcome} iteration {seedpick} for age {ages[age]}")
            print(iterrun)
            
            metrics_list3 <- list()
            
            if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
              
              urines5Test2 <- urines5Test %>% filter(standard_age==ages[age])
              test_matrix2 <- xgb.DMatrix(data = as.matrix(urines5Test2 %>% select(all_of(selected_columns))), 
                                         label = urines5Test2[[outcome]])
          
          pred_prob_test <- predict(xgb_model, newdata = test_matrix2)
          roc_result <- roc(urines5Test2[[outcome]], pred_prob_test)
          auc_value <- auc(roc_result)
          print(paste("AUC-ROC:", auc_value))
          pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
          plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
          dev.off()
          
          pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
          pred_test_class <- factor(pred_test_class,levels=c(1,0))
          actual_test_class <- urines5Test2[[outcome]]
          actual_test_class <- factor(actual_test_class,levels=c(1,0))
          
          confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
          accuracy <- confusion$overall['Accuracy']
          precision <- confusion$byClass['Precision']
          recall <- confusion$byClass['Recall']
          f1_score <- 2 * (precision * recall) / (precision + recall)
          
          metrics_list3[[outcome]] <- list(
            Characteristic = glue("age{ages[[age]]}"),
            AUC = auc_value,
            Accuracy = accuracy,
            Precision = precision,
            Recall = recall,
            F1_Score = f1_score,
            Test_support = nrow(urines5Test2)
          )
          
          metrics_litlist3[[age]] <- metrics_list3
          
        } else {
          
          metrics_litlist3[[age]] <- list(
            Characteristic = glue("age{ages[[age]]}"),
            AUC = NA,
            Accuracy = NA,
            Precision = NA,
            Recall = NA,
            F1_Score = NA,
            Test_support = nrow(urines5Test2)
          )
          
        }
        
      }
      
      metrics_medlist3[[seedpick]] <- metrics_litlist3
      
    }
  }
  
  metrics_biglist3[[outcome]] <- metrics_medlist3
  
}
metrics_df3 <- data.frame(matrix(nrow=0,ncol=9))
for (key in 1:length(names(metrics_biglist3))) {
  for (i in 1:length(metrics_biglist3[[key]])) {
    
    for (j in 1:length(metrics_biglist3[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist3[[key]][[i]][[j]])) {
        
        print(metrics_biglist3[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist3)[key],
          Iteration = i,
          Characteristic = metrics_biglist3[[key]][[i]][[j]][[k]]$Characteristic,
          AUC = metrics_biglist3[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist3[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist3[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist3[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist3[[key]][[i]][[j]][[k]]$F1_Score,
          Test_support = metrics_biglist3[[key]][[i]][[j]][[k]]$Test_support)
        
        colnames(metrics_df3) <- colnames(results)
        
        metrics_df3 <- data.frame(
          rbind(
            metrics_df3,results
          )
        )
        
      }
    }
  }
}
rownames(metrics_df3) <- NULL

###Write fairness metrics to CSV
metrics_protchar_df <- data.frame(rbind(metrics_df2,metrics_df3))
metrics_protchar_df <- metrics_protchar_df %>% mutate(Category = case_when(
  grepl("MALE",Characteristic)~"Gender",
  grepl("marital",Characteristic)~"Marital status",
  grepl("language",Characteristic)~"Language spoken",
  grepl("race",Characteristic)~"Race",
  grepl("^age",Characteristic)~"Age group"
),
Characteristic=case_when(
  grepl("NONMARRIED",Characteristic)~"Non-married",
  grepl("^marital_status_MARRIED$",Characteristic)~"Married",
  grepl("WHITE",Characteristic)~"White",
  grepl("BLACK",Characteristic)~"Black",
  grepl("HISPANIC",Characteristic)~"Hispanic",
  grepl("ASIAN",Characteristic)~"Asian",
  grepl("OTHER",Characteristic)~"Other",
  grepl("MALEFALSE",Characteristic)~"Female",
  grepl("MALETRUE",Characteristic)~"Male",
  grepl("ENGLISH",Characteristic)~"English",
  grepl("NONENG",Characteristic)~"Non-English",
  grepl("18",Characteristic)~"18-29",
  grepl("30",Characteristic)~"30-39",
  grepl("40",Characteristic)~"40-49",
  grepl("50",Characteristic)~"50-59",
  grepl("60",Characteristic)~"60-69",
  grepl("70",Characteristic)~"70-79",
  grepl("80",Characteristic)~"80-89",
  grepl("90",Characteristic)~"≥ 90",
)) %>% relocate(Category,.before = "Characteristic")
metrics_protchar_df <- metrics_protchar_df %>% arrange(Antimicrobial, Iteration, Category, Characteristic)
write_csv(metrics_protchar_df,"fairness_metrics.csv")

##Time sensitivity analysis
patkey <- pats %>% select(subject_id,anchor_year_group)
urines_ref <- urines_ref %>% left_join(patkey)
ind2008 <- which(urines_ref$anchor_year_group=="2008 - 2010")
ind2011 <- which(urines_ref$anchor_year_group=="2011 - 2013")
ind2014 <- which(urines_ref$anchor_year_group=="2014 - 2016")
ind2017 <- which(urines_ref$anchor_year_group=="2017 - 2019")

time_seq <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")
metrics_biglist4 <- list()
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    metrics_medlist4 <- list()
      
      for (seedpick in seq(1,6)) {
        
        set.seed(seedpick)
        
        metrics_medmedlist4 <- list()
        
        for(tim in seq_along(time_seq)) {
          
          metrics_litlist4 <- list()
        
        urines5_filtered1 <- urines5_combined[which(urines_ref$anchor_year_group==time_seq[tim]),]
        trainIndex1 <- createDataPartition(urines5_filtered1[[outcome]], p = 0.8, list = FALSE, times = 1)
        urines5Train <- urines5_filtered1[trainIndex1, ]
        
        for(tim2 in seq_along(time_seq)) {
          
          metrics_list4 <- list()
        
          iterrun <- glue("{outcome} iteration {seedpick} for {time_seq[tim]} training, {time_seq[tim2]} testing")
          print(iterrun)
          
        urines5_filtered2 <- urines5_combined[which(urines_ref$anchor_year_group==time_seq[tim2]),]
        trainIndex2 <- createDataPartition(urines5_filtered2[[outcome]], p = 0.8, list = FALSE, times = 1)
        urines5Test <- urines5_filtered2[-trainIndex2, ]
        
        predictor_columns <- colnames(urines5_predictors)
        selected_columns <- intersect(predictor_columns, colnames(urines5Train))
        train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                    label = urines5Train[[outcome]])
        test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                                   label = urines5Test[[outcome]])
        
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
        
        pred_prob_test <- predict(xgb_model, newdata = test_matrix)
        roc_result <- roc(urines5Test[[outcome]], pred_prob_test)
        auc_value <- auc(roc_result)
        print(paste("AUC-ROC:", auc_value))
        pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
        plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
        dev.off()
        
        pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
        pred_test_class <- factor(pred_test_class,levels=c(1,0))
        actual_test_class <- urines5Test[[outcome]]
        actual_test_class <- factor(actual_test_class,levels=c(1,0))
        
        confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
        accuracy <- confusion$overall['Accuracy']
        precision <- confusion$byClass['Precision']
        recall <- confusion$byClass['Recall']
        f1_score <- 2 * (precision * recall) / (precision + recall)
        
        metrics_list4[[outcome]] <- list(
          AUC = auc_value,
          Accuracy = accuracy,
          Precision = precision,
          Recall = recall,
          F1_Score = f1_score,
          Train_support = nrow(urines5Train),
          Test_support = nrow(urines5Test)
        )
        
        metrics_litlist4[[tim2]] <- metrics_list4
        names(metrics_litlist4[[tim2]]) <- time_seq[tim2]
        
      }
      
      metrics_medmedlist4[[tim]] <- metrics_litlist4
      names(metrics_medmedlist4[[tim]]) <- time_seq[tim]
      
        }
        
        metrics_medlist4[[seedpick]] <- metrics_medmedlist4
        
      }
      
  }
  
  metrics_biglist4[[outcome]] <- metrics_medlist4
  
}
metrics_df4 <- data.frame(matrix(nrow=0,ncol=11))
for (key in 1:length(names(metrics_biglist4))) {
  for (i in 1:length(metrics_biglist4[[key]])) {
    
    for (j in 1:length(metrics_biglist4[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist4[[key]][[i]][[j]])) {
        
        for (l in 1:length(metrics_biglist4[[key]][[i]][[j]][[k]])) {
          
        results <- data.frame(
          Antimicrobial = names(metrics_biglist4)[key],
          Iteration = i,
          Train_year = names(metrics_biglist4[[key]][[i]][[j]][1]),
          Test_year = names(metrics_biglist4[[key]][[i]][[j]][[k]]),
          AUC = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$AUC,
          Accuracy = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Accuracy,
          Precision = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Precision,
          Recall = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Recall,
          F1_score = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$F1_Score,
          Train_support = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Train_support,
          Test_support = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Test_support)
        
        colnames(metrics_df4) <- colnames(results)
        
        metrics_df4 <- data.frame(
          rbind(
            metrics_df4,results
          )
        )
        }
      }
    }
  }
}
rownames(metrics_df4) <- NULL
write_csv(metrics_df4,"time_sens_metrics.csv")

##CDI and toxicity dataset preprocessing
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()
abx <- data.frame(rbind(train_abx,test_abx))
hadm_key <- hadm %>% select(subject_id,race,language,marital_status) %>% 
  distinct(subject_id,.keep_all=T)
age_key <- pats %>% select(subject_id,anchor_age) %>% mutate(age_grp=case_when(
  anchor_age<30~18,anchor_age>=30&anchor_age<40~30,
  anchor_age>=40&anchor_age<50~40,
  anchor_age>=50&anchor_age<60~50,
  anchor_age>=60&anchor_age<70~60,
  anchor_age>=70&anchor_age<80~70,
  anchor_age>=80&anchor_age<90~80,
  anchor_age>=90~90
)) %>% select(-anchor_age) %>% distinct(subject_id,.keep_all = T)
dems_key <- left_join(hadm_key,age_key,by="subject_id")
abx <- abx %>% left_join(dems_key,by="subject_id")
abx <- abx %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                       TRUE~marital_status))
abx_outcomes <- abx %>%
  select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                     overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
abx_predictors <- abx %>% select(pHADM:age65,prAKI:pDIAB,pCARD:curr_service,pICU:pSEPSIS,ob_freq,highCRP,
                                 temperature:dbp,inpatient_7d,
                                 ab_name_Ampicillin_Ceftriaxone:ab_name_Ampicillin,race:age_grp)
set.seed(123)
dummies <- dummyVars(" ~ .", data = abx_predictors)
abx_predictors <- predict(dummies, newdata = abx_predictors)
abx_combined <- as.data.frame(cbind(abx_outcomes, abx_predictors))

###Read-in chosen model parameters (CDI/toxicity)
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

###Model stability analysis (CDI and toxicity prediction)

###Iterate training and testing over training dataset sizes (CDI tox)
tr_size_seq <- c(0.02,0.06,0.10,0.14)
metrics_biglist <- list()
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    metrics_medlist <- list()
    
    for(siz in seq_along(tr_size_seq)) {
      
      metrics_litlist <- list()
      
      for (seedpick in seq(1,6)) {
        
        metrics_list <- list()
        
        iterrun <- glue("{outcome} iteration {seedpick} for {tr_size_seq[siz]}")
        print(iterrun)
        
        set.seed(seedpick)
        trainIndex <- createDataPartition(abx_combined[[outcome]], p = tr_size_seq[siz], list = FALSE, times = 1)
        abxTrain <- abx_combined[trainIndex, ]
        abxTest <- abx_combined[-trainIndex, ]
        
        predictor_columns <- colnames(abx_predictors)
        selected_columns <- intersect(predictor_columns, colnames(abxTrain))
        train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                                    label = abxTrain[[outcome]])
        test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                                   label = abxTest[[outcome]])
        
        params <- list(
          objective = "binary:logistic",
          eval_metric = "auc",
          eta = cdi_tox_final_bestparams[[outcome]]$eta,
          max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
          min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
          subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
          colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
        )
        
        print("Training...")
        
        xgb_model <- xgb.train(
          params = params,
          data = train_matrix,
          nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
        )
        
        pred_prob_test <- predict(xgb_model, newdata = test_matrix)
        roc_result <- roc(abxTest[[outcome]], pred_prob_test)
        auc_value <- auc(roc_result)
        print(paste("AUC-ROC:", auc_value))
        pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
        plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
        dev.off()
        
        pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
        pred_test_class <- factor(pred_test_class,levels=c(1,0))
        actual_test_class <- abxTest[[outcome]]
        actual_test_class <- factor(actual_test_class,levels=c(1,0))
        
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
        
        metrics_litlist[[seedpick]] <- metrics_list
        
      }
      
      metrics_medlist[[siz]] <- metrics_litlist
      names(metrics_medlist[[siz]]) <- tr_size_seq[siz]
      
    }
  }
  
  metrics_biglist[[outcome]] <- metrics_medlist
  
}

###Write stability metrics to CSV (CDI tox)
metrics_df <- data.frame(matrix(nrow=0,ncol=7))
for (key in 1:length(names(metrics_biglist))) {
  for (i in 1:length(metrics_biglist[[key]])) {
    
    for (j in 1:length(metrics_biglist[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist[[key]][[i]][[j]])) {
        
        print(metrics_biglist[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist)[key],
          Training_size = names(metrics_biglist[[key]][[i]][1]),
          AUC = metrics_biglist[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist[[key]][[i]][[j]][[k]]$F1_Score)
        
        colnames(metrics_df) <- colnames(results)
        
        metrics_df <- data.frame(
          rbind(
            metrics_df,results
          )
        )
        
      }
    }
  }
}
rownames(metrics_df) <- NULL
write_csv(metrics_df,"cdi_tox_stability_metrics.csv")

##Model fairness analysis (CDI tox)

###Protected characteristics excluding age (CDI tox)
ur_prot <- abx_combined %>% mutate(
  language_NONENG = case_when(languageENGLISH!=1~1,TRUE~0),
  race_ASIAN = rowSums(select(., contains("raceASIAN"))) > 0,
  race_BLACK = rowSums(select(., contains("raceBLACK"))) > 0,
  race_HISPANIC = rowSums(select(., contains("raceHISPANIC"))) > 0,
  race_OTHER = rowSums(select(., matches("(raceOTHER|racePORTUGUESE|raceSOUTH|raceNATIVE|raceAMERICAN|raceMULTIPLE)"))) > 0,
  race_WHITE = rowSums(select(., contains("raceWHITE"))) > 0,
  marital_status_MARRIED = rowSums(select(., matches("(marital_statusMARRIED)"))) > 0,
  marital_status_NONMARRIED = case_when(marital_statusMARRIED!=1~1,TRUE~0)
)
protchar_index <- which(grepl("(^marital_status_|^language|^MALE|^race_)",colnames(ur_prot))&
                          !grepl("(UNKNOWN|Other|\\?)",colnames(ur_prot)))
protchars <- colnames(ur_prot)[protchar_index]
metrics_biglist2 <- list()
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    metrics_medlist2 <- list()
    
    for (seedpick in seq(1,6)) {
      
      set.seed(seedpick)
      
      trainIndex <- createDataPartition(abx_combined[[outcome]], p = 0.8, list = FALSE, times = 1)
      abxTrain <- abx_combined[trainIndex, ]
      abxTest <- abx_combined[-trainIndex, ]
      abxTest <- abxTest %>% mutate(
        language_NONENG = case_when(languageENGLISH!=1~1,TRUE~0),
        race_ASIAN = rowSums(select(., contains("raceASIAN"))) > 0,
        race_BLACK = rowSums(select(., contains("raceBLACK"))) > 0,
        race_HISPANIC = rowSums(select(., contains("raceHISPANIC"))) > 0,
        race_OTHER = rowSums(select(., matches("(raceOTHER|racePORTUGUESE|raceSOUTH|raceNATIVE|raceAMERICAN|raceMULTIPLE)"))) > 0,
        race_WHITE = rowSums(select(., contains("raceWHITE"))) > 0,
        marital_status_MARRIED = rowSums(select(., matches("(marital_statusMARRIED)"))) > 0,
        marital_status_NONMARRIED = case_when(marital_statusMARRIED!=1~1,TRUE~0)
      )
      
      predictor_columns <- colnames(abx_predictors)
      selected_columns <- intersect(predictor_columns, colnames(abxTrain))
      train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                                  label = abxTrain[[outcome]])
      test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                                 label = abxTest[[outcome]])
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = cdi_tox_final_bestparams[[outcome]]$eta,
        max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
        min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
        subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
        colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
      )
      
      print("Training...")
      
      xgb_model <- xgb.train(
        params = params,
        data = train_matrix,
        nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
      )
      
      metrics_litlist2 <- list()
      
      for(protchar in seq_along(protchar_index)) {
        
        iterrun <- glue("{outcome} iteration {seedpick} for {protchars[protchar]}")
        print(iterrun)
        
        metrics_list2 <- list()
        
        abxTest2 <- abxTest[abxTest[protchar_index[protchar]]==1,]
        
        if (length(abxTest2[[outcome]] %>% unique()) > 1) {
          
          
          test_matrix2 <- xgb.DMatrix(data = as.matrix(abxTest2 %>% select(all_of(selected_columns))), 
                                      label = abxTest2[[outcome]])
          
          pred_prob_test <- predict(xgb_model, newdata = test_matrix2)
          roc_result <- roc(abxTest2[[outcome]], pred_prob_test)
          auc_value <- auc(roc_result)
          print(paste("AUC-ROC:", auc_value))
          pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
          plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
          dev.off()
          
          pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
          pred_test_class <- factor(pred_test_class,levels=c(1,0))
          actual_test_class <- abxTest2[[outcome]]
          actual_test_class <- factor(actual_test_class,levels=c(1,0))
          
          confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
          accuracy <- confusion$overall['Accuracy']
          precision <- confusion$byClass['Precision']
          recall <- confusion$byClass['Recall']
          f1_score <- 2 * (precision * recall) / (precision + recall)
          
          metrics_list2[[outcome]] <- list(
            Characteristic = protchars[[protchar]],
            AUC = auc_value,
            Accuracy = accuracy,
            Precision = precision,
            Recall = recall,
            F1_Score = f1_score,
            Test_support = nrow(abxTest2)
          )
          
          metrics_litlist2[[protchar]] <- metrics_list2
          
        } else {
          
          metrics_litlist2[[protchar]] <- list(
            Characteristic = protchars[[protchar]],
            AUC = NA,
            Accuracy = NA,
            Precision = NA,
            Recall = NA,
            F1_Score = NA,
            Test_support = nrow(abxTest2)
          )
          
        }
        
      }
      
      metrics_medlist2[[seedpick]] <- metrics_litlist2
    }
  }
  
  metrics_biglist2[[outcome]] <- metrics_medlist2
  
}
metrics_df2 <- data.frame(matrix(nrow=0,ncol=9))
for (key in 1:length(names(metrics_biglist2))) {
  for (i in 1:length(metrics_biglist2[[key]])) {
    
    for (j in 1:length(metrics_biglist2[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist2[[key]][[i]][[j]])) {
        
        print(metrics_biglist2[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist2)[[key]],
          Iteration = i,
          Characteristic = metrics_biglist2[[key]][[i]][[j]][[k]]$Characteristic,
          AUC = metrics_biglist2[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist2[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist2[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist2[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist2[[key]][[i]][[j]][[k]]$F1_Score,
          Test_support = metrics_biglist2[[key]][[i]][[j]][[k]]$Test_support)
        
        colnames(metrics_df2) <- colnames(results)
        
        metrics_df2 <- data.frame(
          rbind(
            metrics_df2,results
          )
        )
        
      }
    }
  }
}
rownames(metrics_df2) <- NULL
write_csv(metrics_df2,"cdi_tox_fairness_metrics1.csv")

###Age (CDI tox)
ages <- abx_combined %>% distinct(age_grp) %>% unlist() %>% sort()
metrics_biglist3 <- list()
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    metrics_medlist3 <- list()
    
    for (seedpick in seq(1,6)) {
      
      set.seed(seedpick)
      
      trainIndex <- createDataPartition(abx_combined[[outcome]], p = 0.8, list = FALSE, times = 1)
      abxTrain <- abx_combined[trainIndex, ]
      abxTest <- abx_combined[-trainIndex, ]
      
      predictor_columns <- colnames(abx_predictors)
      selected_columns <- intersect(predictor_columns, colnames(abxTrain))
      train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                                  label = abxTrain[[outcome]])
      test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                                 label = abxTest[[outcome]])
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = cdi_tox_final_bestparams[[outcome]]$eta,
        max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
        min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
        subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
        colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
      )
      
      print("Training...")
      
      xgb_model <- xgb.train(
        params = params,
        data = train_matrix,
        nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
      )
      
      metrics_litlist3 <- list()
      
      for(age in seq_along(ages)) {
        
        iterrun <- glue("{outcome} iteration {seedpick} for age {ages[age]}")
        print(iterrun)
        
        metrics_list3 <- list()
        
        if (sum(!is.na(abx_combined[[outcome]])) > 0) {
          
          abxTest2 <- abxTest %>% filter(age_grp==ages[age])
          test_matrix2 <- xgb.DMatrix(data = as.matrix(abxTest2 %>% select(all_of(selected_columns))), 
                                      label = abxTest2[[outcome]])
          
          pred_prob_test <- predict(xgb_model, newdata = test_matrix2)
          roc_result <- roc(abxTest2[[outcome]], pred_prob_test)
          auc_value <- auc(roc_result)
          print(paste("AUC-ROC:", auc_value))
          pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
          plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
          dev.off()
          
          pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
          pred_test_class <- factor(pred_test_class,levels=c(1,0))
          actual_test_class <- abxTest2[[outcome]]
          actual_test_class <- factor(actual_test_class,levels=c(1,0))
          
          confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
          accuracy <- confusion$overall['Accuracy']
          precision <- confusion$byClass['Precision']
          recall <- confusion$byClass['Recall']
          f1_score <- 2 * (precision * recall) / (precision + recall)
          
          metrics_list3[[outcome]] <- list(
            Characteristic = glue("age{ages[[age]]}"),
            AUC = auc_value,
            Accuracy = accuracy,
            Precision = precision,
            Recall = recall,
            F1_Score = f1_score,
            Test_support = nrow(abxTest2)
          )
          
          metrics_litlist3[[age]] <- metrics_list3
          
        } else {
          
          metrics_litlist3[[age]] <- list(
            Characteristic = glue("age{ages[[age]]}"),
            AUC = NA,
            Accuracy = NA,
            Precision = NA,
            Recall = NA,
            F1_Score = NA,
            Test_support = nrow(abxTest2)
          )
          
        }
        
      }
      
      metrics_medlist3[[seedpick]] <- metrics_litlist3
      
    }
  }
  
  metrics_biglist3[[outcome]] <- metrics_medlist3
  
}
metrics_df3 <- data.frame(matrix(nrow=0,ncol=9))
for (key in 1:length(names(metrics_biglist3))) {
  for (i in 1:length(metrics_biglist3[[key]])) {
    
    for (j in 1:length(metrics_biglist3[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist3[[key]][[i]][[j]])) {
        
        print(metrics_biglist3[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist3)[key],
          Iteration = i,
          Characteristic = metrics_biglist3[[key]][[i]][[j]][[k]]$Characteristic,
          AUC = metrics_biglist3[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist3[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist3[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist3[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist3[[key]][[i]][[j]][[k]]$F1_Score,
          Test_support = metrics_biglist3[[key]][[i]][[j]][[k]]$Test_support)
        
        colnames(metrics_df3) <- colnames(results)
        
        metrics_df3 <- data.frame(
          rbind(
            metrics_df3,results
          )
        )
        
      }
    }
  }
}
rownames(metrics_df3) <- NULL

###Write fairness metrics to CSV (CDI tox)
metrics_protchar_df <- data.frame(rbind(metrics_df2,metrics_df3))
metrics_protchar_df <- metrics_protchar_df %>% mutate(Category = case_when(
  grepl("MALE",Characteristic)~"Gender",
  grepl("marital",Characteristic)~"Marital status",
  grepl("language",Characteristic)~"Language spoken",
  grepl("race",Characteristic)~"Race",
  grepl("^age",Characteristic)~"Age group"
),
Characteristic=case_when(
  grepl("NONMARRIED",Characteristic)~"Non-married",
  grepl("^marital_status_MARRIED$",Characteristic)~"Married",
  grepl("WHITE",Characteristic)~"White",
  grepl("BLACK",Characteristic)~"Black",
  grepl("HISPANIC",Characteristic)~"Hispanic",
  grepl("ASIAN",Characteristic)~"Asian",
  grepl("OTHER",Characteristic)~"Other",
  grepl("MALEFALSE",Characteristic)~"Female",
  grepl("MALETRUE",Characteristic)~"Male",
  grepl("ENGLISH",Characteristic)~"English",
  grepl("NONENG",Characteristic)~"Non-English",
  grepl("18",Characteristic)~"18-29",
  grepl("30",Characteristic)~"30-39",
  grepl("40",Characteristic)~"40-49",
  grepl("50",Characteristic)~"50-59",
  grepl("60",Characteristic)~"60-69",
  grepl("70",Characteristic)~"70-79",
  grepl("80",Characteristic)~"80-89",
  grepl("90",Characteristic)~"≥ 90",
)) %>% relocate(Category,.before = "Characteristic")
metrics_protchar_df <- metrics_protchar_df %>% arrange(Antimicrobial, Iteration, Category, Characteristic)
write_csv(metrics_protchar_df,"cdi_tox_fairness_metrics.csv")

##Time sensitivity analysis (CDI tox)
patkey <- pats %>% select(subject_id,anchor_year_group)
urines_ref <- urines_ref %>% left_join(patkey)
ind2008 <- which(urines_ref$anchor_year_group=="2008 - 2010")
ind2011 <- which(urines_ref$anchor_year_group=="2011 - 2013")
ind2014 <- which(urines_ref$anchor_year_group=="2014 - 2016")
ind2017 <- which(urines_ref$anchor_year_group=="2017 - 2019")

time_seq <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")
metrics_biglist4 <- list()
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    metrics_medlist4 <- list()
    
    for (seedpick in seq(1,6)) {
      
      set.seed(seedpick)
      
      metrics_medmedlist4 <- list()
      
      for(tim in seq_along(time_seq)) {
        
        metrics_litlist4 <- list()
        
        abx_filtered1 <- abx_combined[which(urines_ref$anchor_year_group==time_seq[tim]),]
        trainIndex1 <- createDataPartition(abx_filtered1[[outcome]], p = 0.8, list = FALSE, times = 1)
        abxTrain <- abx_filtered1[trainIndex1, ]
        
        for(tim2 in seq_along(time_seq)) {
          
          metrics_list4 <- list()
          
          iterrun <- glue("{outcome} iteration {seedpick} for {time_seq[tim]} training, {time_seq[tim2]} testing")
          print(iterrun)
          
          abx_filtered2 <- abx_combined[which(urines_ref$anchor_year_group==time_seq[tim2]),]
          trainIndex2 <- createDataPartition(abx_filtered2[[outcome]], p = 0.8, list = FALSE, times = 1)
          abxTest <- abx_filtered2[-trainIndex2, ]
          
          predictor_columns <- colnames(abx_predictors)
          selected_columns <- intersect(predictor_columns, colnames(abxTrain))
          train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(all_of(selected_columns))), 
                                      label = abxTrain[[outcome]])
          test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(all_of(selected_columns))), 
                                     label = abxTest[[outcome]])
          
          params <- list(
            objective = "binary:logistic",
            eval_metric = "auc",
            eta = cdi_tox_final_bestparams[[outcome]]$eta,
            max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
            min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
            subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
            colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
          )
          
          print("Training...")
          
          xgb_model <- xgb.train(
            params = params,
            data = train_matrix,
            nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
          )
          
          pred_prob_test <- predict(xgb_model, newdata = test_matrix)
          roc_result <- roc(abxTest[[outcome]], pred_prob_test)
          auc_value <- auc(roc_result)
          print(paste("AUC-ROC:", auc_value))
          pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
          plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
          dev.off()
          
          pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
          pred_test_class <- factor(pred_test_class,levels=c(1,0))
          actual_test_class <- abxTest[[outcome]]
          actual_test_class <- factor(actual_test_class,levels=c(1,0))
          
          confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
          accuracy <- confusion$overall['Accuracy']
          precision <- confusion$byClass['Precision']
          recall <- confusion$byClass['Recall']
          f1_score <- 2 * (precision * recall) / (precision + recall)
          
          metrics_list4[[outcome]] <- list(
            AUC = auc_value,
            Accuracy = accuracy,
            Precision = precision,
            Recall = recall,
            F1_Score = f1_score,
            Train_support = nrow(abxTrain),
            Test_support = nrow(abxTest)
          )
          
          metrics_litlist4[[tim2]] <- metrics_list4
          names(metrics_litlist4[[tim2]]) <- time_seq[tim2]
          
        }
        
        metrics_medmedlist4[[tim]] <- metrics_litlist4
        names(metrics_medmedlist4[[tim]]) <- time_seq[tim]
        
      }
      
      metrics_medlist4[[seedpick]] <- metrics_medmedlist4
      
    }
    
  }
  
  metrics_biglist4[[outcome]] <- metrics_medlist4
  
}
metrics_df4 <- data.frame(matrix(nrow=0,ncol=11))
for (key in 1:length(names(metrics_biglist4))) {
  for (i in 1:length(metrics_biglist4[[key]])) {
    
    for (j in 1:length(metrics_biglist4[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist4[[key]][[i]][[j]])) {
        
        for (l in 1:length(metrics_biglist4[[key]][[i]][[j]][[k]])) {
          
          results <- data.frame(
            Antimicrobial = names(metrics_biglist4)[key],
            Iteration = i,
            Train_year = names(metrics_biglist4[[key]][[i]][[j]][1]),
            Test_year = names(metrics_biglist4[[key]][[i]][[j]][[k]]),
            AUC = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$AUC,
            Accuracy = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Accuracy,
            Precision = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Precision,
            Recall = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Recall,
            F1_score = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$F1_Score,
            Train_support = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Train_support,
            Test_support = metrics_biglist4[[key]][[i]][[j]][[k]][[l]]$Test_support)
          
          colnames(metrics_df4) <- colnames(results)
          
          metrics_df4 <- data.frame(
            rbind(
              metrics_df4,results
            )
          )
        }
      }
    }
  }
}
rownames(metrics_df4) <- NULL
write_csv(metrics_df4,"cdi_tox_time_sens_metrics.csv")

metrics_dfA <- read_csv("stability_metrics.csv")
metrics_df<- data.frame(rbind(metrics_dfA,metrics_df))
metrics_df <- metrics_df %>% rename(Model="Antimicrobial")
write_csv(metrics_df,"overall_stability_metrics.csv")

metrics_protchar_dfA <- read_csv("fairness_metrics.csv")
metrics_protchar_df<- data.frame(rbind(metrics_protchar_dfA,metrics_protchar_df))
metrics_protchar_df <- metrics_protchar_df %>% rename(Model="Antimicrobial")
write_csv(metrics_protchar_df,"overall_fairness_metrics.csv")

metrics_df4a <- read_csv("time_sens_metrics.csv")
metrics_df4 <- read_csv("cdi_tox_time_sens_metrics.csv")
metrics_df4 <- data.frame(rbind(metrics_df4a,metrics_df4))
metrics_df4 <- metrics_df4 %>% rename(Model="Antimicrobial")
write_csv(metrics_df4,"overall_time_sens_metrics.csv")

##Model testing plots

###Stability analysis
stability_plot(metrics_df,AUC,"AUC-ROC")
stability_plot(metrics_df,Accuracy,"Accuracy")
stability_plot(metrics_df,Precision,"Precision")
stability_plot(metrics_df,Recall,"Recall")
stability_plot(metrics_df,F1_score,"F1 score")

###Fairness analysis
protchar_plot(metrics_protchar_df,"Age group",AUC,"AUC-ROC","age groups")
protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",AUC,"AUC-ROC",
              "genders, languages, and marital statuses")
protchar_plot(metrics_protchar_df,"Race",AUC,"AUC-ROC","racial groups")
protchar_plot(metrics_protchar_df,"Age group",Accuracy,"Accuracy","age groups")
protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",Accuracy,"Accuracy",
              "genders, languages, and marital statuses")
protchar_plot(metrics_protchar_df,"Race",Accuracy,"Accuracy","racial groups")
protchar_plot(metrics_protchar_df,"Age group",Precision,"Precision","age groups")
protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",Precision,"Precision",
              "genders, languages, and marital statuses")
protchar_plot(metrics_protchar_df,"Race",Precision,"Precision","racial groups")
protchar_plot(metrics_protchar_df,"Age group",Recall,"Recall","age groups")
protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",Recall,"Recall",
              "genders, languages, and marital statuses")
protchar_plot(metrics_protchar_df,"Race",Recall,"Recall","racial groups")
protchar_plot(metrics_protchar_df,"Age group",F1_score,"F1 score","age groups")
protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",F1_score,"F1 score",
              "genders, languages, and marital statuses")
protchar_plot(metrics_protchar_df,"Race",F1_score,"F1 score","racial groups")

###Time sensitivity analysis
yearlist <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")
for (i in seq_along(yearlist)) {

  timesens_plot(metrics_df4,yearlist[i],AUC,"AUC-ROC")
  timesens_plot(metrics_df4,yearlist[i],Accuracy,"Accuracy")
  timesens_plot(metrics_df4,yearlist[i],Precision,"Precision")
  timesens_plot(metrics_df4,yearlist[i],Recall,"Recall")
  timesens_plot(metrics_df4,yearlist[i],F1_score,"F1 score")
  
}

##Hyperparameter dataframe
hyp_names = c("Learning rate","Maximum tree depth",
              "Minimum child weight","Row subsample",
              "Column subsample","N training rounds")

hyp_tab <- tibble(AMP=final_bestparams$AMP %>% unlist(),
       SAM=final_bestparams$SAM %>% unlist(),
       TZP=final_bestparams$TZP %>% unlist(),
       CZO=final_bestparams$CZO %>% unlist(),
       CRO=final_bestparams$CRO %>% unlist(),
       CAZ=final_bestparams$CAZ %>% unlist(),
       FEP=final_bestparams$FEP %>% unlist(),
       MEM=final_bestparams$MEM %>% unlist(),
       CIP=final_bestparams$CIP %>% unlist(),
       GEN=final_bestparams$GEN %>% unlist(),
       SXT=final_bestparams$SXT %>% unlist(),
       NIT=final_bestparams$NIT %>% unlist(),
       VAN=final_bestparams$VAN %>% unlist(),
       AMP_CRO=final_bestparams$AMP_CRO %>% unlist(),
       AMP_GEN=final_bestparams$AMP_GEN %>% unlist(),
       AMP_VAN=final_bestparams$AMP_VAN %>% unlist(),
       SAM_VAN=final_bestparams$SAM_VAN %>% unlist(),
       TZP_CIP=final_bestparams$TZP_CIP %>% unlist(),
       TZP_SXT=final_bestparams$TZP_SXT %>% unlist(),
       TZP_VAN=final_bestparams$TZP_VAN %>% unlist(),
       CZO_CIP=final_bestparams$CZO_CIP %>% unlist(),
       CZO_SXT=final_bestparams$CZO_SXT %>% unlist(),
       CZO_VAN=final_bestparams$CZO_VAN %>% unlist(),
       CRO_FEP=final_bestparams$CRO_FEP %>% unlist(),
       CRO_CIP=final_bestparams$CRO_CIP %>% unlist(),
       CRO_SXT=final_bestparams$CRO_SXT %>% unlist(),
       CRO_VAN=final_bestparams$CRO_VAN %>% unlist(),
       CAZ_VAN=final_bestparams$CAZ_VAN %>% unlist(),
       FEP_CIP=final_bestparams$FEP_CIP %>% unlist(),
       FEP_SXT=final_bestparams$FEP_SXT %>% unlist(),
       FEP_VAN=final_bestparams$FEP_VAN %>% unlist(),
       MEM_SXT=final_bestparams$MEM_SXT %>% unlist(),
       MEM_VAN=final_bestparams$MEM_VAN %>% unlist(),
       CIP_SXT=final_bestparams$CIP_SXT %>% unlist(),
       CIP_VAN=final_bestparams$CIP_VAN %>% unlist(),
       GEN_VAN=final_bestparams$GEN_VAN %>% unlist(),
       SXT_VAN=final_bestparams$SXT_VAN %>% unlist(),
       CDI=cdi_tox_final_bestparams$CDI %>% unlist(),
       Toxicity=cdi_tox_final_bestparams$overall_tox %>% unlist()) %>% 
  dplyr::slice(-c(1:2))

replace_values <- function(column, map) {
  flipped_map <- setNames(names(map), map)
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(flipped_map)) flipped_map[[x]] else x)
}
colnames(hyp_tab) <- replace_values(colnames(hyp_tab),combined_antimicrobial_map) %>% 
  str_replace_all("_"," & ") %>% str_replace_all("-","/")
  
hyp_tab <- hyp_tab %>% mutate(Hyperparameter=hyp_names) %>% 
  relocate(Hyperparameter,.before=1) %>% t() %>% data.frame() 
hyp_tab$Model <- rownames(hyp_tab)
hyp_tab <- hyp_tab %>% relocate(Model,.before = 1)
colnames(hyp_tab) <- hyp_tab[1,]
hyp_tab <- hyp_tab %>% rename(Model="Hyperparameter") %>% 
  dplyr::slice(-1)
hyp_tab <- hyp_tab %>% tibble()
write_csv(hyp_tab,"hyp_tab.csv")
hyp_abs <- c(c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN") %>% ab_name(),"CDI","Toxicity")
hyp_tab_singles <- hyp_tab %>% filter(Model %in% hyp_abs)
write_csv(hyp_tab_singles,"hyp_tab_singles.csv")

##Values for results

###Stability analysis
stabmets <- read_csv("overall_stability_metrics.csv")
max_stabdifs <- stabmets %>% group_by(Model,Training_size) %>%
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% 
  summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
            max_sd=max(sd_AUC)) %>% ungroup()
max_stabdifs %>% arrange(desc(maxAUC_meandif))
max_stabdifs %>% arrange(desc(max_sd))
stabmets %>% group_by(Model,Training_size) %>% 
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model=="CDI")
stabmets %>% group_by(Model,Training_size) %>% 
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model=="VAN")

###Year group cluster analysis
timemets <- read_csv("overall_time_sens_metrics.csv")
max_timemets <- timemets %>% group_by(Model,Train_year,Test_year) %>%
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% group_by(Model) %>%  
  summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
            max_sd=max(sd_AUC)) %>% ungroup()
max_timemets %>% arrange(desc(maxAUC_meandif))
max_timemets %>% arrange(desc(max_sd))
timemets %>% group_by(Model,Train_year,Test_year) %>% 
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model=="overall_tox") %>% 
  arrange(desc(mean_AUC))
timemets %>% group_by(Model,Train_year,Test_year) %>% 
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model=="CDI") %>% 
  arrange(desc(sd_AUC))

###Fairness analysis
fairmets <- read_csv("overall_fairness_metrics.csv")
fairnessprinter1 <- function(cat) {
max_fairmets <- fairmets %>% filter(Category==cat) %>% group_by(Model,Characteristic) %>%
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% 
  summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
            max_sd=max(sd_AUC)) %>% ungroup()
max_fairmets %>% arrange(desc(maxAUC_meandif)) %>% print()
max_fairmets %>% arrange(desc(max_sd)) %>% print()
}
fairnessprinter2 <- function(cat,cat1,cat2) {
fairmets %>% filter(Category==cat) %>% group_by(Model,Characteristic) %>% 
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model==cat1) %>% 
  arrange(desc(mean_AUC)) %>% print()
fairmets %>% filter(Category==cat) %>% group_by(Model,Characteristic) %>%
  summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% filter(Model==cat2) %>% 
  arrange(desc(sd_AUC)) %>% print()
}
fairnessprinter1("Age group")
fairnessprinter2("Age group","TZP","TZP")
fairnessprinter1("Race")
fairnessprinter2("Race","TZP","TZP")
fairnessprinter1("Gender")
fairnessprinter2("Gender","TZP","TZP")
