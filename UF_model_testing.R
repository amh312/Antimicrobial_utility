#MODEL ADDITIONAL TESTING

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
  
  df <- df %>% mutate(Model=case_when(Model=="overall_tox"~"Toxicity",
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

##Read in

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
tr_size_seq <- c(0.02,0.06,0.08,0.12)
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
          actual_test_class <- urines5Test2[[outcome]]
          
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
          actual_test_class <- urines5Test2[[outcome]]
          
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
        actual_test_class <- urines5Test[[outcome]]
        
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
age_key <- pats %>% select(subject_id,anchor_age) %>% mutate(anchor_age=case_when(
  anchor_age<30~18,anchor_age>=30&anchor_age<40~30,
  anchor_age>=40&anchor_age<50~40,
  anchor_age>=50&anchor_age<60~50,
  anchor_age>=60&anchor_age<70~60,
  anchor_age>=70&anchor_age<80~70,
  anchor_age>=80&anchor_age<90~80,
  anchor_age>=90~90
)) %>% distinct(subject_id,.keep_all = T)
dems_key <- left_join(hadm_key,age_key,by="subject_id")
abx <- abx %>% left_join(dems_key,by="subject_id")
abx <- abx %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                       TRUE~marital_status))
abx_outcomes <- abx %>%
  select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                     overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
abx_predictors <- abx %>% select(pHADM:age65,prAKI:pDIAB,pCARD:pSEPSIS,temperature:dbp,pc_dyspnea:pc_fever,
                                 abx_name_Ampicillin_Ceftriaxone:ob_freq,highCRP,
                                 race:anchor_age)
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
tr_size_seq <- c(0.02,0.06,0.08,0.12)
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
        actual_test_class <- abxTest[[outcome]]
        
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
          actual_test_class <- abxTest2[[outcome]]
          
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
ages <- abx_combined %>% distinct(anchor_age) %>% unlist() %>% sort()
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
          
          abxTest2 <- abxTest %>% filter(anchor_age==ages[age])
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
          actual_test_class <- abxTest2[[outcome]]
          
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
          actual_test_class <- abxTest[[outcome]]
          
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

##Fairness compensation hyperparameter tuning

###MALE - cephalosporins, meropenem, gentamicin
chosen_outcomes <- c("TZP","CZO","CRO","CAZ","FEP","MEM","GEN")

num_samples <- 10
max_depth_range <- c(2,9)
min_child_weight_range <- c(1, 10)
lhs_sample <- randomLHS(num_samples, 2)
max_depth <- round(lhs_sample[, 1] * (max_depth_range[2] - max_depth_range[1]) + max_depth_range[1])
min_child_weight <- round(lhs_sample[, 2] * (min_child_weight_range[2] - min_child_weight_range[1]) + min_child_weight_range[1])
parameter_grid <- data.frame(max_depth = max_depth, min_child_weight = min_child_weight)
max_child_bestparams_male <- c()
best_auc <- 0
for (outcome in chosen_outcomes) {
  
  if (sum(!is.na(urines5_outcomes[[outcome]])) > 0) {
    
    set.seed(123)
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    urines5Test <- urines5Test %>% filter(MALETRUE==1)
    
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
    
    max_child_bestparams_male[[outcome]] <- best_params
    
  }
}
for (i in 1:length(max_child_bestparams_male)) {
  
  maxy_male <- data.frame(max_child_bestparams_male[i])
  
  write_csv(maxy_male,glue("max_child_{chosen_outcomes[i]}_male.csv"))
  
}

num_samples <- 10
subsample_range <- c(0.5,1)
colsample_bytree_range <- c(0.5, 1)
lhs_sample <- randomLHS(num_samples, 2)
subsample <- round(lhs_sample[, 1] * (subsample_range[2] - subsample_range[1]) + subsample_range[1],2)
colsample_bytree <- round(lhs_sample[, 2] * (colsample_bytree_range[2] - colsample_bytree_range[1]) + colsample_bytree_range[1],2)
parameter_grid <- data.frame(subsample = subsample, colsample_bytree = colsample_bytree)
print(parameter_grid)
col_sub_bestparams_male <- c()
best_auc <- 0
for (outcome in chosen_outcomes) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    set.seed(123)
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    urines5Test <- urines5Test %>% filter(MALETRUE==1)
    
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
        max_depth = max_child_bestparams_male[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams_male[[outcome]]$min_child_weight,
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
    
    col_sub_bestparams_male[[outcome]] <- best_params
    
  }
}
for (i in 1:length(col_sub_bestparams_male)) {
  
  coly <- data.frame(col_sub_bestparams_male[i])
  
  write_csv(coly,glue("col_sub_male_{chosen_outcomes[i]}.csv"))
  
}

parameter_list <- c(0.1,0.05,0.01,0.001)
best_auc <- 0
final_bestparams_male <- c()
for (outcome in chosen_outcomes) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    set.seed(123)
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    urines5Test <- urines5Test %>% filter(MALETRUE==1)
    
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
        max_depth = max_child_bestparams_male[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams_male[[outcome]]$min_child_weight,
        subsample = col_sub_bestparams_male[[outcome]]$subsample,
        colsample_bytree = col_sub_bestparams_male[[outcome]]$colsample_bytree
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
    
    final_bestparams_male[[outcome]] <- best_params
    final_bestparams_male[[outcome]]$best_nrounds <- best_nrounds
    
  }
}
for (i in 1:length(final_bestparams_male)) {
  
  param <- data.frame(final_bestparams_male[i])
  
  write_csv(param,glue("final_params_{chosen_outcomes[i]}.csv"))
  
}

test_probs_df <- data.frame(matrix(nrow=floor(nrow(urines5_combined)*0.2),ncol=0))
micro_probs_df <- data.frame(matrix(nrow=nrow(ur_xg_combined),ncol=0))
aucs <- data.frame(matrix(nrow=1,ncol=0))
shap_summary_tables <- list()
metrics_list <- list() 
for (outcome in chosen_outcomes) {
  
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
    pdf(glue("{outcome}_xg_roc_male.pdf"), width = 10, height = 10)
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
  
  shappy_male <- data.frame(shap_summary_tables[i]) %>% mutate(across(2, ~ round(., 3))) %>% filter(if_any(2, ~ . != 0))
  
  colnames(shappy_male) <- c("Feature","Variable")
  
  write_csv(shappy_male,glue("SHAP_male_{chosen_outcomes[i]}.csv"))
  
}
for (i in 1:length(metrics_list)) {
  
  metricky_male <- data.frame(metrics_list[i])
  
  write_csv(metricky_male,glue("metrics_male_{chosen_outcomes[i]}.csv"))
  
}

replace_values <- function(column, map) {
  flipped_map <- setNames(names(map), map)
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(flipped_map)) flipped_map[[x]] else x)
}
micro_probs_df$micro_specimen_id <- ur_xg$micro_specimen_id
micro_probs_df$subject_id <- ur_xg$subject_id
write_csv(micro_probs_df,"male_micro_probs_df.csv")
probs_df_male <- micro_probs_df %>% melt(id.vars=c("micro_specimen_id","subject_id")) %>% 
  mutate(variable=replace_values(variable, combined_antimicrobial_map)) %>% 
  rename(S="value",Antimicrobial="variable") %>% mutate(I=NA,R=1-S,NT=NA)
probs_df_male <- probs_df_male %>% relocate(micro_specimen_id,.after = "NT") %>% 
  relocate(subject_id,.after="micro_specimen_id") %>% 
  relocate(I,.before="S") %>% mutate(id_no=1:nrow(probs_df_male)) %>% 
  relocate(id_no,.before="Antimicrobial")
write_csv(probs_df_male,"probs_df_male.csv")

util_probs_df <- read_csv("pre_util_probs_df.csv")
ur_util <- read_csv("ur_util.csv")
male_key <- ur_util %>% select(subject_id,MALE)
util_probs_df <- util_probs_df %>% left_join(male_key,by="subject_id")
male_rows <- which(util_probs_df$MALE==TRUE)
tzp_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Piperacillin-tazobactam"))
cefaz_rows <-intersect(male_rows,which(util_probs_df$Antimicrobial=="Cefazolin"))
ceftri_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Ceftriaxone"))
ceftaz_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Ceftazidime"))
fep_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Cefepime"))
mem_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Meropenem"))
gent_rows <- intersect(male_rows,which(util_probs_df$Antimicrobial=="Gentamicin"))
final_rows <- c(tzp_rows,cefaz_rows,ceftri_rows,ceftaz_rows,fep_rows,mem_rows,gent_rows)
util_probs_df$S[final_rows] <- probs_df_male$S[final_rows]
write_csv(util_probs_df,"male_pre_util_probs_df.csv")

