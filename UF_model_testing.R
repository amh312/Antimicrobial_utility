#MODEL ADDITIONAL TESTING

##Functions


##Read-in and mapping lists

###Read-in
train_abx <- read_csv("train_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
urines_ref <- read_csv("urines_ref.csv")
pats <- read_csv("patients.csv")

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

###Read-in chosen model parameters
file_names <- glue("final_params_{combined_antimicrobial_map}.csv")
final_bestparams <- lapply(file_names, function(file) {
  read_csv(file)
})
namelist <- c("eta","max_depth","min_child_weight","subsample","colsample_bytree",
              "best_nrounds")
for (i in 1:length(final_bestparams)) {
  names(final_bestparams[[i]])[3:8] <- namelist
}

###Model stability analysis

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
metrics_df4 <- data.frame(matrix(nrow=0,ncol=10))
for (key in 1:length(names(metrics_biglist4))) {
  for (i in 1:length(metrics_biglist4[[key]])) {
    
    for (j in 1:length(metrics_biglist4[[key]][[i]])) {
      
      for (k in 1:length(metrics_biglist4[[key]][[i]][[j]])) {
        metrics_biglist4[[1]][[1]][[1]]
        print(metrics_biglist4[[key]][[i]][[j]][[k]])
        
        results <- data.frame(
          Antimicrobial = names(metrics_biglist4)[key],
          Iteration = i,
          Train_year = names(metrics_biglist4[[key]][[i]][[j]]),
          AUC = metrics_biglist4[[key]][[i]][[j]][[k]]$AUC,
          Accuracy = metrics_biglist4[[key]][[i]][[j]][[k]]$Accuracy,
          Precision = metrics_biglist4[[key]][[i]][[j]][[k]]$Precision,
          Recall = metrics_biglist4[[key]][[i]][[j]][[k]]$Recall,
          F1_score = metrics_biglist4[[key]][[i]][[j]][[k]]$F1_Score,
          Train_support = metrics_biglist4[[key]][[i]][[j]][[k]]$Train_support,
          Test_support = metrics_biglist4[[key]][[i]][[j]][[k]]$Test_support)
        
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
rownames(metrics_df4) <- NULL

