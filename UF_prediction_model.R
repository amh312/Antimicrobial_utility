#PREDICTION MODEL AND FEATURE TUNING

##Functions

###Search for previous event across multiple variables
apply_prev_event <- function(df, param,organism) {
  df %>%
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, 365, 1)
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  params <- paste0("p", antibiotics, suffix)
  
  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_event(df, params[i], antibiotics[i])
  }, .init = df1) %>%
    ungroup()
  
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  params <- paste0("p", antibiotics, suffix)
  
  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_event(df, params[i], antibiotics[i])
  }, .init = df1) %>%
    ungroup()
  
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(event_type, event)) %>%
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Assigning previous treatment
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx_groupvar <- enquo(abx_groupvar)
  
  drug_df %>%
    select('subject_id', ab_name,charttime='starttime') %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    bind_rows(ur_df) %>% 
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df,icd_df,prefix,codes,time_frame) {
  
  apply_prev_event_assignments <- function(df, code) {
    param_name <- paste0(prefix, code)
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, time_frame, 1)
  }
  
  df <- reduce(codes, function(df, code) {
    apply_prev_event_assignments(df, code)
  }, .init = df) %>%
    mutate(pDIAG_U = FALSE) %>%
    ungroup()
  
}

###Assigning previous event feature variable
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(!is.na(event)) %>% 
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Search for previous event across multiple variables
apply_prev_event <- function(df, param,organism,time_frame=365,no_events=1) {
  df %>%
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, time_frame, no_events)
}

###Checking for previous care events
care_event_assigner <- function(df,search_df,search_term,search_column,feature_name,event_date_col,timeframe,n_events=1) {
  
  feature_name <- enquo(feature_name)
  search_column <- enquo(search_column)
  
  care_event <- search_df %>% filter(grepl(search_term,!!search_column,ignore.case=T)) %>% mutate(
    !!search_column:=search_term) %>% rename(admittime=event_date_col)
  df %>% 
    prev_event_type_assign(!!feature_name,care_event,!!search_column,search_term,timeframe,n_events) %>%
    ungroup()
  
}

##Read-in

train_abx <- read_csv("train_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")

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

##Feature tuning

###Tuning time to last resistance
urines_ref <- read_csv("urines_ref.csv")
ur_util <- ur_xg_combined
micro3 <- micro %>% rename(admittime = "charttime")
time_list <- c(30,180,720,1e4)
for (i in 1:length(time_list)) {
  
  antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                   "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                   "CLI", "LVX", "AMK", "TOB")
  urines_ref <- prev_AST_applier(urines_ref,micro3,glue("r_{as.character(time_list[i])}"),"R",timeframe=time_list[i])
  ur_util <- prev_AST_applier(ur_util,micro3,glue("r_{as.character(time_list[i])}"),"R",timeframe=time_list[i])
  
}
new_r_cols <- urines_ref %>% select(pAMPr_30:pTOBr_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_r_cols))
time_list <- c(7,30,180,720,1e4)

###Tuning time to last susceptibility
for (i in 1:length(time_list)) {
  
  antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                   "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                   "CLI", "LVX", "AMK", "TOB")
  urines_ref <- prev_AST_applier(urines_ref,micro3,glue("s_{as.character(time_list[i])}"),"S",timeframe=time_list[i])
  ur_util <- prev_AST_applier(ur_util,micro3,glue("s_{as.character(time_list[i])}"),"S",timeframe=time_list[i])
  
}
new_s_cols <- urines_ref %>% select(pAMPs_30:pTOBs_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_s_cols))

###Tuning time to last treatment
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin", "Cefalexin", "Cefpodoxime proxetil",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- c("AMPrx", "AMXrx", "AMCrx", "SAMrx", "TZPrx", "CZOrx", "CZOrx", "CZOrx",
              "CROrx", "CAZrx", "FEPrx", "MEMrx", "ETPrx", "ATMrx", "CIPrx", "CIPrx",
              "GENrx", "TOBrx", "AMKrx", "RIFrx", "SXTrx", "NITrx", "ERYrx", "CLRrx",
              "AZMrx", "CLIrx", "VANrx", "MTRrx", "LNZrx", "DAPrx", "DOXrx")
apply_prev_rx <- function(df, suffix, antibiotic,time_to_event=365) {
  param_name <- paste0("p", suffix,"_",time_list[j])
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, time_to_event, 1)
}
time_list <- c(30,180,720,1e4)
drugs$ab_name <- drugs$abx_name
for (j in 1:length(time_list)) {
  urines_ref <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i],time_list[j])
  }, .init = urines_ref) %>%
    ungroup()
  ur_util <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_rx(df, suffixes[i], antibiotics[i],time_list[j])
  }, .init = ur_util) %>%
    ungroup()
}
new_rx_cols <- urines_ref %>% select(pAMPrx_30:pDOXrx_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_rx_cols))

###2 episodes of resistance in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                 "CLI", "LVX", "AMK", "TOB")
urines_ref <- prev_AST_applier(urines_ref,micro3,"r2","R",timeframe=365,n_events=2)
ur_util <- prev_AST_applier(ur_util,micro3,"r2","R",timeframe=365,n_events = 2)

new_2r_cols <- urines_ref %>% select(pAMPr2:pTOBr2)
urines5_combined <- data.frame(cbind(urines5_combined,new_2r_cols))

###2 antimicrobial courses in the last year
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin", "Cefalexin", "Cefpodoxime proxetil",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- c("AMPrx", "AMXrx", "AMCrx", "SAMrx", "TZPrx", "CZOrx", "CZOrx", "CZOrx",
              "CROrx", "CAZrx", "FEPrx", "MEMrx", "ETPrx", "ATMrx", "CIPrx", "CIPrx",
              "GENrx", "TOBrx", "AMKrx", "RIFrx", "SXTrx", "NITrx", "ERYrx", "CLRrx",
              "AZMrx", "CLIrx", "VANrx", "MTRrx", "LNZrx", "DAPrx", "DOXrx")
apply_prev_rx <- function(df, suffix, antibiotic,time_to_event=365,no_events=2) {
  param_name <- paste0("p", suffix,"2")
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, time_to_event, no_events)
}
urines_ref <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],365,2)
}, .init = urines_ref) %>%
  ungroup()
ur_util <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],365,2)
}, .init = ur_util) %>%
  ungroup()
new_rx_cols <- urines_ref %>% select(pAMPrx2:pDOXrx2)
urines5_combined <- data.frame(cbind(urines5_combined,new_rx_cols))

###Tuning hospital admission
time_list <- c(30,180,720,1e4)
for (i in 1:length(time_list)) {
  urines_ref <- urines_ref %>% 
    prev_event_assign(pHADM2,hadm,hadm_id,time_list[i],1) %>%
    ungroup()
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pHADM_{time_list[i]}")
  ur_util <- ur_util %>% 
    prev_event_assign(pHADM2,hadm,hadm_id,time_list[i],1) %>%
    ungroup()
  colnames(ur_util)[ncol(urines_ref)] <- glue("pHADM_{time_list[i]}")
}
new_hadm_cols <- urines_ref %>% select(pHADM_10000:pHADM_720)
urines5_combined <- data.frame(cbind(urines5_combined,new_hadm_cols))
urines_ref <- urines_ref %>% select(-pHADM_10000)
write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###2 hospital admissions in the last year
ur_util <- ur_util %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,365,2) %>%
  ungroup()
urines_ref <- urines_ref %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,365,2) %>%
  ungroup()
hadm2 <- urines_ref %>% select(pHADM2)
urines5_combined <- data.frame(cbind(urines5_combined,hadm2))

###Tuning coded ICD-10 diagnosis
diag_codes <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
                "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                "V", "W", "X", "Y", "Z")
for (i in 1:length(time_list)) {
  urines_ref <- urines_ref %>% prev_ICD_applier(diagnoses,glue("pDIAG{time_list[i]}_"),diag_codes,time_frame = time_list[i])
  ur_util <- ur_util %>% prev_ICD_applier(diagnoses,glue("pDIAG{time_list[i]}_"),diag_codes,time_frame = time_list[i])
}
new_icd_cols <- urines_ref %>% select(pDIAG30_A:pDIAG10000_Z)
urines5_combined <- data.frame(cbind(urines5_combined,new_icd_cols))

###Tuning coded ICD-10 procedure
proc_codes <- c("0", "3", "8", "5", "T", "4", "S", "A", "9", 
                "H", "I", "B", "7", "G", "1", "R", "J", "Q", 
                "K", "6", "M", "P", "L", "D", "F", "2", "N", 
                "C", "E", "X", "O")
for (i in 1:length(time_list)) {
  urines_ref <- urines_ref %>% prev_ICD_applier(procedures,glue("pPROC{time_list[i]}_"),proc_codes,time_frame = time_list[i])
  ur_util <- ur_util %>% prev_ICD_applier(procedures,glue("pPROC{time_list[i]}_"),proc_codes,time_frame = time_list[i])
}
new_proc_cols <- urines_ref %>% select(pPROC30_0:pPROC10000_O)
urines5_combined <- data.frame(cbind(urines5_combined,new_proc_cols))

###Tuning recent organism growth
urine_df <- urine_df %>% mutate(admittime=charttime)
organisms <- urine_df %>% count(org_fullname) %>% arrange(desc(n)) %>% 
  dplyr::slice(1:10) %>% pull(org_fullname)
for (j in 1:length(time_list)) {
  urines_ref <- reduce(seq_along(organisms), function(df, i) {
    apply_prev_event(df, paste0("pG", organisms[i],"Urine",as.character(time_list[j])), organisms[i],time_frame=time_list[j])
  }, .init = urines_ref) %>%
    ungroup()
  ur_util <- reduce(seq_along(organisms), function(df, i) {
    apply_prev_event(df, paste0("pG", organisms[i],"Urine",as.character(time_list[j])), organisms[i],time_frame=time_list[j])
  }, .init = ur_util) %>%
    ungroup()
}
new_org_cols <- urines_ref %>% select(`pGEscherichia coliUrine30`:`pGMorganella morganiiUrine10000`)
urines5_combined <- data.frame(cbind(urines5_combined,new_org_cols))

###Tuning other previous care events
time_list1 <- c(7,365,720,1e4)
time_list2 <- c(7,30,720,1e4)
for (i in 1:length(time_list1)) {
  
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"DNR",field_value,pDNR2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pDNR_{time_list2[i]}")
  urines_ref <- urines_ref %>%   
    care_event_assigner(poe,"Psychiatry",field_value,pPsych2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pPsych_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Nephrostomy",field_value,pNeph2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pNeph_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Surgery",field_value,pSURG2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pSURG_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pNUTR_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pPhysio_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Restraints",order_subtype,pRestr2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pRestr_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pOT_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Central TPN",order_subtype,pTPN2,"ordertime",time_list2[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pTPN_{time_list2[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"cath",field_value,pCATH2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pCATH_{time_list1[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Discharge",field_value,pDISC2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pDISC_{time_list1[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"ICU",field_value,pICU2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pICU_{time_list1[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"NGT",field_value,pNGT2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pNGT_{time_list1[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Hydration",field_value,pHyd2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pHyd_{time_list1[i]}")
  urines_ref <- urines_ref %>% 
    care_event_assigner(poe,"Chemo",field_value,pChemo2,"ordertime",time_list1[i])
  colnames(urines_ref)[ncol(urines_ref)] <- glue("pChemo_{time_list1[i]}")
  
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"DNR",field_value,pDNR2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pDNR_{time_list2[i]}")
  ur_util <- ur_util %>%   
    care_event_assigner(poe,"Psychiatry",field_value,pPsych2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pPsych_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Nephrostomy",field_value,pNeph2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pNeph_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Surgery",field_value,pSURG2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pSURG_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pNUTR_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pPhysio_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Restraints",order_subtype,pRestr2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pRestr_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pOT_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Central TPN",order_subtype,pTPN2,"ordertime",time_list2[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pTPN_{time_list2[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"cath",field_value,pCATH2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pCATH_{time_list1[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Discharge",field_value,pDISC2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pDISC_{time_list1[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"ICU",field_value,pICU2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pICU_{time_list1[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"NGT",field_value,pNGT2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pNGT_{time_list1[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Hydration",field_value,pHyd2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pHyd_{time_list1[i]}")
  ur_util <- ur_util %>% 
    care_event_assigner(poe,"Chemo",field_value,pChemo2,"ordertime",time_list1[i])
  colnames(ur_util)[ncol(ur_util)] <- glue("pChemo_{time_list1[i]}")
  
}
new_poe_cols <- urines_ref %>% select(pDNR_7:pChemo_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_poe_cols))

###Saving interim CSVs
write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")
write_csv(urines_ref,"urines_ref_combined.csv")


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

  

##Final clinical prediction model

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

###Saving interim CSVs
write_csv(micro_probs_df,"micro_xg_probs_df.csv")
write_csv(test_probs_df,"test_xg_probs_df.csv")
write_csv(aucs,"xg_aucs.csv")

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

##CDI prediction model

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
  eta = final_bestparams[[1]]$eta,
  max_depth = final_bestparams[[1]]$max_depth,
  min_child_weight = final_bestparams[[1]]$min_child_weight,
  subsample = final_bestparams[[1]]$subsample,
  colsample_bytree = final_bestparams[[1]]$colsample_bytree
)

print("Training...")

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = final_bestparams[[1]]$best_nrounds
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

###Feature selection
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(shap_summary %>% pull(Feature))), 
                            label = abxTrain[['CDI']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(shap_summary %>% pull(Feature))), 
                           label = abxTest[['CDI']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(shap_summary %>% pull(Feature))), 
                            label = ur_abx_combined[['CDI']])

###Run again with selected features
xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = final_bestparams[[1]]$best_nrounds,
)

pred_prob_test <- predict(xgb_model, newdata = test_matrix)
roc_result <- roc(abxTest[['CDI']], pred_prob_test)
auc_value <- auc(roc_result)
print(paste("AUC-ROC:", auc_value))
pdf(glue("CDI_xg_roc.pdf"), width = 10, height = 10)
plot(roc_result, main = glue("CDI ROC Curve"), col = "blue")
dev.off()
cdi_util_probs <- predict(xgb_model, newdata = micro_matrix)

pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
actual_test_class <- abxTest[['CDI']]

cdi_confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
cdi_accuracy <- confusion$overall['Accuracy']
cdi_precision <- confusion$byClass['Precision']
cdi_recall <- confusion$byClass['Recall']
cdi_f1_score <- 2 * (precision * recall) / (precision + recall)
probs_df_overall$prob_CDI <- cdi_util_probs
write_csv(probs_df_overall,"probs_df_overall.csv")

###Toxicity
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
  eta = final_bestparams[[1]]$eta,
  max_depth = final_bestparams[[1]]$max_depth,
  min_child_weight = final_bestparams[[1]]$min_child_weight,
  subsample = final_bestparams[[1]]$subsample,
  colsample_bytree = final_bestparams[[1]]$colsample_bytree
)

print("Training...")

xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = final_bestparams[[1]]$best_nrounds
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

###Feature selection
train_matrix <- xgb.DMatrix(data = as.matrix(abxTrain %>% select(shap_summary %>% pull(Feature))), 
                            label = abxTrain[['overall_tox']])
test_matrix <- xgb.DMatrix(data = as.matrix(abxTest %>% select(shap_summary %>% pull(Feature))), 
                           label = abxTest[['overall_tox']])
micro_matrix <- xgb.DMatrix(data = as.matrix(ur_abx_combined %>% select(shap_summary %>% pull(Feature))), 
                            label = ur_abx_combined[['overall_tox']])

###Run again with selected features
xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = final_bestparams[[1]]$best_nrounds,
)

pred_prob_test <- predict(xgb_model, newdata = test_matrix)
roc_result <- roc(abxTest[['overall_tox']], pred_prob_test)
auc_value <- auc(roc_result)
print(paste("AUC-ROC:", auc_value))
pdf(glue("overall_tox_xg_roc.pdf"), width = 10, height = 10)
plot(roc_result, main = glue("overall_tox ROC Curve"), col = "blue")
dev.off()
overall_tox_util_probs <- predict(xgb_model, newdata = micro_matrix)

pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
actual_test_class <- abxTest[['overall_tox']]

overall_tox_confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
overall_tox_accuracy <- confusion$overall['Accuracy']
overall_tox_precision <- confusion$byClass['Precision']
overall_tox_recall <- confusion$byClass['Recall']
overall_tox_f1_score <- 2 * (precision * recall) / (precision + recall)
probs_df_overall$prob_tox <- overall_tox_util_probs
write_csv(probs_df_overall,"probs_df_overall.csv")
