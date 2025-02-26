###COVERAGE MODEL HYPERPARAMETER TUNING - STAGE 1

set.seed(123)

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
reverse_values <- function(column, map) {
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% map) names(map)[map == x] else x)
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
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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

for (outcome in colnames(urines5_outcomes)) {
  
  best_auc <- 0
  
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



