###PRESCRIPTION OUTCOME PREDICTION MODEL FINAL TRAINING AND VALIDATION

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
         urdftrain,envir=.GlobalEnv)
  assign(chosnametest,
         urdftest,envir=.GlobalEnv)
  
}

###Line up dummy features between dataframes for class absence
lineup_features <- function(microsim_df,pred_df,train_df){
  
  #get column names of predictors
  predictor_columns <- colnames(pred_df)
  
  #get columns in common with training df
  selected_columns <- intersect(predictor_columns, colnames(train_df))
  
  #find out which columns are missing in microsim df
  missing_cols <- setdiff(selected_columns, colnames(microsim_df))
  
  #fill missing microsim df columns with 0
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

###Feature importance calculation
shapper <- function(trainmat,model,outc) {
  
  #update message
  print(glue("Checking feature importances for {outc}"))
  
  #get shap values
  ur_shapvals <- predict(xgb_abxmodel, newdata = trainmat, predcontrib = TRUE)
  
  #make dataframe from shap list
  shap_df <- as.data.frame(ur_shapvals)
  
  #remove redundant last column
  shap_df <- shap_df[, -ncol(shap_df)]
  
  #make new dataframe with mean shaps for each feature
  shap_ursmry <- data.frame(
    Feature = colnames(shap_df),
    meanshap = colMeans(shap_df)
  )
  
  #remove zero features
  shap_ursmry <- shap_ursmry %>% filter(meanshap!=0)
  shap_ursmry <- shap_ursmry[order(-shap_ursmry$meanshap), ]
  
  shap_ursmry
  
}

###Feature selection using SHAP
mat_feat_selector <- function(dataset,shapsum,outc) {
  
  #select only non-zero features left in the shap dataframe
  xgb.DMatrix(data = as.matrix(dataset %>% select(shapsum %>% pull(Feature))), 
              label = dataset[[outc]])
  
}

###Retrieving predicted probabilities/classes and actual class
probclassactual <- function(testdf,testmat,outc,
                            probnam,classnam,actnam) {
  
  #get predicted probabilities
  ur_predprobs <- predict(xgb_abxmodel, newdata = urtest_matrix)
  
  #get predicted classes (S >0.5, 5 ≤0.5)
  ur_predclass <- ifelse(ur_predprobs > 0.5, 1, 0)
  
  #set 1 as positive
  ur_predclass <- relevel(factor(ur_predclass), ref = "1")
  
  #get actual class
  ur_actclass <- testdf[[outc]]
  
  #set 1 as positive
  ur_actclass <- relevel(factor(ur_actclass), ref = "1")
  
  #assign to global environment
  assign(probnam,ur_predprobs,envir=.GlobalEnv)
  assign(classnam,ur_predclass,envir=.GlobalEnv)
  assign(actnam,ur_actclass,envir=.GlobalEnv)
  
}

###AUROC value and ROC curve
roc_maker <- function(actclass,predpr,outc,aurocnam,
                      abmap){
  
  #ger roc curve
  urroc <- roc(actclass, predpr,levels=c(0,1))
  
  #get auroc
  ur_auroc_value <- auc(urroc)
  
  #update message
  print(glue("Validation AUROC for {outc} = {round(ur_auroc_value,2)}"))
  
  #assign to global env
  assign(aurocnam,ur_auroc_value,envir=.GlobalEnv)
  
  #plot roc curve
  ggroc(urroc,color = "blue3") + 
    
    #titles
    ggtitle(glue("{abcombo_replace(outc,abmap)}\nROC curve"))+
    labs(x = "False positive rate", y = "True positive rate")+
    
    #zero effect line
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
                 color = "grey", linetype = "dashed")+
    
    #auc annotation
    geom_rect(aes(xmin = 0.1, xmax = 0.4, ymin = 0.1, ymax = 0.2), 
              fill = "white", alpha = 0.2, color = "grey9") +
    annotate("text", x = 0.25, y = 0.15, label = glue("AUC: {round(ur_auroc_value,2)}"),
             size = 10, color = "grey9") +
    
    #theme, text and ticks
    theme_minimal()+
    theme(
      plot.title = element_text(hjust = 0.5,size=30),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
      
    )
  
}

###Calibration curve and slope value
calibmaker <- function(actc,predp,outc){
  
  #predicted probs and actual class into dataframe
  ur_calib_df <- data.frame(
    pred_probs = predp,
    act_probs = as.numeric(as.character(actc))
  ) %>% 
    
    #put predicted probs into 10 bins
    mutate(probs_bin=cut(pred_probs,
                         breaks=quantile(pred_probs,
                                         probs=seq(0.1,1,by=0.1),
                                         na.rm=T),
                         labels=F)) %>% 
    
    #get means and n samples
    group_by(probs_bin) %>%
    summarise(meanpp = mean(pred_probs),
              act_prop = mean(act_probs),
              nsamp=n()) %>% 
    ungroup()
  
  #loess smoothed values for actual probabilities
  loesspreds <- predict(loess(ur_calib_df$act_prop~ur_calib_df$meanpp),
                        span=1,se=T)
  ur_calib_df$sm_act <- loesspreds$fit
  
  #smoothing 95% confidence intervals
  ur_calib_df$upperci <- loesspreds$fit+1.96*loesspreds$se.fit
  ur_calib_df$lowerci <- loesspreds$fit-1.96*loesspreds$se.fit
  
  #actual means 95% confidence intervals
  ur_calib_df$grupci <- ur_calib_df$act_prop+(sqrt((ur_calib_df$act_prop*1-ur_calib_df$act_prop)/ur_calib_df$nsamp)*1.96)
  ur_calib_df$grloci <- ur_calib_df$act_prop-((sqrt(ur_calib_df$act_prop*1-ur_calib_df$act_prop)/ur_calib_df$nsamp)*1.96)
  
  #slope from linear model
  urcalib_model <- lm(ur_calib_df$act_prop~ur_calib_df$meanpp)
  ur_calslope <- coef(urcalib_model)[2]
  
  #values for annotation box
  xrange <- range(ur_calib_df$meanpp)
  xrangesize <- xrange[2]-xrange[1]
  yrange <- range(ur_calib_df$act_prop)
  yrangesize <- yrange[2]-yrange[1]
  max_xbox <- xrange[1]+xrangesize*0.65
  min_xbox <- xrange[1]+xrangesize*0.95
  x_text <- mean(c(min_xbox,max_xbox))
  min_y <- yrange[1]+yrangesize*0.01
  max_y <- yrange[1]+yrangesize*0.12
  y_text <- mean(c(min_y,max_y))
  
  #calibration plot
  ggplot(ur_calib_df, aes(x = meanpp, y = sm_act)) +
    
    #loess smoothed line
    geom_line(color = "#00BFC4", linetype = "solid") +
    
    #ideal line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey20") +
    
    #loess confidence intervals
    geom_ribbon(data = ur_calib_df, aes(x = meanpp, ymin = lowerci, ymax = upperci), fill = "#00BFC4", alpha = 0.2)+
    
    #actual means
    geom_point(data=ur_calib_df,aes(x=meanpp,y=act_prop),col="#F8766D",size=3)+
    
    #confidence intervals of means
    geom_errorbar(data=ur_calib_df,aes(ymin=grloci,ymax=grupci),col="#F8766D",width=0.005)+
    
    #theme and titles
    theme_minimal() +
    labs(x = "Mean predicted probability", y = "Actual proportion of positives",
         title = glue("{abcombo_replace(outc,combined_antimicrobial_map)}\ncalibration curve")) +
    theme(plot.title = element_text(hjust = 0.5,size = 30),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))+
    
    #x and y limits
    ylim(min(ur_calib_df$lowerci),max(ur_calib_df$upperci))+
    xlim(min(ur_calib_df$meanpp),max(ur_calib_df$meanpp))+
    
    #annotation with slope
    geom_rect(aes(xmin = min_xbox, xmax = max_xbox, ymin = min_y, ymax = max_y), 
              fill = "white", alpha = 0.2, color = "grey9") +
    annotate("text", x = x_text, y = y_text, label = glue("Slope: {round(ur_calslope,2)}"),
             size = 10, color = "grey9")
  
  
}

###Classification report
ur_perf_mets <- function(df, indexrows,bootstr=T) {
  
  #if bootstrapping
  if (bootstr==T) {
    
    #subset by selected indices
    ur_act <- df$act_val[indexrows]
    ur_probs <- df$pred_probs[indexrows]
    ur_class <- df$pred_class[indexrows]
    
    #confusion matrix
    ur_confmat <- confusionMatrix(factor(ur_class), factor(ur_act))
    
    #accuracy
    acc <- ur_confmat$overall['Accuracy']
    
    #precision
    prec <- ur_confmat$byClass['Precision']
    
    #recall
    rec <- ur_confmat$byClass['Recall']
    
    #f1 score
    f1 <- 2 * (prec * rec) / (prec + rec)
    
    #auroc
    auroc <- auc(roc(ur_act, ur_probs,levels=c(0,1)))
    
    #return vector
    c(auroc = auroc, precision = prec, recall = rec, accuracy = acc, f1 = f1)
    
    #if not bootstrapping
  } else {
    
    #use whole vectors
    ur_act <- df$act_val
    ur_probs <- df$pred_probs
    ur_class <- df$pred_class
    
    #confusion matrix
    ur_confmat <- confusionMatrix(factor(ur_class), factor(ur_act))
    
    #accuracy
    acc <- ur_confmat$overall['Accuracy']
    
    #precision
    prec <- ur_confmat$byClass['Precision']
    
    #recall
    rec <- ur_confmat$byClass['Recall']
    
    #f1 score
    f1 <- 2 * (prec * rec) / (prec + rec)
    
    #auroc
    auroc <- auc(roc(ur_act, ur_probs,levels=c(0,1)))
    
    #return list
    list(
      AUC = auroc,
      Accuracy = acc,
      Precision = prec,
      Recall = rec,
      F1_Score = f1
    )
    
  }
  
}

###Bootstrapping of performance metrics for CIs
bootstrap_perfmets <- function(pred_df,classrep_func,n_samps,outc) {
  
  #update message
  print(glue("Bootstrapping for {outc}"))
  
  #set seed
  set.seed(123)
  
  #get classification reports on 1000 bootstrapped samples
  ur_bootmetrics <- boot(data = pred_df, statistic = classrep_func, R = n_samps)
  
  #replace any NAs with 0
  ur_bootmetrics[[1]][is.na(ur_bootmetrics[[1]])] <- 0
  ur_bootmetrics[[2]][is.na(ur_bootmetrics[[2]])] <- 0
  
  #get 95% confidence intervals for each metric
  ur_cis <- lapply(1:5, function(i) boot.ci(ur_bootmetrics, type = "perc", index = i))
  
  #name metrics in dataframe and return
  names(ur_cis) <- c("AUC", "Precision", "Recall", "Accuracy", "F1_Score")
  ur_cis
  
}

###Writing feature importances to csv
shapwriter <- function(shaptable,abmap) {
  
  #iterate over shap tables list
  for (i in 1:length(shaptable)) {
    
    #make dataframe from list and add names
    shappy <- data.frame(shaptable[i])
    colnames(shappy) <- c("Feature","Shapley value")
    
    shappy <- shappy %>% 
      
      #round mean shapley value to 3 decimal places
      mutate(`Shapley value`=round(`Shapley value`,3)) %>% 
      
      #filter to only non-zero values
      filter(`Shapley value`!= 0)
    
    
    #write to csv
    write_csv(shappy,glue("SHAP_{abmap[i]}.csv"))
    
  }
  
}

###Writing performance metrics to CSV
metricwriter <- function(metlist,abmap) {
  
  #iterate over antibiotics
  for (i in 1:length(metlist)) {
    
    #make df from list
    metricky <- data.frame(metlist[i])
    
    #write to csv
    write_csv(metricky,glue("metrics_{abmap[i]}.csv"))
    
  }
  
}

###Writing confidence intervals on performance metrics to csv
ci_writer <- function(conflist,outcomedf,filename){
  
  #empty dataframe
  ci_df <- data.frame(matrix(nrow=length(outcomedf),ncol=5))
  colnames(ci_df) <- c("AUROC_CI","Precision_CI","Recall_CI","F1_CI","Accuracy_CI")
  
  #iterate over antibiotics
  for (i in 1: length(conflist)) {
    
    #paste AUC confidence intervals into brackets
    ci_df$AUROC_CI[i] <- ifelse(
      !is.null(conflist[[i]]$AUC$percent[4])&!is.null(conflist[[i]]$AUC$percent[5]),
      glue(" ({round(conflist[[i]]$AUC$percent[4],3)}-{round(conflist[[i]]$AUC$percent[5],3)})"),
      NA)
    
    #paste precision confidence intervals into brackets
    ci_df$Precision_CI[i] <- ifelse(
      !is.null(conflist[[i]]$Precision$percent[4])&!is_null(conflist[[i]]$Precision$percent[5]),
      glue(" ({round(conflist[[i]]$Precision$percent[4],3)}-{round(conflist[[i]]$Precision$percent[5],3)})"),
      NA)
    
    #paste recall confidence intervals into brackets
    ci_df$Recall_CI[i] <- ifelse(
      !is.null(conflist[[i]]$Recall$percent[4])&!is_null(conflist[[i]]$Recall$percent[5]),
      glue(" ({round(conflist[[i]]$Recall$percent[4],3)}-{round(conflist[[i]]$Recall$percent[5],3)})"),
      NA)
    
    #paste accuracy confidence intervals into brackets
    ci_df$Accuracy_CI[i] <- ifelse(
      !is.null(conflist[[i]]$Accuracy$percent[4])&!is_null(conflist[[i]]$Accuracy$percent[5]),
      glue(" ({round(conflist[[i]]$Accuracy$percent[4],3)}-{round(conflist[[i]]$Accuracy$percent[5],3)})"),
      NA)
    
    #paste F1 score confidence intervals into brackets
    ci_df$F1_CI[i] <- ifelse(
      !is.null(conflist[[i]]$F1_Score$percent[4])&!is_null(conflist[[i]]$F1_Score$percent[5]),
      glue(" ({round(conflist[[i]]$F1_Score$percent[4],3)}-{round(conflist[[i]]$F1_Score$percent[5],3)})"),
    )
  }
  
  #write to csv
  write_csv(ci_df,filename)
  
}

##Read-in

train_abx <- read_csv("train_abx.csv")
urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("interim_ur_util.csv")
urines5_combined <- read_csv("urines5_combined.csv")
ur_xg_combined <- read_csv("ur_xg_combined.csv")
abx_combined <- read_csv("abx_combined.csv")
ur_abx_combined <- read_csv("ur_abx_combined.csv")

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

###Outcome mapping list
abx_outcome_map <- list("CDI","overall_tox")
names(abx_outcome_map) <- c("C. difficile infection","Antibiotic toxicity")

###Big combined map
overall_map <- append(combined_antimicrobial_map,abx_outcome_map)

##Model training

###Read in hyperparameters
cdi_tox_final_bestparams <- hypparamreader("cdi_tox_final_params_",colnames(abx_outcomes))


###Iterate over prescription outcomes
for (outcome in colnames(abx_outcomes)) {
  
  ###Check no NAs
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    ###Set seed
    set.seed(123)
    
    ###Split df to train and test
    abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
    
    ###Ensure features line up between dataframes
    ur_abx_combined <- ur_abx_combined %>%
      lineup_features(abx_predictors,abxTrain)
    
    ###Make xgboost training matrices
    abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
    abxtest_matrix <- abxTest %>% model_matrixmaker(abx_predictors,outcome)
    abxmicro_matrix <- ur_abx_combined %>% model_matrixmaker(abx_predictors,outcome)
    
    ###Set parameters
    abxparams <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = cdi_tox_final_bestparams[[outcome]]$eta,
      max_depth = cdi_tox_final_bestparams[[outcome]]$max_depth,
      min_child_weight = cdi_tox_final_bestparams[[outcome]]$min_child_weight,
      subsample = cdi_tox_final_bestparams[[outcome]]$subsample,
      colsample_bytree = cdi_tox_final_bestparams[[outcome]]$colsample_bytree
    )
    
    ###First training with all features
    print(glue("First training for {outcome}"))
    xgb_abxmodel <- xgb.train(
      params = abxparams,data = abxtrain_matrix,
      nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
    )
    
    ###Feature importances
    shap_abx_summary <- abxtrain_matrix %>% shapper(xgb_abxmodel,outcome)
    
    ###Feature selection
    abxtrain_matrix <- abxTrain %>% mat_feat_selector(shap_abx_summary,outcome)
    abxtest_matrix <- abxTest %>% mat_feat_selector(shap_abx_summary,outcome)
    abxmicro_matrix <- ur_abx_combined %>% mat_feat_selector(shap_abx_summary,outcome)
    
    ###Second training with only selected features
    print(glue("Second training for {outcome}"))
    xgb_abxmodel <- xgb.train(
      params = abxparams,data = abxtrain_matrix,
      nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
    )
    
    ##Model validation
    
    ###Get predicted probability/class and actual class
    abxtestmatrix %>% probclassactual(abxTest,outcome,'abx_predprobs',
                                     'abx_predclass','abx_actualclass')
    
    ###ROC curve and AUROC
    roc_plot_abx <- roc_maker(abx_actualclass,abx_predprobs,outcome,
                             "abx_auroc",abx_outcome_map)
    
    ###Calibration plot and slope value
    abx_calibplot <- calibmaker(abx_actualclass,abx_predprobs,outcome)
    
    ###Save plots to file
    ggsave(filename = glue("{outcome}_xg_roc.pdf"), plot = roc_plot_abx,
           width = 10, height = 10)
    ggsave(filename = glue("{outcome}_calib_plot.pdf"), plot = abx_calibplot, width = 10, height = 10)
    
    ###Other performance metrics
    abx_predact_df <- data.frame(act_val = abx_actualclass, pred_class = abx_predclass,pred_probs = abx_predprobs)
    abx_metrics <- ur_perf_mets(abx_predact_df,NULL,bootstr = F)
    
    ###Bootstrapping for confidence intervals
    abx_confints <- abx_predact_df %>%
      bootstrap_perfmets(abx_perf_mets,1000,outcome)
    
    ###Make probability predictions on microsimulation dataset
    pred_prob_micro <- predict(xgb_abxmodel, newdata = abxmicro_matrix)
    
    ###Add plots and values to lists
    ur_aucs[[outcome]] <- abx_auroc
    roc_plots[[outcome]] <- roc_plot_abx
    ur_calibplots[[outcome]] <- abx_calibplot
    test_probs_df[[outcome]] <- abx_predprobs
    micro_probs_df[[outcome]] <- pred_prob_micro
    shap_ur_summary_tables[[outcome]] <- shap_abx_summary
    ur_metrics_list[[outcome]] <- abx_metrics
    confidence_biglist[[outcome]] <- abx_confints
    
  }
  
}

###Feature importances to CSV
shap_ur_summary_tables %>% shapwriter(overall_map)

###Performance metrics to CSV
ur_metrics_list %>% metricwriter(overall_map)

###Confidence intervals to csv
confidence_biglist %>% ci_writer(overall_map,"final_ci_df.csv")

###Add microsimulation probability predictions to overall probs df and save
probs_df_overall$prob_CDI <- micro_probs_df$CDI
probs_df_overall$prob_tox <- micro_probs_df$overall_tox
write_csv(probs_df_overall,"probs_df_overall.csv")

