###PRESCRIPTION OUTCOME PREDICTION MODEL FINAL TRAINING AND VALIDATION

set.seed(123)

##Functions

  ###Factorise training and testing datasets
  factorise <- function(df) {
    df %>% mutate(CDI = factor(CDI),
                  overall_tox = factor(overall_tox),
                  sepsis_ae=factor(sepsis_ae))
  }

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
  
  ###Preprocessing of model and microsimulation dataframes
  abx_ur_datpreproc <- function(ab_df,ab_predname,ab_outcname,ab_combname,
                                abur_df,abur_predname,abur_outcname,abur_combname,
                                probs_df) {
    
    set.seed(123)
    
    #convert prescription df outcomes to 1s and 0s
    ab_outcomes <- ab_df %>%
      select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                         overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
    
    #select predictor variables of choice
    ab_predictors <- ab_df %>% select(pHADM:age65,prAKI:pDIAB,pCARD:curr_service,pICU:pSEPSIS,ob_freq,highCRP,
                                      pc_dyspnea:SIRS,
                                      ab_name_Ampicillin_Ceftriaxone:ab_name_Ampicillin,race:language,anchor_age)
    
    #make boolean dummy features from predictors
    ab_dummy <- dummyVars(" ~ .", data = ab_predictors)
    ab_predictors <- predict(ab_dummy, newdata = ab_predictors)
    
    #put outcomes and predictors back together
    ab_combined <- as.data.frame(cbind(ab_outcomes, ab_predictors))
    
    #convert microsim df outcomes to 1s and 0s
    ur_outcomes <- abur_df %>%
      select(micro_specimen_id,CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                                           overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
    
    #join microsim outcomes to antimicrobials in probs df
    ur_outcomes <- probs_df %>%
      left_join(ur_outcomes, relationship = 'many-to-one') %>% 
      select(-c(id_no,I:subject_id))
    
    #filter microsim to predictor columns based on those in prescriptions df
    common_columns <- intersect(names(ab_predictors),names(abur_df))
    ur_predictors <- abur_df %>% select(micro_specimen_id,all_of(common_columns))
    
    #join predictor vars to antimicrobials in probs df
    ur_predictors <- probs_df %>% left_join(ur_predictors,
                                            relationship = 'many-to-one') %>% 
      select(-c(id_no,I:subject_id)) %>% rename(ab_name="Antimicrobial")
    
    #make dummy features from microsim predictors
    dummy_ur <- dummyVars(" ~ .", data = ur_predictors)
    ur_predictors <- predict(dummy_ur, newdata = ur_predictors)
    
    #put microsim predictors and outcomes back together
    ur_combined <- as.data.frame(cbind(ur_outcomes, ur_predictors))
    
    #standardise naming of ab name feature variables
    colnames(ur_combined)[grepl("ab_name",colnames(ur_combined))] <- colnames(ur_combined)[grepl("ab_name",colnames(ur_combined))] %>% 
      str_replace("ab_name","ab_name_") %>% 
      str_replace_all("-",".")
    
    assign(ab_predname,ab_predictors,envir = .GlobalEnv)
    assign(ab_outcname,ab_outcomes,envir = .GlobalEnv)
    assign(ab_combname,ab_combined,envir = .GlobalEnv)
    assign(abur_predname,ur_predictors,envir = .GlobalEnv)
    assign(abur_outcname,ur_outcomes,envir = .GlobalEnv)
    assign(abur_combname,ur_combined,envir = .GlobalEnv)
    
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
    
    assign("selected_columns",selected_columns,envir=.GlobalEnv)
    
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
    ur_shapvals <- predict(model, newdata = trainmat, predcontrib = TRUE)
    
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
  probclassactual <- function(testdf,model,testmat,outc,
                              probnam,classnam,actnam) {
    
    #get predicted probabilities
    ur_predprobs <- predict(model, newdata = testmat)
    
    #get predicted classes (S >0.5, 5 ≤0.5)
    ur_predclass <- ifelse(ur_predprobs > 0.5, 1, 0)
    
    #set 1 as positive
    ur_predclass <- factor(ur_predclass,levels=c(0,1))
    
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
    
    #plot title
    plottit <- ifelse(abcombo_replace(outc,combined_antimicrobial_map)=="overall_tox",
      "Toxicity",abcombo_replace(outc,combined_antimicrobial_map))
    
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
      geom_errorbar(data=ur_calib_df,aes(ymin=grloci,ymax=grupci),col="#F8766D",width=0)+
      
      #theme and titles
      theme_minimal() +
      labs(x = "Mean predicted probability", y = "Actual proportion of positives",
           title = glue("{plottit}\ncalibration curve")) +
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
  
  ###Clean up, split, and write confidence interval tables for manuscript
  cidf_cleanup_write <- function(cidf,sing_writnam,comb_writnam){
    
    #sub-function to replace antibiotic names
    repcomb <- function(df){
      
      df %>% mutate(Model = case_when(
        Model %in% combined_antimicrobial_map ~ abcombo_replace(
          Model, combined_antimicrobial_map
        ),TRUE~Model
      ))
      
    }
    
    #add model name and arrange
    cidf$Model <- c(shortmap,"CDI","overall_tox")
    cidf <- cidf %>% relocate(Model,.before=1)
    
    #get full list of shortnames for outcomes
    ci_filter_list <- c(all_singles,"CDI","overall_tox")
    
    #ensure model name is character
    cidf <- cidf %>% mutate(Model=as.character(Model))
    
    #filter to main analysis and replace overall_tox model name
    ci_singles_table <- cidf %>% filter(Model %in% ci_filter_list) %>% 
      mutate(Model=case_when(grepl("tox",Model)~"Toxicity",TRUE~Model))
    
    #filter for antibiotic combination analysis
    ci_combos_table <- cidf %>% filter(!Model %in% ci_filter_list)
    
    #replace antibiotic shortnames with longnames
    ci_singles_table <- ci_singles_table %>% repcomb()
    ci_combos_table <- ci_combos_table %>% repcomb() %>%
      mutate(Model = str_replace_all(Model,"_"," & "))
    
    write_csv(ci_singles_table,sing_writnam)
    write_csv(ci_combos_table,comb_writnam)
    
  }
  
  ###Make grids of plot pdfs and save to pdfs
  pdf_gridmaker <- function(plot_list,nam_main,nam_other){
    
    plots_main <- c(plot_list[1:13],plot_list[38:39])
    plots_other <- plot_list[14:37]
    
    main_grid_plot <- plot_grid(plotlist = plots_main, ncol = 3)
    other_grid_plot <- plot_grid(plotlist = plots_other, ncol = 3)
    
    ggsave(glue(nam_main), plot = main_grid_plot, width = 30, height = 60,limitsize=F)
    ggsave(glue(nam_other), plot = other_grid_plot, width = 30, height = 60,limitsize=F)
    
  }
  
  ###Clean up, split, and write feature contributions table for manuscript
  shap_tabler <- function(probs_df,full_nam,sing_nam){
    
    other_outcome_map <- list("CDI","toxicity")
    names(other_outcome_map) <- c("CDI","Toxicity")
    ab_filter_list <- probs_df %>% distinct(Antimicrobial) %>% unlist()
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
        grepl("Alanine Aminotransferase",Feature)~str_replace(Feature,"Alanine Aminotransferase","Most recent Alanine Aminotransferase"),
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
        Feature = str_replace(Feature,"Benzylpenicillins","Benzylpenicillin susceptibility in the last year"),
        Feature = str_replace(Feature,"AMKs","Amikacin susceptibility in the last year"),
        Feature = str_replace(Feature,"AMK","Amikacin"),
        Feature = str_replace(Feature,"Benzylpenicillinnt", "Benzylpenicillin not tested in the last year"),
        Feature = str_replace(Feature,"Tobramycins", "Tobramycin susceptibility in the last year"),
        Feature = str_replace(Feature,"Tobramycini", "Tobramycin intermediate in the last year")
      ) %>% 
      mutate(Model=str_replace_all(Model,"_"," & "))
    
    write_csv(feat_table,full_nam)
    
    feat_abs <- c(c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                    "MEM","CIP","GEN","SXT","NIT","VAN") %>% ab_name() %>% str_replace_all("/","-"),
                  "CDI","Toxicity")
    feat_table_singles <- feat_table %>% filter(Model %in% feat_abs) %>% 
      mutate(Model = factor(Model, levels = unique(Model))) %>%
      group_by(Model) %>% arrange(desc(abs(`Shapley value`)),.by_group = T) %>% ungroup() %>% 
      mutate(`Shapley value`=case_when(`Shapley value`<0.001&`Shapley value`>-0.001 ~ as.character("<0.001"),
                                       TRUE~as.character(`Shapley value`)))
    write_csv(feat_table_singles,sing_nam)
    
  }
  
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

##Read-in

  train_abx <- read_csv("train_abx.csv")
  urines5 <- read_csv("urines5.csv")
  ur_xg <- read_csv("interim_ur_util.csv")
  urines5_combined <- read_csv("urines5_combined.csv")
  ur_xg_combined <- read_csv("ur_xg_combined.csv")
  probs_df_overall <- read_csv("probs_df_overall_postur.csv")

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
  
  
##Preprocessing
  
  ###Make full abx dataframe, factorise outcomes, and add missing demographics
  abx <- data.frame(rbind(train_abx,test_abx)) %>% factorise() %>% binddems()
  ur_xg <- ur_xg %>% left_join(age_key,by="subject_id")
  
  ###Preprocessing outcomes and predictors
  abx_ur_datpreproc(abx,"abx_predictors","abx_outcomes","abx_combined",
                    ur_xg,"ur_abx_predictors","ur_abx_outcomes","ur_abx_combined",
                    probs_df_overall)

##Model training

  ###Read in hyperparameters
  cdi_tox_final_bestparams <- hypparamreader("cdi_tox_final_params_",abx_outcome_map)
  
  ###Make blank microsim probabilities dataframe
  micro_probs_df2 <- data.frame(matrix(nrow=nrow(probs_df_overall),ncol=0))
  
  ###Iterate over prescription outcomes
  for (outcome in abx_outcome_map) {
    
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
      
      ###Retrieving predicted probabilities/classes and actual class
      abxTest %>% probclassactual(xgb_abxmodel,abxtest_matrix,outcome,'abx_predprobs',
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
        bootstrap_perfmets(ur_perf_mets,1000,outcome)
      
      ###Make probability predictions on microsimulation dataset
      pred_prob_micro <- predict(xgb_abxmodel, newdata = abxmicro_matrix)
      
      ###Add plots and values to lists
      ur_aucs[[outcome]] <- abx_auroc
      roc_plots[[outcome]] <- roc_plot_abx
      ur_calibplots[[outcome]] <- abx_calibplot
      micro_probs_df2[[outcome]] <- pred_prob_micro
      shap_ur_summary_tables[[outcome]] <- shap_abx_summary
      ur_metrics_list[[outcome]] <- abx_metrics
      confidence_biglist[[outcome]] <- abx_confints
      
    }
    
  }
  
  ###Feature importances to CSV
  shap_ur_summary_tables %>% shapwriter(overall_map)
  
  ###Performance metrics to CSV
  ur_metrics_list %>% metricwriter(overall_map)
  
  ###Confidence intervals to csvs
  confidence_biglist %>% ci_writer(overall_map,"final_ci_df.csv")
  ci_df <- read_csv("final_ci_df.csv")
  ci_df %>% cidf_cleanup_write("ci_singles_table.csv","ci_combos_table.csv")
  
  ###Add microsimulation probability predictions to overall probs df and save
  probs_df_overall$prob_CDI <- micro_probs_df2$CDI
  probs_df_overall$prob_tox <- micro_probs_df2$overall_tox
  write_csv(probs_df_overall,"probs_df_overall.csv")
  
  ###Save ROC and calibration plots to grid pdfs
  roc_plots %>% pdf_gridmaker("roc_plots_main.pdf","roc_plots_other.pdf")
  ur_calibplots %>% pdf_gridmaker("calibration_plots_main.pdf","calibration_plots_other.pdf")

  ###Feature contributions dataframe
  probs_df_overall %>% shap_tabler("feat_table.csv","feat_table_singles.csv")

  
