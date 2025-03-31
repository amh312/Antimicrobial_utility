#MODEL ADDITIONAL TESTING

set.seed(123)

##Functions

  ###Factorise training and testing datasets
  factorise <- function(df) {
    
    #convert cdi and toxicity to factors
    df %>% mutate(CDI = factor(CDI),
                  overall_tox = factor(overall_tox),
                  sepsis_ae=factor(sepsis_ae))
  }
  
  ###Time sens plot
  stability_plot <- function(df,metric,perf_metric) {
    
    #quosure
    metric <- enquo(metric)
    
    df <- df %>% 
      
      #filter to non-nas
      filter(!is.na(!!metric)) %>% 
      
      #clean prescription model names
      mutate(Model=abcombo_replace(Model,overall_map))
    
    #clean and add missing names
    df <- df %>% rename(`Training dataset size`="Training_size")
    model_levels <- names(overall_map) %>% rev()
    
    #factorise and characterise
    df$Model <- factor(df$Model,levels=model_levels)
    df$`Training dataset size` <- as.character(df$`Training dataset size`)
    
    #plot
    dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=`Training dataset size`,color=`Training dataset size`)) +
      
      #points
      geom_point()+
      
      #theme and titles
      theme_minimal()+
      ggtitle(glue("{perf_metric} after 50 XGBoost training rounds using different\ntraining dataset proportions"))+
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #limits
      xlim(0,1)
    
    #save and print
    ggsave(glue("stability_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(dfplot)
    
  }
  
  ###Fairness plot
  protchar_plot <- function(df,prot_char,metric,perf_metric,title_bit) {
    
    #quosure
    metric <- enquo(metric)
    
    #filter to relevant characteristics
    df <- df %>% filter(grepl(prot_char,Category)) %>% 
      
      #remove nas
      filter(!is.na(!!metric)) %>%
      
      #clean and add names
      mutate(Model=abcombo_replace(Model,overall_map))
      model_levels <- names(overall_map) %>% rev()
    
    #factors/characters
    df$Model <- factor(df$Model,levels=model_levels)
    df$Characteristic <- factor(df$Characteristic,
                                levels=df %>% filter(grepl(prot_char,Category)) %>% 
                                  distinct(Characteristic) %>% unlist())
    
    #plot
    dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=Characteristic,color=Characteristic)) +
      
      #points
      geom_point()+
      
      #theme and titles
      theme_minimal()+
      ggtitle(glue("{perf_metric} after 50 XGBoost training rounds for different\n{title_bit}"))+
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      xlim(0,1)
    
    #save and print
    ggsave(glue("protchar_{prot_char}_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(dfplot)
    
  }
  
  ###Time plot
  timesens_plot <- function(df,tr_yr,metric,perf_metric) {
    
    #quosures
    metric <- enquo(metric)
    
    #filter to training years, remove nas, clean and add names
    df <- df %>% filter(grepl(tr_yr,Train_year)) %>% 
      filter(!is.na(!!metric)) %>%
      mutate(Model=abcombo_replace(Model,overall_map))
      model_levels <- names(overall_map) %>% rev()
    
    #characters and factors
    df$Model <- factor(df$Model,levels=model_levels)
    df$Test_year <- factor(df$Test_year,
                           levels=df %>% filter(grepl(tr_yr,Train_year)) %>% 
                             distinct(Test_year) %>% unlist())
    
    #rename test year for plot
    df <- df %>% rename(`Test year range`="Test_year")
    
    #plot
    dfplot <- ggplot(df, aes(x=!!metric,y=Model,group=`Test year range`,color=`Test year range`)) +
      
      #points
      geom_point()+
      
      #theme and titles
      theme_minimal()+
      ggtitle(glue("{perf_metric} after 50 XGBoost training rounds for different\ntesting timeframes when trained on {tr_yr}"))+
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      
      #limits
      xlim(0,1)
    
    #save and print
    ggsave(glue("timesens_{tr_yr}_{perf_metric}.pdf"), plot = dfplot, device = "pdf", width = 12, height = 8,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(dfplot)
    
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
    
    abcombos_map <- combn(names(ab_map), 2, simplify = FALSE)
    c(antimicrobial_map,
      setNames(
        lapply(abcombos_map, function(x) paste(ab_map[x], collapse = "_")),
        sapply(abcombos_map, function(x) paste(x, collapse = "_"))
      )
    )
    
  }
  
  ###Make filtered, full and short antimicrobial maps
  abcombo_variants <- function(abmap,fullname,mainname,shortname,abreflist) {
    
    fullmap <- abmap %>% unlist()
    names(fullmap) <- NULL
    abmap2 <- abmap[names(abmap) %in% abreflist]
    shortmap <- abmap2 %>% unlist()
    names(shortmap) <- NULL
    
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
  
  ###Put metrics into dataframes
  metric_cleaner <- function(metric){
    
    #empty df
    metricdf <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
    
    #add antibiotics
    metricdf$Antimicrobial <- metrics_ablist
    
    #read in metrics csvs into df
    for (i in 1:length(metrics_ablist)) {
      
      metricdf[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains(metric))
      
    }
    
    #add names
    colnames(metricdf) <- c("Antimicrobial",metric)
    metricdf$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
    
    #add prescription models
    metricdf[38,1] <- "CDI"
    metricdf[39,1] <- "toxicity"
    
    #characterise
    metricdf$Antimicrobial <- as.character(metricdf$Antimicrobial)
    
    metricdf
    
  }
  
  ###Compile metric table including confidence intervals
  metric_tablecompiler <- function(aucdf,precdf,recalldf,f1df,accdf,
                                   singcitab,combcitab,outputsingtab,
                                   outputcombotab){
    
    #sub-function for cleaning table
    tablecleaner <- function(df,cidf){
      
      #put metrics and confidence intervals together
      df %>% mutate(Model=case_when(Model=="toxicity"~"Toxicity",TRUE~Model)) %>%
        left_join(cidf) %>% 
        mutate(
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
      
    }
    
    #join metrics together
    metrics_manus_table <- aucdf %>% left_join(precdf,by="Antimicrobial") %>% 
      left_join(recalldf,by="Antimicrobial") %>%
      left_join(f1df,by="Antimicrobial") %>%
      left_join(accdf,by="Antimicrobial") %>% 
      
      #rename columns
      rename(Model="Antimicrobial",AUROC="AUC",`F1 score`="F1") %>% 
      
      #rounding and clean combo names
      mutate(across(AUROC:`Accuracy`, ~ round(.x,3)),
             Model=str_replace_all(Model,"_"," & "))
    
    #read in ci dfs
    ci_singles_table <- read_csv(singcitab)
    ci_combos_table <- read_csv(combcitab)
    
    #single agent table cleaning
    metrics_singles_table <- metrics_manus_table %>% 
      dplyr::slice(c(1:13,38:39)) %>% tablecleaner(ci_singles_table)
    
    #combo agent table cleaning
    metrics_combos_table <- metrics_manus_table %>% 
      dplyr::slice(c(14:39)) %>% tablecleaner(ci_combos_table)
    
    #write to csv
    write_csv(metrics_singles_table,outputsingtab)
    write_csv(metrics_combos_table,outputcombotab)
    
  }
  
  ###Preprocessing of urine model and microsimulation dataframes
  ur_datpreproc <- function(ur_df,ur_outcomename,ur_predictorname,ur_comboname,
                            microsim_df,microsim_outcomename,microsim_predictorname,
                            microsim_comboname) {
    
    #sub-function for filling in marital status NAs
    mutmar <- function(dfx) {
      
      dfx %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                              TRUE~marital_status))
      
    }
    
    #sub-function to select and binarising ast results to 1s and 0s
    binariseast <- function(dfx) {
      
      dfx %>% select(all_of(shortmap)) %>%
        mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                       ifelse(. == "S" | . == "I", 1, NA))))
      
    }
    
    #sub-function to make dummy variables
    dummyer <- function(dfx) {
      
      urfeatdummies <- dummyVars(" ~ .", data = dfx)
      predict(urfeatdummies, newdata = dfx)
      
    }
    
    #fill marital status and select/binarise ast results
    ur_df_outcomes <- ur_df %>% mutmar() %>% binariseast()
    microsim_df_outcomes <- microsim_df %>% mutmar() %>% binariseast()
    
    #make predictor dataframes
    ur_df_predictors <- ur_df %>% select(!all_of(fullmap))
    microsim_df_predictors <- microsim_df %>%
      select(any_of(colnames(ur_df_predictors)))
    
    #make dummy variables
    ur_df_predictors <- ur_df_predictors %>% dummyer()
    microsim_df_predictors <- microsim_df_predictors %>% dummyer()
    
    #make combined dataframes
    ur_df_combined <- as.data.frame(cbind(ur_df_outcomes, ur_df_predictors))
    microsim_df_combined <- as.data.frame(cbind(microsim_df_outcomes, microsim_df_predictors))
    
    #assign resulting dfs to global environment
    assign(ur_outcomename,ur_df_outcomes,.GlobalEnv)
    assign(ur_predictorname,ur_df_predictors,.GlobalEnv)
    assign(ur_comboname,ur_df_combined,.GlobalEnv)
    assign(microsim_outcomename,microsim_df_outcomes,.GlobalEnv)
    assign(microsim_predictorname,microsim_df_predictors,.GlobalEnv)
    assign(microsim_comboname,microsim_df_combined,.GlobalEnv)
    
    
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
  
  ###Stability metrics to data frame
  stabmets_writer <- function(big_list){
    
    #empty df
    metdf <- data.frame(matrix(nrow=0,ncol=7))
    
    #iteration top level
    for (key in 1:length(names(big_list))) {
      
      #iteration level 3
      for (i in 1:length(big_list[[key]])) {
        
        #iteration level 2
        for (j in 1:length(big_list[[key]][[i]])) {
          
          #iteration level 1
          for (k in 1:length(big_list[[key]][[i]][[j]])) {
            
            #print current iteration level
            print(big_list[[key]][[i]][[j]][[k]])
            
            #add metrics to dataframe
            stabresults <- data.frame(
              Antimicrobial = names(big_list)[key],
              Training_size = names(big_list[[key]][[i]][1]),
              AUC = big_list[[key]][[i]][[j]][[k]]$AUC,
              Accuracy = big_list[[key]][[i]][[j]][[k]]$Accuracy,
              Precision = big_list[[key]][[i]][[j]][[k]]$Precision,
              Recall = big_list[[key]][[i]][[j]][[k]]$Recall,
              F1_score = big_list[[key]][[i]][[j]][[k]]$F1_Score)
            
            colnames(metdf) <- colnames(stabresults)
            
            #bind rows to large dataframe
            metdf <- data.frame(
              rbind(
                metdf,stabresults
              )
            )
            
          }
        }
      }
    }
    rownames(metdf) <- NULL
    metdf
    
  }
  
  ###Protected characteristic dataframe assembler
  protchar_assembler <- function(comb_df){
    
    comb_df %>% mutate(
      
      #non-english speaker
      language_NONENG = case_when(languageENGLISH!=1~1,TRUE~0),
      
      #race types
      race_ASIAN = rowSums(select(., contains("raceASIAN"))) > 0,
      race_BLACK = rowSums(select(., contains("raceBLACK"))) > 0,
      race_HISPANIC = rowSums(select(., contains("raceHISPANIC"))) > 0,
      race_OTHER = rowSums(select(., matches("(raceOTHER|racePORTUGUESE|raceSOUTH|raceNATIVE|raceAMERICAN|raceMULTIPLE)"))) > 0,
      race_WHITE = rowSums(select(., contains("raceWHITE"))) > 0,
      
      #marital status types
      marital_status_MARRIED = rowSums(select(., matches("(marital_statusMARRIED)"))) > 0,
      marital_status_NONMARRIED = case_when(marital_statusMARRIED!=1~1,TRUE~0)
    )
    
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
      ur_confmat <- confusionMatrix(factor(ur_class,levels=c(1,0)), factor(ur_act,levels=c(1,0)))
      
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
  
  ###Fairness dataframe maker
  protchar_writer <- function(big_list){
    
    #empty df
    metdf <- data.frame(matrix(nrow=0,ncol=9))
    
    #iteration top level
    for (key in 1:length(names(big_list))) {
      
      #iteration level 3
      for (i in 1:length(big_list[[key]])) {
        
        #iteration level 2
        for (j in 1:length(big_list[[key]][[i]])) {
          
          #iteration level 1
          for (k in 1:length(big_list[[key]][[i]][[j]])) {
            
            #print iteration level
            print(big_list[[key]][[i]][[j]][[k]])
            
            #add metrics to df
            results <- data.frame(
              Antimicrobial = names(big_list)[[key]],
              Iteration = i,
              Characteristic = big_list[[key]][[i]][[j]][[k]]$Characteristic,
              AUC = big_list[[key]][[i]][[j]][[k]]$AUC,
              Accuracy = big_list[[key]][[i]][[j]][[k]]$Accuracy,
              Precision = big_list[[key]][[i]][[j]][[k]]$Precision,
              Recall = big_list[[key]][[i]][[j]][[k]]$Recall,
              F1_score = big_list[[key]][[i]][[j]][[k]]$F1_Score,
              Test_support = big_list[[key]][[i]][[j]][[k]]$Test_support)
            
            #sync up colnames
            colnames(metdf) <- colnames(results)
            
            #add df rows to large df
            metdf <- data.frame(
              rbind(
                metdf,results
              )
            )
            
          }
        }
      }
    }
    rownames(metdf) <- NULL
    metdf
    
  }
  
  ###Combined age with other fairness metrics into dataframe
  agefairmettodf <- function(big_list,fairmetdf) {
    
    #empty df
    metdf <- data.frame(matrix(nrow=0,ncol=9))
    
    #iteration top level
    for (key in 1:length(names(big_list))) {
      
      #iteration level 3
      for (i in 1:length(big_list[[key]])) {
        
        #iteration level 2
        for (j in 1:length(big_list[[key]][[i]])) {
          
          #iteration level 1
          for (k in 1:length(big_list[[key]][[i]][[j]])) {
            
            #print iteration
            print(big_list[[key]][[i]][[j]][[k]])
            
            #add metrics to df
            results <- data.frame(
              Antimicrobial = names(big_list)[key],
              Iteration = i,
              Characteristic = big_list[[key]][[i]][[j]][[k]]$Characteristic,
              AUC = big_list[[key]][[i]][[j]][[k]]$AUC,
              Accuracy = big_list[[key]][[i]][[j]][[k]]$Accuracy,
              Precision = big_list[[key]][[i]][[j]][[k]]$Precision,
              Recall = big_list[[key]][[i]][[j]][[k]]$Recall,
              F1_score = big_list[[key]][[i]][[j]][[k]]$F1_Score,
              Test_support = big_list[[key]][[i]][[j]][[k]]$Test_support)
            
            #sync colnames
            colnames(metdf) <- colnames(results)
            
            #add df to large df
            metdf <- data.frame(
              rbind(
                metdf,results
              )
            )
            
          }
        }
      }
    }
    
    rownames(metdf) <- NULL
    
    #bind metrics dfs together
    data.frame(rbind(fairmetdf,metdf)) %>%
      
      #clean up names for table
      mutate(Category = case_when(
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
      )) %>% relocate(Category,.before = "Characteristic") %>%
      
      #arrange by levels
      arrange(Antimicrobial, Iteration, Category, Characteristic)
    
  }
  
  ###Time sensitivity metrics to dataframe
  timesenstodf <- function(big_list){
    
    #empty df
    metdf4 <- data.frame(matrix(nrow=0,ncol=11))
    
    #iteration top level
    for (key in 1:length(names(metrics_biglist4))) {
      
      #iteration level 4
      for (i in 1:length(metrics_biglist4[[key]])) {
        
        #iteration level 3
        for (j in 1:length(metrics_biglist4[[key]][[i]])) {
          
          #iteration level 2
          for (k in 1:length(metrics_biglist4[[key]][[i]][[j]])) {
            
            #iteration level 1
            for (l in 1:length(metrics_biglist4[[key]][[i]][[j]][[k]])) {
              
              #add results to df
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
              
              #sync colnames
              colnames(metdf4) <- colnames(results)
              
              #add df to larger df
              metdf4 <- data.frame(
                rbind(
                  metdf4,results
                )
              )
            }
          }
        }
      }
    }
    rownames(metdf4) <- NULL
    metdf4
    
  }
  
  ###Prescription data preprocessing
  modeltest_preproc <- function(df,outcnam,prednam,combnam){
    
    #set seed
    set.seed(123)
    
    #admissions key
    hadm_key <- hadm %>% select(subject_id,race,language,marital_status) %>% 
      distinct(subject_id,.keep_all=T)
    
    #age key
    age_key <- pats %>% select(subject_id,anchor_age) %>% mutate(age_grp=case_when(
      anchor_age<30~18,anchor_age>=30&anchor_age<40~30,
      anchor_age>=40&anchor_age<50~40,
      anchor_age>=50&anchor_age<60~50,
      anchor_age>=60&anchor_age<70~60,
      anchor_age>=70&anchor_age<80~70,
      anchor_age>=80&anchor_age<90~80,
      anchor_age>=90~90
    )) %>% select(-anchor_age) %>% distinct(subject_id,.keep_all = T)
    
    #bind age and admissions into demographics key
    dems_key <- left_join(hadm_key,age_key,by="subject_id")
    
    #join demographics to df
    df <- df %>% left_join(dems_key,by="subject_id")
    
    #add missing marital status variable as unknown
    df <- df %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                 TRUE~marital_status))
    
    #1s and 0s for prescription outcomes
    df_outcomes <- df %>%
      select(CDI,overall_tox) %>% mutate(CDI=case_when(CDI==TRUE~1,TRUE~0),
                                         overall_tox=case_when(overall_tox==TRUE~1,TRUE~0))
    
    #select predictors of interest
    df_predictors <- df %>% select(pHADM:age65,prAKI:pDIAB,pCARD:curr_service,pICU:pSEPSIS,ob_freq,highCRP,
                                   temperature:dbp,inpatient_7d,
                                   ab_name_Ampicillin_Ceftriaxone:ab_name_Ampicillin,race:age_grp)
    
    #make dummy variables for predictors
    urdummies <- dummyVars(" ~ .", data = df_predictors)
    df_predictors <- predict(urdummies, newdata = df_predictors)
    df_combined <- as.data.frame(cbind(df_outcomes, df_predictors))
    
    assign(outcnam,df_outcomes,envir = .GlobalEnv)
    assign(prednam,df_predictors,envir = .GlobalEnv)
    assign(combnam,df_combined,envir = .GlobalEnv)
    
  }
  
  ###Compiling hyperparameter dataframe
  hyptabber <- function(map1,map2,params1,params2){
    
    #names of hyperparameters
    hyp_names = c("Learning rate","Maximum tree depth",
                  "Minimum child weight","Row subsample",
                  "Column subsample","N training rounds")
    
    #empty df
    hyptab <- data.frame(matrix(nrow=length(params1[[1]]),ncol=0))
    
    #put hyperparameters in for each urine model
    hyptab <- reduce(seq_along(map1),function(df,i){
      
      df %>% 
        mutate(!!sym((map1[i] %>% unlist())) := params1[[map1[i] ]] %>% unlist())
      
    },.init=hyptab)
    
    #add hyperparameters for prescription model
    hyptab <- hyptab %>% mutate(
      CDI=params2$CDI %>% unlist(),
      Toxicity=params2$overall_tox %>% unlist()) %>% 
      dplyr::slice(-c(1:2))
    
    #set, clean, and arrange columns
    colnames(hyptab) <- abcombo_replace(colnames(hyptab),map2) %>% 
      str_replace_all("_"," & ") %>% str_replace_all("-","/")
    hyptab <- hyptab %>% mutate(Hyperparameter=hyp_names) %>% 
      relocate(Hyperparameter,.before=1) %>% t() %>% data.frame() 
    hyptab$Model <- rownames(hyptab)
    hyptab <- hyptab %>% relocate(Model,.before = 1)
    colnames(hyptab) <- hyptab[1,]
    hyptab %>% rename(Model="Hyperparameter") %>% 
      dplyr::slice(-1) %>% tibble()
    
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
    ur_actclass <- factor(ur_actclass, levels=c(0,1))
    
    #assign to global environment
    assign(probnam,ur_predprobs,envir=.GlobalEnv)
    assign(classnam,ur_predclass,envir=.GlobalEnv)
    assign(actnam,ur_actclass,envir=.GlobalEnv)
    
  }
  
  ###Getting maximum mean and SDs of AUCs in fairness analysis
  fairmetfilter <- function(df,cat){
    
    fairmets %>% filter(Category==cat) %>% group_by(Model,Characteristic) %>%
      summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% 
      summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
                max_sd=max(sd_AUC)) %>% ungroup()
    
  }
  
  ###Putting together fairness metrics
  fair_subfunc <- function(df,char){
    
    top_charac_amd <- max_characdifs %>% arrange(desc(maxAUC_meandif)) %>% dplyr::slice(1)
    top_charac_amdmod <- top_charac_amd %>% select(Model) %>% abcombo_replace(overall_map) %>% tolower()
    top_charac_amdcase <- fairmets %>% filter(Model==names(top_charac_amdmod)&Category==char) %>% arrange(desc(AUC))
    top_charac_amdtop <- top_charac_amdcase %>% dplyr::slice(1) %>% select(AUC) %>% round(3) %>% unlist()
    top_charac_amdtoplab <- top_charac_amdcase %>% dplyr::slice(1) %>% select(Characteristic) %>% unlist()
    top_charac_amdbot <- top_charac_amdcase %>% dplyr::slice(nrow(top_charac_amdcase)) %>% select(AUC) %>% round(3) %>% unlist()
    top_charac_amdbotlab <- top_charac_amdcase %>% dplyr::slice(nrow(top_charac_amdcase)) %>% select(Characteristic) %>% unlist()
    
    top_charac_asd <- max_characdifs %>% arrange(desc(max_sd)) %>% dplyr::slice(1)
    top_charac_asdmod <- top_charac_asd %>% select(Model) %>% abcombo_replace(overall_map) %>% tolower()
    top_charac_asdval <- top_charac_asd %>% select(max_sd) %>% round(3) %>% unlist()
    top_charac_asdlab <- fairmets %>% filter(Category==char) %>% group_by(Model,Characteristic) %>%
      summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% ungroup() %>% arrange(desc(sd_AUC)) %>% 
      select(Characteristic) %>% dplyr::slice(1) %>% unlist()
    
    data.frame(char=char %>% tolower(),top_aucmodel=top_charac_amdmod,top_auclab=top_charac_amdtoplab,top_auc=top_charac_amdtop,
               bot_auclab=top_charac_amdbotlab,bot_auc=top_charac_amdbot,top_sdmodel=top_charac_asdmod,
               top_sdlab=top_charac_asdlab,top_sd=top_charac_asdval) %>% remove_rownames()
    
  }
  
  ###Manuscript stability metrics into table
  stabmetter <- function(df) {
    
    max_stabdifs <- df %>% group_by(Model,Training_size) %>%
      summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% 
      summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
                max_sd=max(sd_AUC)) %>% ungroup()
    
    top_stab_amd <- max_stabdifs %>% arrange(desc(maxAUC_meandif)) %>% dplyr::slice(1)
    top_stab_amdmod <- top_stab_amd %>% select(Model) %>% abcombo_replace(overall_map) %>% unlist()
    top_stab_amdcase <- df %>% filter(Model==names(top_stab_amdmod)) %>% group_by(Training_size) %>% 
      mutate(meanauc=mean(AUC)) %>% ungroup() %>% arrange(desc(meanauc))
    top_stab_amdtop <- top_stab_amdcase %>% dplyr::slice(1) %>% select(meanauc) %>% round(3) %>% unlist()
    top_stab_amdtoplab <- top_stab_amdcase %>% dplyr::slice(1) %>% select(Training_size) %>% unlist()
    top_stab_amdbot <- top_stab_amdcase %>% dplyr::slice(nrow(top_stab_amdcase)) %>% select(meanauc) %>% round(3) %>% unlist()
    top_stab_amdbotlab <- top_stab_amdcase %>% dplyr::slice(nrow(top_stab_amdcase)) %>% select(Training_size) %>% unlist()
    
    top_stab_asd <- max_stabdifs %>% arrange(desc(max_sd)) %>% dplyr::slice(1)
    top_stab_asdmod <- top_stab_asd %>% select(Model) %>% abcombo_replace(overall_map) %>% tolower()
    top_stab_asdval <- top_stab_asd %>% select(max_sd) %>% round(3) %>% unlist()
    top_stab_asdcase <- df %>% filter(Model==names(top_stab_asdmod)) %>% group_by(Training_size) %>% 
      mutate(sdauc=sd(AUC)) %>% ungroup() %>% arrange(desc(sdauc)) %>% dplyr::slice(1) %>% 
      select(Training_size) %>% unlist()
    
    data.frame(top_stab_amdmod,top_stab_amdtoplab,top_stab_amdtop,
               top_stab_amdbotlab,top_stab_amdbot,top_stab_asdmod,
               top_stab_asdcase,top_stab_asdval)
    
  }
  
  ###Manuscript time metrics into table
  timemetter <- function(df) {
    
    max_timemets <- df %>% group_by(Model,Train_year,Test_year) %>%
      summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% group_by(Model) %>%  
      summarise(maxAUC_meandif=max(mean_AUC)-min(mean_AUC),
                max_sd=max(sd_AUC)) %>% ungroup()
    
    top_time_amd <- max_timemets %>% arrange(desc(maxAUC_meandif)) %>% dplyr::slice(1)
    top_time_amdmod <- top_time_amd %>% select(Model) %>% abcombo_replace(overall_map)
    top_time_amdcase <- df %>% filter(Model==names(top_time_amdmod)) %>%
      group_by(Train_year,Test_year) %>% mutate(meanauc=mean(AUC)) %>% ungroup() %>% arrange(desc(meanauc))
    
    top_time_amdtop <- top_time_amdcase %>% dplyr::slice(1) %>% select(meanauc) %>% round(3) %>% unlist()
    top_time_amdtoptrain <- top_time_amdcase %>% dplyr::slice(1) %>% select(Train_year) %>% unlist()
    top_time_amdtoptest <- top_time_amdcase %>% dplyr::slice(1) %>% select(Test_year) %>% unlist()
    top_time_amdbot <- top_time_amdcase %>% dplyr::slice(nrow(top_time_amdcase)) %>% select(meanauc) %>% round(3) %>% unlist()
    top_time_amdbottrain <- top_time_amdcase %>% dplyr::slice(nrow(top_time_amdcase)) %>% select(Train_year) %>% unlist()
    top_time_amdbottest <- top_time_amdcase %>% dplyr::slice(nrow(top_time_amdcase)) %>% select(Test_year) %>% unlist()
    
    top_time_asd <- max_timemets %>% arrange(desc(max_sd)) %>% dplyr::slice(1)
    top_time_asdmod <- top_time_asd %>% select(Model) %>% abcombo_replace(overall_map)
    top_time_asdval <- top_time_asd %>% select(max_sd) %>% round(3) %>% unlist()
    top_time_asdlabs <- df %>% filter(Model==names(top_time_asdmod)) %>% group_by(Train_year,Test_year) %>%
      summarise(mean_AUC=mean(AUC),sd_AUC=sd(AUC)) %>% ungroup() %>% arrange(desc(sd_AUC)) %>% 
      select(Train_year,Test_year) %>% dplyr::slice(1)
    top_time_sdtrainlab <- top_time_asdlabs %>% select(1) %>% unlist()
    top_time_sdtestlab <- top_time_asdlabs %>% select(2) %>% unlist()
    
    data.frame(top_time_amdmod,top_time_amdtoptrain,top_time_amdtoptest,
               top_time_amdtop,top_time_amdbottrain,top_time_amdbottest,
               top_time_amdbot,top_time_asdmod,top_time_sdtrainlab,
               top_time_sdtestlab,top_time_asdval)
    
  }
  
  ###Fairness metrics into table
  fair_subfunc <- function(df1,df2,char){
    
    top_charac_amd <- df2 %>% arrange(desc(maxAUC_meandif)) %>% dplyr::slice(1)
    top_charac_amdmod <- top_charac_amd %>% select(Model) %>% abcombo_replace(overall_map) %>% tolower()
    top_charac_amdcase <- df1 %>% filter(Model==names(top_charac_amdmod)&Category==char) %>% group_by(Characteristic) %>% 
      mutate(meanauc=mean(AUC)) %>% ungroup() %>% arrange(desc(meanauc))
    top_charac_amdtop <- top_charac_amdcase %>% dplyr::slice(1) %>% select(meanauc) %>% round(3) %>% unlist()
    top_charac_amdtoplab <- top_charac_amdcase %>% dplyr::slice(1) %>% select(Characteristic) %>% unlist()
    top_charac_amdbot <- top_charac_amdcase %>% dplyr::slice(nrow(top_charac_amdcase)) %>% select(meanauc) %>% round(3) %>% unlist()
    top_charac_amdbotlab <- top_charac_amdcase %>% dplyr::slice(nrow(top_charac_amdcase)) %>% select(Characteristic) %>% unlist()
    
    top_charac_asd <- df2 %>% arrange(desc(max_sd)) %>% dplyr::slice(1)
    top_charac_asdmod <- top_charac_asd %>% select(Model) %>% abcombo_replace(overall_map) %>% tolower()
    top_charac_asdval <- top_charac_asd %>% select(max_sd) %>% round(3) %>% unlist()
    top_charac_asdlab <- df1 %>% filter(Model==names(top_charac_asdmod)&Category==char) %>% group_by(Characteristic) %>%
      mutate(sd_AUC=sd(AUC)) %>% ungroup() %>% arrange(desc(sd_AUC)) %>% 
      select(Characteristic) %>% dplyr::slice(1) %>% unlist()
    
    data.frame(char=char %>% tolower(),top_aucmodel=top_charac_amdmod,top_auclab=top_charac_amdtoplab,top_auc=top_charac_amdtop,
               bot_auclab=top_charac_amdbotlab,bot_auc=top_charac_amdbot,top_sdmodel=top_charac_asdmod,
               top_sdlab=top_charac_asdlab,top_sd=top_charac_asdval) %>% remove_rownames()
    
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
  metrics_ablist <- c(combined_antimicrobial_map,"CDI","toxicity")
  
  ###Big combined map
  overall_map <- append(combined_antimicrobial_map,abx_outcome_map)

##Performance metrics summary

  ###Metric dataframes
  auc_df <- metric_cleaner("AUC")
  recall_df <- metric_cleaner("Recall")
  precision_df <- metric_cleaner("Precision")
  accuracy_df <- metric_cleaner("Accuracy")
  f1_score_df <- metric_cleaner("F1")
  
  ###Write dataframes to csvs
  write_csv(auc_df,"auc_df.csv")
  write_csv(recall_df,"recall_df.csv")
  write_csv(precision_df,"precision_df.csv")
  write_csv(accuracy_df,"accuracy_df.csv")
  write_csv(f1_score_df,"f1_score_df.csv")
  
  ###Compile and write metric table including confidence intervals
  metric_tablecompiler(auc_df,precision_df,recall_df,f1_score_df,accuracy_df,
                       "ci_singles_table.csv","ci_combos_table.csv",
                       "metrics_singles_table.csv","metrics_combos_table.csv")
  
##Urine model stability analysis

###Final preprocessing and dataset train/test splits - urines dataset
ur_datpreproc(urines5,"urines5_outcomes","urines5_predictors","urines5_combined",
              ur_xg,"ur_xg_outcomes","ur_xg_predictors","ur_xg_combined")

###Read-in chosen model parameters (urine susceptibility prediction)
final_bestparams <- hypparamreader("final_params_",combined_antimicrobial_map)

###Model stability analysis (susceptibility prediction)

###Training dataset size list
tr_size_seq <- c(0.02,0.06,0.10,0.14)

###Empty df for time sens
metrics_biglist <- list()

###Iteration top level: antibiotic models
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    ###Empty top-level list
    metrics_medlist <- list()
    
    ###Iteration level 3: training dataset size
    for(siz in seq_along(tr_size_seq)) {
      
      ###Empty level 3 list
      metrics_litlist <- list()
      
      ###Iteration level 2: seed
      for (seedpick in seq(1,6)) {
        
        ###Empty level 2 list
        metrics_list <- list()
      
        ###Status update
        iterrun <- glue("{outcome} interation {seedpick} for {tr_size_seq[siz]}")
        print(iterrun)
        
        ###Iteration seed
        set.seed(seedpick)
    
        ###Split df to train and test
        urines5_combined %>% TTsplitter(outcome,tr_size_seq[siz],"urines5Train","urines5Test")
        
        if (length(urines5Test[[outcome]] %>% unique()) > 1) {
        
        ###Ensure features line up between dataframes
        ur_xg_combined <- ur_xg_combined %>%
          lineup_features(urines5_predictors,urines5Train)
    
        ###Make xgboost training matrices
        urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
        urtest_matrix <- urines5Test %>% model_matrixmaker(urines5_predictors,outcome)
        urmicro_matrix <- ur_xg_combined %>% model_matrixmaker(urines5_predictors,outcome)
        
        ###Set parameters
        urparams <- list(
          objective = "binary:logistic",
          eval_metric = "auc",
          eta = final_bestparams[[outcome]]$eta,
          max_depth = final_bestparams[[outcome]]$max_depth,
          min_child_weight = final_bestparams[[outcome]]$min_child_weight,
          subsample = final_bestparams[[outcome]]$subsample,
          colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
        )
        
        ###Training
        print(glue("Training for {outcome}"))
        xgb_urinemodel <- xgb.train(
          params = urparams,data = urtrain_matrix,
          nrounds = final_bestparams[[outcome]]$best_nrounds
        )
        
        ###Get predicted probability/class and actual class
        urines5Test %>% probclassactual(xgb_urinemodel,urtest_matrix,outcome,'ur_predprobs',
                                         'ur_predclass','ur_actualclass')
        
        ###Performance metrics
        ur_predact_df <- data.frame(act_val = ur_actualclass, pred_class = ur_predclass,pred_probs = ur_predprobs)
        ur_metrics <- ur_perf_mets(ur_predact_df,NULL,bootstr = F)
        
        ###Add metrics list to level 1 list
        metrics_list[[outcome]] <- ur_metrics
        
        ###Add model metrics to level 2 seed list
        metrics_litlist[[seedpick]] <- metrics_list
        
        } else {
          
          ###If insufficient variance for metrics, add NAs
          metrics_litlist[[seedpick]] <- list(
            AUC = NA,
            Accuracy = NA,
            Precision = NA,
            Recall = NA,
            F1_Score = NA
          )
          
          }
    
        }
    
      ###Add seed list to level 3 training dataset size list
      metrics_medlist[[siz]] <- metrics_litlist
      names(metrics_medlist[[siz]]) <- tr_size_seq[siz]
    
      }
    
    }
  
  ###Add training dataset size list to top level model list
  metrics_biglist[[outcome]] <- metrics_medlist

}

###Convert stability metrics to dataframe and write to csv
metrics_df <- metrics_biglist %>% stabmets_writer()
write_csv(metrics_df,"stability_metrics.csv")

##Model fairness analysis

###Protected characteristics excluding age into list
ur_prot <- urines5_combined %>% protchar_assembler()
protchar_index <- which(grepl("(^marital_status_|^language|^MALE|^race_)",colnames(ur_prot))&
                          !grepl("(UNKNOWN|Other|\\?)",colnames(ur_prot)))
protchars <- colnames(ur_prot)[protchar_index]

###Empty top-level performance metrics list
metrics_biglist2 <- list()

###Top-level iteration: urine prediction models
for (outcome in colnames(urines5_outcomes)[1:13]) {
    
    if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
      
      ###Empty level 3 list
      metrics_medlist2 <- list()
      
      ###Level 3 iteration: seeds
      for (seedpick in seq(1,6)) {
        
        ###Set iteration seed
        set.seed(seedpick)
          
        ###Split df to train and test
        urines5_combined %>% TTsplitter(outcome,0.8,"urines5Train","urines5Test")
        
        ###Ensure protected characteristics mutated
        urines5Test <- urines5Test %>% protchar_assembler()
        
        ###Ensure features line up between dataframes
        ur_xg_combined <- ur_xg_combined %>%
          lineup_features(urines5_predictors,urines5Train)
        
        ###Make xgboost training matrices
        urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
        urtest_matrix <- urines5Test %>% model_matrixmaker(urines5_predictors,outcome)
        urmicro_matrix <- ur_xg_combined %>% model_matrixmaker(urines5_predictors,outcome)
        
        ###Set parameters
        urparams <- list(
          objective = "binary:logistic",
          eval_metric = "auc",
          eta = final_bestparams[[outcome]]$eta,
          max_depth = final_bestparams[[outcome]]$max_depth,
          min_child_weight = final_bestparams[[outcome]]$min_child_weight,
          subsample = final_bestparams[[outcome]]$subsample,
          colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
        )
        
        ###Training
        print(glue("Training for {outcome}"))
        xgb_urinemodel <- xgb.train(
          params = urparams,data = urtrain_matrix,
          nrounds = final_bestparams[[outcome]]$best_nrounds
        )
        
        ###Empty level 2 df
        metrics_litlist2 <- list()
        
        ###Level 2 iteration: protected characteristic
        for(protchar in seq_along(protchar_index)) {
          
          ###Iteration update
          iterrun <- glue("{outcome} iteration {seedpick} for {protchars[protchar]}")
          print(iterrun)
            
          ###Level 1 empty df
          metrics_list2 <- list()
          
          ###Filter testing dataframe to characteristic of interest
          urines5Test2 <- urines5Test[urines5Test[protchar_index[protchar]]==1,]
          
          if (length(urines5Test2[[outcome]] %>% unique()) > 1) {
            
            ###Make xgboost training matrix
            test_matrix2 <- urines5Test2 %>% model_matrixmaker(urines5_predictors,outcome)
            
            ###Get predicted probability/class and actual class
            urines5Test2 %>% probclassactual(xgb_urinemodel,test_matrix2,outcome,'ur_predprobs',
                                             'ur_predclass','ur_actualclass')
            
            ###Performance metrics
            ur_predact_df <- data.frame(act_val = ur_actualclass, pred_class = ur_predclass,pred_probs = ur_predprobs)
            ur_metrics2 <- ur_perf_mets(ur_predact_df,NULL,bootstr = F)
            charlist <- list(protchars[[protchar]])
            names(charlist) <- "Characteristic"
            testsup <- list(nrow(urines5Test2))
            names(testsup) <- "Test_support"
            ur_metrics2 <- c(charlist,ur_metrics2,testsup)
            
            ###Populate level 1 list with metrics
            metrics_list2[[outcome]] <- ur_metrics2
          
            ###Add level 1 metrics list to level 2 characteristic list
            metrics_litlist2[[protchar]] <- metrics_list2
          
            } else {
            
              ###If insufficient variance for metrics, add NAs
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
        
        ###Add level 2 list to level 3 seed list
        metrics_medlist2[[seedpick]] <- metrics_litlist2
        
        }
      
      }
  
  #Add level 3 list to top-level model list
  metrics_biglist2[[outcome]] <- metrics_medlist2
  
  }

###Make fairness dataframe and write to csv
metrics_df2 <- metrics_biglist2 %>% protchar_writer()
write_csv(metrics_df2,"fairness_metrics1.csv")

###Age

###Age list
ages <- urines5_combined %>% distinct(standard_age) %>% unlist() %>% sort()

###Empty stop-level modelperformance metrics list
metrics_biglist3 <- list()

###Top-level iteration: urine models
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    ###Empty level 3 list
    metrics_medlist3 <- list()
    
    ###Level 3 iteration: seeds
    for (seedpick in seq(1,6)) {
      
      ###Iteration seed
      set.seed(seedpick)
        
      ###Split df to train and test
      urines5_combined %>% TTsplitter(outcome,0.8,"urines5Train","urines5Test")
      
      ###Ensure protected characteristics mutated
      urines5Test <- urines5Test %>% protchar_assembler()
      
      ###Ensure features line up between dataframes
      ur_xg_combined <- ur_xg_combined %>%
        lineup_features(urines5_predictors,urines5Train)
      
      ###Make xgboost training matrices
      urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
      urtest_matrix <- urines5Test %>% model_matrixmaker(urines5_predictors,outcome)
      urmicro_matrix <- ur_xg_combined %>% model_matrixmaker(urines5_predictors,outcome)
      
      ###Set parameters
      urparams <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = final_bestparams[[outcome]]$eta,
        max_depth = final_bestparams[[outcome]]$max_depth,
        min_child_weight = final_bestparams[[outcome]]$min_child_weight,
        subsample = final_bestparams[[outcome]]$subsample,
        colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
      )
      
      ###Training
      print(glue("Training for {outcome}"))
      xgb_urinemodel <- xgb.train(
        params = urparams,data = urtrain_matrix,
        nrounds = final_bestparams[[outcome]]$best_nrounds
      )
      
      ###Empty level 2 list
      metrics_litlist3 <- list()
      
      ###Level 2 iteration: ages
      for(age in seq_along(ages)) {
        
        ###Iteration update 
        iterrun <- glue("{outcome} iteration {seedpick} for age {ages[age]}")
        print(iterrun)
        
        ###Empty level 1 list  
        metrics_list3 <- list()
            
        if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
          
          ###Filter test df to age of interest
          urines5Test2 <- urines5Test %>% filter(standard_age==ages[age])
              
          ###Make xgboost training matrix
          test_matrix2 <- urines5Test2 %>% model_matrixmaker(urines5_predictors,outcome)
          
          ###Get predicted probability/class and actual class
          urines5Test2 %>% probclassactual(xgb_urinemodel,test_matrix2,outcome,'ur_predprobs',
                                           'ur_predclass','ur_actualclass')
          
          ###Performance metrics
          ur_predact_df <- data.frame(act_val = ur_actualclass, pred_class = ur_predclass,pred_probs = ur_predprobs)
          ur_metrics3 <- ur_perf_mets(ur_predact_df,NULL,bootstr = F)
          charlist <- list(glue("age{ages[[age]]}"))
          names(charlist) <- "Characteristic"
          testsup <- list(nrow(urines5Test2))
          names(testsup) <- "Test_support"
          ur_metrics3 <- c(charlist,ur_metrics3,testsup)
          
          ###Populate level 1 list with metrics
          metrics_list3[[outcome]] <- ur_metrics3
          
          ###Populate level 2 age list with metrics list
          metrics_litlist3[[age]] <- metrics_list3
          
          } else {
          
            ###If insufficient variance for analysis, add NAs
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
      
      ###Add level 2 age list to level 3 seed list
      metrics_medlist3[[seedpick]] <- metrics_litlist3
      
      }
    
    }

  ###Add level 3 seed list to top-level model list
  metrics_biglist3[[outcome]] <- metrics_medlist3
  
  }

###Write fairness metrics to df, then to CSV
metrics_protchar_df <- metrics_biglist3 %>% 
  agefairmettodf(metrics_df2)
write_csv(metrics_protchar_df,"fairness_metrics.csv")

##Time sensitivity analysis

###Year group indexes in testing dataframe
patkey <- pats %>% select(subject_id,anchor_year_group)
urines_ref <- urines_ref %>% left_join(patkey)

###Set time sequence
time_seq <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")

###Set empty list
metrics_biglist4 <- list()

###Top-level iteration: urine models
for (outcome in colnames(urines5_outcomes)[1:13]) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    ###Empty level 4 list
    metrics_medlist4 <- list()
    
    ###Level 4 iteration: seed
    for (seedpick in seq(1,6)) {
    
      ###Iteration seed
      set.seed(seedpick)
      
      ###Empty level 3 list
      metrics_medmedlist4 <- list()
      
      ###Level 3 iteration: training dataset
      for(tim in seq_along(time_seq)) {
        
        ###Empty level 2 list
        metrics_litlist4 <- list()
        
        ###Filter df to training dataset year and partition
        urines5_filtered1 <- urines5_combined[which(urines_ref$anchor_year_group==time_seq[tim]),]
        trainIndex1 <- createDataPartition(urines5_filtered1[[outcome]], p = 0.8, list = FALSE, times = 1)
        urines5Train <- urines5_filtered1[trainIndex1, ]
        
        ###Level 2 iteration: testing dataset
        for(tim2 in seq_along(time_seq)) {
          
          ###Empty level 1 list
          metrics_list4 <- list()
        
          ###Iteration update
          iterrun <- glue("{outcome} iteration {seedpick} for {time_seq[tim]} training, {time_seq[tim2]} testing")
          print(iterrun)
          
          ###Filter to testing dataset year and index
          urines5_filtered2 <- urines5_combined[which(urines_ref$anchor_year_group==time_seq[tim2]),]
          trainIndex2 <- createDataPartition(urines5_filtered2[[outcome]], p = 0.8, list = FALSE, times = 1)
        
          ###If the same train-test year, ensure 1st holdout set used for testing
          if (time_seq[tim]==time_seq[tim2]) {
          
            urines5Test <- urines5_filtered1[-trainIndex1, ]
          
            print("same time period")
          
            ###If different time period, use second filtered dataset
            } else {
              
              urines5Test <- urines5_filtered2[-trainIndex2, ]
        
              print("different time period")
        
              }
          
          ###Sync predictor columns
          predictor_columns <- colnames(urines5_predictors)
          selected_columns <- intersect(predictor_columns, colnames(urines5Train))
          
          ###XGBoost mitrices
          urtrain_matrix <- urines5Train %>% model_matrixmaker(urines5_predictors,outcome)
          urtest_matrix <- urines5Test %>% model_matrixmaker(urines5_predictors,outcome)
        
          ###Set parameters
          urparams <- list(
            objective = "binary:logistic",
            eval_metric = "auc",
            eta = final_bestparams[[outcome]]$eta,
            max_depth = final_bestparams[[outcome]]$max_depth,
            min_child_weight = final_bestparams[[outcome]]$min_child_weight,
            subsample = final_bestparams[[outcome]]$subsample,
            colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
          )
          
          ###Training
          print(glue("Training for {outcome}"))
          xgb_urinemodel <- xgb.train(
            params = urparams,data = urtrain_matrix,
            nrounds = final_bestparams[[outcome]]$best_nrounds
          )
          
          ###Get predicted probability/class and actual class
          urines5Test %>% probclassactual(xgb_urinemodel,urtest_matrix,outcome,'ur_predprobs',
                                           'ur_predclass','ur_actualclass')
          
          ###Performance metrics
          ur_predact_df <- data.frame(act_val = ur_actualclass, pred_class = ur_predclass,pred_probs = ur_predprobs)
          ur_metrics4<- ur_perf_mets(ur_predact_df,NULL,bootstr = F)
          trainsup <- list(nrow(urines5Train))
          names(trainsup) <- "Train_support"
          testsup <- list(nrow(urines5Test))
          names(testsup) <- "Test_support"
          ur_metrics4 <- c(ur_metrics4,trainsup,testsup)
          
          ###Populate level 1 list with metrics
          metrics_list4[[outcome]] <- ur_metrics4
          
          ###Add level 1 metrics list to level 2 testing data list
          metrics_litlist4[[tim2]] <- metrics_list4
          names(metrics_litlist4[[tim2]]) <- time_seq[tim2]
        
          }
        
        ###Add level 2 testing data list to level 3 training data list
        metrics_medmedlist4[[tim]] <- metrics_litlist4
        names(metrics_medmedlist4[[tim]]) <- time_seq[tim]
      
        }
      
      ###Add level 3 training data list to level 4 seed list
      metrics_medlist4[[seedpick]] <- metrics_medmedlist4
        
      }
      
    }
  
  ###Add level 4 seed list to top-level model list
  metrics_biglist4[[outcome]] <- metrics_medlist4
  
  }

###Time sensitivity metrics to dataframe
metrics_df4 <- metrics_biglist4 %>% timesenstodf()
write_csv(metrics_df4,"time_sens_metrics.csv")

##CDI and toxicity dataset preprocessing

###Joining up abx dataframes
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()
abx <- data.frame(rbind(train_abx,test_abx))

###Preprocessing
abx %>% modeltest_preproc("abx_outcomes","abx_predictors","abx_combined")

###Read in hyperparameters
cdi_tox_final_bestparams <- hypparamreader("cdi_tox_final_params_",colnames(abx_outcomes))

###Model stability analysis (CDI and toxicity prediction)

###Training dataset size list
tr_size_seq <- c(0.02,0.06,0.10,0.14)

###Empty top-level list
metrics_biglist <- list()

###Top-level iteration: prescription outcome models
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    ###Empty level 3 list
    metrics_medlist <- list()
    
    ###Level 3 iteration: training dataset size
    for(siz in seq_along(tr_size_seq)) {
      
      ###Empty level 2 list
      metrics_litlist <- list()
      
      ###Level 2 iteration: seeds
      for (seedpick in seq(1,6)) {
        
        ###Empty level 1 list
        metrics_list <- list()
        
        ###Iteration update
        iterrun <- glue("{outcome} iteration {seedpick} for {tr_size_seq[siz]}")
        print(iterrun)
        
        ###Iteration seed
        set.seed(seedpick)
        
        ###Split df to train and test
        abx_combined %>% TTsplitter(outcome,tr_size_seq[siz],"abxTrain","abxTest")
        
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
        
        ###Training
        print(glue("Training for {outcome}"))
        xgb_abxmodel <- xgb.train(
          params = abxparams,data = abxtrain_matrix,
          nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
        )
        
        ###Get predicted probability/class and actual class
        abxTest %>% probclassactual(xgb_abxmodel,abxtest_matrix,outcome,'abx_predprobs',
                                          'abx_predclass','abx_actualclass')
        
        ###Performance metrics
        abx_predact_df <- data.frame(act_val = abx_actualclass, pred_class = abx_predclass,pred_probs = abx_predprobs)
        abx_metrics <- ur_perf_mets(abx_predact_df,NULL,bootstr = F)
        
        ###Populate level 1 list with metrics
        metrics_list[[outcome]] <- abx_metrics
        
        ###Poulate level 2 seed list with level 1 metric list
        metrics_litlist[[seedpick]] <- metrics_list
        
        }
      
      ###Populate level 3 training size list with level 2 seed list
      metrics_medlist[[siz]] <- metrics_litlist
      names(metrics_medlist[[siz]]) <- tr_size_seq[siz]
      
      }
    
    }
  
  ###Populate top-level prescription model list with level 3 training size list
  metrics_biglist[[outcome]] <- metrics_medlist
  
  }

###Write stability metrics to CSV (CDI tox)
metrics_df <- metrics_biglist %>% stabmets_writer()
write_csv(metrics_df,"cdi_tox_stability_metrics.csv")

##Model fairness analysis (CDI tox)

###Protected characteristics excluding age list
abx_prot <- abx_combined %>% protchar_assembler()
protchar_index <- which(grepl("(^marital_status_|^language|^MALE|^race_)",colnames(abx_prot))&
                          !grepl("(UNKNOWN|Other|\\?)",colnames(abx_prot)))
protchars <- colnames(abx_prot)[protchar_index]

###Empty top-level list
metrics_biglist2 <- list()

###Top-level iteration: prescription outcome model
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    ###Empty level 3 list
    metrics_medlist2 <- list()
    
    ###Level 3 iteration: seed
    for (seedpick in seq(1,6)) {
      
      ###Iteration seed
      set.seed(seedpick)
      
      ###Split df to train and test
      abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
      
      ###Ensure protected characteristics mutated
      abxTest <- abxTest %>% protchar_assembler()
      
      ###Ensure features line up between dataframes
      abx_combined <- abx_combined %>%
        lineup_features(abx_predictors,abxTrain)
      
      ###Make xgboost training matrices
      abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
      abxtest_matrix <- abxTest %>% model_matrixmaker(abx_predictors,outcome)
      abxmicro_matrix <- abx_combined %>% model_matrixmaker(abx_predictors,outcome)
      
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
      
      ###Training
      print(glue("Training for {outcome}"))
      xgb_abxmodel <- xgb.train(
        params = abxparams,data = abxtrain_matrix,
        nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
      )
      
      ###Empty level 2 list
      metrics_litlist2 <- list()
      
      ###Level 2 iteration: protected characteristics
      for(protchar in seq_along(protchar_index)) {
        
        ###Iteration update
        iterrun <- glue("{outcome} iteration {seedpick} for {protchars[protchar]}")
        print(iterrun)
        
        ###Empty level 1 list
        metrics_list2 <- list()
        
        ###Filter testing df by protected characteristic
        abxTest2 <- abxTest[abxTest[protchar_index[protchar]]==1,]
        
        if (length(abxTest2[[outcome]] %>% unique()) > 1) {
          
          ###Make xgboost training matrices
          test_matrix2 <- abxTest2 %>% model_matrixmaker(abx_predictors,outcome)
          
          ###Get predicted probability/class and actual class
          abxTest2 %>% probclassactual(xgb_abxmodel,test_matrix2,outcome,'abx_predprobs',
                                           'abx_predclass','abx_actualclass')
          
          ###Performance metrics
          abx_predact_df <- data.frame(act_val = abx_actualclass, pred_class = abx_predclass,pred_probs = abx_predprobs)
          abx_metrics2 <- ur_perf_mets(abx_predact_df,NULL,bootstr = F)
          charlist <- list(protchars[[protchar]])
          names(charlist) <- "Characteristic"
          testsup <- list(nrow(abxTest2))
          names(testsup) <- "Test_support"
          abx_metrics2 <- c(charlist,abx_metrics2,testsup)
          
          ###Populate level 1 list with metrics
          metrics_list2[[outcome]] <- abx_metrics2
          
          ###Populate level 2 characteristic list with level 1 metric list
          metrics_litlist2[[protchar]] <- metrics_list2
          
          } else {
          
            ###If insufficient variance for metrics, add NAs
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
      
      ###Add level 2 characteristic list to level 3 seed list
      metrics_medlist2[[seedpick]] <- metrics_litlist2
      
      }
    
    }

  ###Add level 3 seed list to top-level prescription model list
  metrics_biglist2[[outcome]] <- metrics_medlist2
  
  }

###Write abx outcome fairness metrics to csv
metrics_df2 <- metrics_biglist2 %>% protchar_writer()
write_csv(metrics_df2,"cdi_tox_fairness_metrics1.csv")

###List of ages
ages <- abx_combined %>% distinct(age_grp) %>% unlist() %>% sort()

###Empty top-level list
metrics_biglist3 <- list()

###Top-level iteration: prescription outcome models
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    ###Empty level 3 list
    metrics_medlist3 <- list()
    
    ###Iteration level 3: seed
    for (seedpick in seq(1,6)) {
      
      ###Iteration seed
      set.seed(seedpick)
      
      ###Split df to train and test
      abx_combined %>% TTsplitter(outcome,0.8,"abxTrain","abxTest")
      
      ###Ensure protected characteristics mutated
      abxTest <- abxTest %>% protchar_assembler()
      
      ###Ensure features line up between dataframes
      abx_combined <- abx_combined %>%
        lineup_features(abx_predictors,abxTrain)
      
      ###Make xgboost training matrices
      abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
      abxtest_matrix <- abxTest %>% model_matrixmaker(abx_predictors,outcome)
      abxmicro_matrix <- abx_combined %>% model_matrixmaker(abx_predictors,outcome)
      
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
      
      ###Training
      print(glue("Training for {outcome}"))
      xgb_abxmodel <- xgb.train(
        params = abxparams,data = abxtrain_matrix,
        nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
      )
      
      ###Empty level 2 list
      metrics_litlist3 <- list()
      
      ###Iteration level 2: age
      for(age in seq_along(ages)) {
        
        ###Iteration update
        iterrun <- glue("{outcome} iteration {seedpick} for age {ages[age]}")
        print(iterrun)
        
        ###Empty level 1 list
        metrics_list3 <- list()
        
        if (sum(!is.na(abx_combined[[outcome]])) > 0) {
          
          ###Filter testing df to age of interest
          abxTest2 <- abxTest %>% filter(age_grp==ages[age])
          
          ###Make xgboost test matrix
          test_matrix2 <- abxTest2 %>% model_matrixmaker(abx_predictors,outcome)
          
          ###Get predicted probability/class and actual class
          abxTest2 %>% probclassactual(xgb_abxmodel,test_matrix2,outcome,'abx_predprobs',
                                            'abx_predclass','abx_actualclass')
          
          ###Performance metrics
          abx_predact_df <- data.frame(act_val = abx_actualclass, pred_class = abx_predclass,pred_probs = abx_predprobs)
          abx_metrics3 <- ur_perf_mets(abx_predact_df,NULL,bootstr = F)
          charlist <- list(glue("age{ages[[age]]}"))
          names(charlist) <- "Characteristic"
          testsup <- list(nrow(abxTest2))
          names(testsup) <- "Test_support"
          abx_metrics3 <- c(charlist,abx_metrics3,testsup)
          
          ###Add metrics to level 1 list
          metrics_list3[[outcome]] <- abx_metrics3
          
          ###Add level 1 metrics list to level 2 age list
          metrics_litlist3[[age]] <- metrics_list3
          
          } else {
          
            ###If insufficient variance for metrics, add NAs
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
      
      ###Add level 2 age list to level 3 seed list
      metrics_medlist3[[seedpick]] <- metrics_litlist3
      
      }
    
    }
  
  ###Add level 3 seed list to top-level prescription model list
  metrics_biglist3[[outcome]] <- metrics_medlist3

}

###Write fairness metrics to df, then to CSV
metrics_protchar_df <- metrics_biglist3 %>% 
  agefairmettodf(metrics_df2)
write_csv(metrics_protchar_df,"cdi_tox_fairness_metrics.csv")

##List of time periods
patkey <- pats %>% select(subject_id,anchor_year_group)
urines_ref <- urines_ref %>% left_join(patkey)
time_seq <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")

###Empty top-level list
metrics_biglist4 <- list()

###Top-level iteration: prescription model outcome
for (outcome in colnames(abx_outcomes)) {
  
  if (sum(!is.na(abx_combined[[outcome]])) > 0) {
    
    ###Empty level 4 list
    metrics_medlist4 <- list()
    
    ###Level 4 iteration: seed
    for (seedpick in seq(1,6)) {
      
      ###Iteration seed
      set.seed(seedpick)
      
      ###Empty level 3 list
      metrics_medmedlist4 <- list()
      
      ###Level 3 iteration: training dataset time period
      for(tim in seq_along(time_seq)) {
        
        ###Empty level 2 list
        metrics_litlist4 <- list()
        
        ###Filter to training dataset time period of interest
        abx_filtered1 <- abx_combined[which(urines_ref$anchor_year_group==time_seq[tim]),]
        trainIndex1 <- createDataPartition(abx_filtered1[[outcome]], p = 0.8, list = FALSE, times = 1)
        abxTrain <- abx_filtered1[trainIndex1, ]
        
        ###Iteration level 2: testing dataset time period
        for(tim2 in seq_along(time_seq)) {
          
          ###Empty level 1 list
          metrics_list4 <- list()
          
          ###Iteration update
          iterrun <- glue("{outcome} iteration {seedpick} for {time_seq[tim]} training, {time_seq[tim2]} testing")
          print(iterrun)
          
          ###Filter to testing dataset time period of interest
          abx_filtered2 <- abx_combined[which(urines_ref$anchor_year_group==time_seq[tim2]),]
          trainIndex2 <- createDataPartition(abx_filtered2[[outcome]], p = 0.8, list = FALSE, times = 1)
          
          ###If train and test are from same time period, use holdout set from 1
          if (time_seq[tim]==time_seq[tim2]) {
            
            abxTest <- abx_filtered1[-trainIndex1, ]
            
            ###If train and test are from different periods, use dataset 2
            } else {
            
              abxTest <- abx_filtered2[-trainIndex2, ]
            
              }
          
          ###Sync up predictor solumns
          predictor_columns <- colnames(abx_predictors)
          selected_columns <- intersect(predictor_columns, colnames(abxTrain))
          
          ###XGBoost matrices
          abxtrain_matrix <- abxTrain %>% model_matrixmaker(abx_predictors,outcome)
          abxtest_matrix <- abxTest %>% model_matrixmaker(abx_predictors,outcome)
          
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
          
          ###Training
          print(glue("Training for {outcome}"))
          xgb_abxmodel <- xgb.train(
            params = abxparams,data = abxtrain_matrix,
            nrounds = cdi_tox_final_bestparams[[outcome]]$best_nrounds
          )
          
          ###Get predicted probability/class and actual class
          abxTest %>% probclassactual(xgb_abxmodel,abxtest_matrix,outcome,'abx_predprobs',
                                       'abx_predclass','abx_actualclass')
          
          ###Performance metrics
          abx_predact_df <- data.frame(act_val = abx_actualclass, pred_class = abx_predclass,pred_probs = abx_predprobs)
          abx_metrics4<- ur_perf_mets(abx_predact_df,NULL,bootstr = F)
          trainsup <- list(nrow(abxTrain))
          names(trainsup) <- "Train_support"
          testsup <- list(nrow(abxTest))
          names(testsup) <- "Test_support"
          abx_metrics4 <- c(abx_metrics4,trainsup,testsup)
          
          ###Populate level 1 list with metrics
          metrics_list4[[outcome]] <- abx_metrics4
          
          ###Add level 1 metrics list to level 2 testing dataset list
          metrics_litlist4[[tim2]] <- metrics_list4
          names(metrics_litlist4[[tim2]]) <- time_seq[tim2]
          
          }
        
        ###Add level 2 testing data list to level 3 training data list
        metrics_medmedlist4[[tim]] <- metrics_litlist4
        names(metrics_medmedlist4[[tim]]) <- time_seq[tim]
        
        }
      
      ###Add level 3 training data list to level 4 seed list
      metrics_medlist4[[seedpick]] <- metrics_medmedlist4
      
      }
    
    }
  
  ###Add level 4 seed list to top-level prescription outcome model list
  metrics_biglist4[[outcome]] <- metrics_medlist4
  
  }


###Write time sens metrics to df then csv
metrics_df4 <- metrics_biglist4 %>% timesenstodf()
write_csv(metrics_df4,"cdi_tox_time_sens_metrics.csv")

##Join urine and prescription model dataframes together

###Stability metrics
metrics_dfA <- read_csv("stability_metrics.csv")
metrics_df<- data.frame(rbind(metrics_dfA,metrics_df))
metrics_df <- metrics_df %>% rename(Model="Antimicrobial")
write_csv(metrics_df,"overall_stability_metrics.csv")

###Fairness metrics
metrics_protchar_dfA <- read_csv("fairness_metrics.csv")
metrics_protchar_df<- data.frame(rbind(metrics_protchar_dfA,metrics_protchar_df))
metrics_protchar_df <- metrics_protchar_df %>% rename(Model="Antimicrobial")
write_csv(metrics_protchar_df,"overall_fairness_metrics.csv")

###Time metrics
metrics_df4a <- read_csv("time_sens_metrics.csv")
metrics_df4 <- read_csv("cdi_tox_time_sens_metrics.csv")
metrics_df4 <- data.frame(rbind(metrics_df4a,metrics_df4))
metrics_df4 <- metrics_df4 %>% rename(Model="Antimicrobial")
write_csv(metrics_df4,"overall_time_sens_metrics.csv")

##Model testing plots

###Stability analysis
allmetrics <- c("AUC-ROC","Accuracy","Precision","Recall","F1 score")
allmet_vars <- c("AUC","Accuracy","Precision","Recall","F1_score")
for (i in seq_along(allmetrics)) {
  
  stability_plot(metrics_df,!!sym(allmet_vars[i]),allmetrics[i])
  
  }



###Fairness analysis
for (i in seq_along(allmetrics)){
  
  protchar_plot(metrics_protchar_df,"Age group",!!sym(allmet_vars[i]),
                allmetrics[i],"age groups")
  protchar_plot(metrics_protchar_df,"(Gender|Language|Marital)",!!sym(allmet_vars[i]),
                allmetrics[i],"genders, languages, and marital statuses")
  protchar_plot(metrics_protchar_df,"Race",!!sym(allmet_vars[i]),
                allmetrics[i],"racial groups")
  
  }

###Time sensitivity analysis
yearlist <- c("2008 - 2010","2011 - 2013","2014 - 2016","2017 - 2019")
for (i in seq_along(yearlist)) {

  for (j in seq_along(allmetrics)){
    
    timesens_plot(metrics_df4,yearlist[i],!!sym(allmet_vars[i]),allmetrics[i])
    
    }

  }

##Hyperparameter dataframe
hyp_tab <- fullmap %>% hyptabber(combined_antimicrobial_map,
                                 final_bestparams,
                                 cdi_tox_final_bestparams)
write_csv(hyp_tab,"hyp_tab.csv")
hyp_abs <- c(names(antimicrobial_map),"CDI","Toxicity")
hyp_tall_singles <- hyp_tab %>% filter(Model %in% hyp_abs)
write_csv(hyp_tall_singles,"hyp_tall_singles.csv")

##Values for results

###Stability analysis
stabmets <- read_csv("overall_stability_metrics.csv")
stabmet_df <- stabmetter(stabmets)
write_csv(stabmet_df,"stabmet_df.csv")

###Year group cluster analysis
timemets <- read_csv("overall_time_sens_metrics.csv")
timemet_df <- timemetter(timemets)
write_csv(timemet_df,"timemet_df.csv")

###Fairness analysis
fairmets <- read_csv("overall_fairness_metrics.csv")
max_racedifs <- fairmets %>% fairmetfilter("Race")
max_agedifs <- fairmets %>% fairmetfilter("Age group")
max_sexdifs <- fairmets %>% fairmetfilter("Gender")
fairness_writedf <- rbind(
fairmets %>% fair_subfunc(max_racedifs,"Race"),
fairmets %>% fair_subfunc(max_agedifs,"Age group"),
fairmets %>% fair_subfunc(max_sexdifs,"Gender"))
write_csv(fairness_writedf,"fairness_writedf.csv")


