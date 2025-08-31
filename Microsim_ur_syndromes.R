#SENSITIVITY ANALYSIS (URINARY SYMPTOMS OR DIAGNOSES ONLY)

set.seed(123)

##Functions

###Utility calculation
U_calculation <- function(df,b=1) {
  
  df %>%
    
    #calculate U
    mutate(U=S*(util_tox+util_CDI+util_uti+util_access+util_oral+
                  util_reserve+util_highcost+(util_iv*exp(acuity*b)))) %>% 
    
    mutate(
      
      #intravenous U
      Urosepsis_U = case_when(
        util_iv ==0 ~0,
        TRUE~U),
      
      #oral U
      Outpatient_U = case_when(
        util_oral ==0 ~0,
        TRUE~U))
  
}

###Factorise training and testing datasets
factorise <- function(df) {
  
  #factors for prescription outcomes
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Utility data visualisation
utility_plot <- function(df, variable,application,modification="") {
  
  #quosure
  variable <- enquo(variable)
  
  #colour coding of formulary agents
  axiscols <- if_else(
    df %>% pull(Antimicrobial) %in% formulary_agents,
    "seagreen", "black")
  
  #clean antibiotic combination names
  df <- df %>% mutate(Antimicrobial = str_replace_all(Antimicrobial,"_"," & "))
  
  #filter for ivs
  if (application=="Intravenous treatment") {
    
    df <- df %>% filter(Urosepsis_U != min(Urosepsis_U))
    
    #filter for orals
  } else if (application=="Oral treatment") {
    
    df <- df %>% filter(Outpatient_U != min(Outpatient_U))
    
    #filter for ast (not used in this study)
  } else if (application=="AST") {
    
    df <- df %>% filter(AST_utility != min(AST_utility) & single_agent)
    
  }
  
  #filter for single agents
  if (grepl("single",modification,ignore.case=T)) {
    
    df <- df %>% filter(single_agent)
    
    #filter for combinations
  } else if (grepl("combination",modification,ignore.case=T)) {
    
    df <- df %>% filter(!single_agent)
    
  }
  
  #plot
  thisplot <- ggplot(df %>% mutate(Antimicrobial = 
                                     factor(Antimicrobial,
                                            levels = df %>% group_by(Antimicrobial) %>% 
                                              summarise(Median_util=median(!!variable)) %>% 
                                              arrange(Median_util) %>% select(Antimicrobial) %>% unlist())), 
                     aes(x=!!variable,y=Antimicrobial,fill=Antimicrobial)) +
    
    #box plot
    geom_boxplot(outlier.color = NULL,
                 outlier.alpha = 0.3) +
    
    #theme and titles
    theme_minimal() +
    ggtitle(glue("{application} value distribution on\nsimulation dataset {modification}\n(patients with coded UTI diagnoses)")) +
    theme(legend.position = "None",axis.text.y = element_text(
      colour = axiscols))+
    xlab(glue("{application} treatment value"))
  
  #save to pdf
  ggsave(glue("ur_synd_utility_{application}_{modification}.pdf"), plot = thisplot, device = "pdf", width = 10, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(thisplot)
  
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
  
  fullmap <- abmap %>% unlist()
  names(fullmap) <- NULL
  abmap2 <- abmap[names(abmap) %in% abreflist]
  shortmap <- abmap2 %>% unlist()
  names(shortmap) <- NULL
  
  assign(fullname,fullmap,envir = .GlobalEnv)
  assign(mainname,abmap2,envir = .GlobalEnv)
  assign(shortname,shortmap,envir = .GlobalEnv)
  
}

###Barplot of pathogen coverage
overall_bar_plot <- function(dfs1,dfs2,filter_nonterms,pos_term,invert_term,n_cats,chart_cat) {
  
  #sub-function to get opposites of categories
  inverter <- function(df,n_cat,inverse,type) {
    
    #get percentage of inverse case
    df2 <- df %>% mutate(Perc2=Percentage-lead(Percentage,n=4)) %>% 
      
      #slice to number of classes
      select(Perc2) %>% dplyr::slice(1:n_cat) %>% 
      
      #clean df and add illness severity weights
      rename(Percentage="Perc2") %>% mutate(Metric=inverse,
                                            Weight =seq(0,n_cat-1),
                                            Type=type)
    
    #bind dfs together
    rbind(df,df2) %>% tibble() %>% filter(Metric!="All agents")
  }
  
  #filter to cases of interest
  dfs1a <- dfs1 %>% filter(!grepl(filter_nonterms,Metric,ignore.case=T))
  dfs2a <- dfs2 %>% filter(!grepl(filter_nonterms,Metric,ignore.case=T))
  
  #bind dfs together with inverse cases added
  combined_df <- rbind(dfs1a %>% inverter(n_cats,invert_term,"AUF"),
                       dfs2a %>% inverter(n_cats,invert_term,"Human")) %>% tibble()
  
  #factorise metric
  combined_df$Metric=factor(combined_df$Metric,levels=c(invert_term,pos_term))
  
  #plot
  barplot <- ggplot(combined_df, aes(x = Weight, y = Percentage, group=Metric,fill = Metric)) +
    
    #bars for auf
    geom_bar(stat = "identity", position = position_stack(), data = subset(combined_df, Type == "AUF"),aes(x = Weight - 0.125),width=0.25) +
    
    #bars for human prescriptions
    geom_bar(stat = "identity", position = position_stack(), data = subset(combined_df, Type == "Human"), aes(x = Weight + 0.125),width=0.25)+
    
    #theme and titles
    theme_minimal()+
    ggtitle(glue("{chart_cat}\nof empirical treatment choices and coverage of urine culture isolates, by illness severity"))+
    xlab("Illness severity score")+
    ylab("Isolates susceptible to treatment (%)")+
    
    #space under axis for labels
    ylim(-5,100)+
    
    #no grids
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )+
    
    #add auf and human bar labels
    annotate("text", x = - 0.125, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 0.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x =  0.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 1.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 1.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 2.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x =  2.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 3.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") 
  
  #save to pdf
  ggsave(glue("bar_{chart_cat}.pdf"), plot = barplot, device = "pdf", width = 10, height = 3,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
}

###Barplot of pathogen resistance
overallr_bar_plot <- function(dfs1,dfs2,filter_nonterms,pos_term,invert_term,n_cats,chart_cat) {
  
  #sub-function to get opposites of categories
  inverter <- function(df,n_cat,inverse,type) {
    
    #get percentage of inverse case
    df2 <- df %>% mutate(Perc2=Percentage-lead(Percentage,n=4)) %>% 
      
      #slice to number of classes
      select(Perc2) %>% dplyr::slice(1:n_cat) %>% 
      
      #clean df and add illness severity weights
      rename(Percentage="Perc2") %>% mutate(Metric=inverse,
                                            Weight =seq(0,n_cat-1),
                                            Type=type)
    
    #bind dfs together
    rbind(df,df2) %>% tibble() %>% filter(Metric!="All agents")
  }
  
  #filter to cases of interest
  dfs1a <- dfs1 %>% filter(!grepl(filter_nonterms,Metric,ignore.case=T))
  dfs2a <- dfs2 %>% filter(!grepl(filter_nonterms,Metric,ignore.case=T))
  
  #bind dfs together with inverse cases added
  combined_df <- rbind(dfs1a %>% inverter(n_cats,invert_term,"AUF"),
                       dfs2a %>% inverter(n_cats,invert_term,"Human")) %>% tibble()
  
  #factorise metric
  combined_df$Metric=factor(combined_df$Metric,levels=c(invert_term,pos_term))
  
  #plot
  barplot <- ggplot(combined_df, aes(x = Weight, y = Percentage, group=Metric,fill = Metric)) +
    
    #bars for auf
    geom_bar(stat = "identity", position = position_stack(), data = subset(combined_df, Type == "AUF"),aes(x = Weight - 0.125),width=0.25) +
    
    #bars for human prescriptions
    geom_bar(stat = "identity", position = position_stack(), data = subset(combined_df, Type == "Human"), aes(x = Weight + 0.125),width=0.25)+
    
    #theme and titles
    theme_minimal()+
    ggtitle(glue("{chart_cat}\nof empirical treatment choices and resistance of urine culture isolates, by illness severity"))+
    xlab("Illness severity score")+
    ylab("Isolates resistant to treatment (%)")+
    
    #space under axis for labels
    ylim(-5,100)+
    
    #no grids
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )+
    
    #add auf and human bar labels
    annotate("text", x = - 0.125, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 0.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x =  0.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 1.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 1.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 2.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x =  2.875, y = -5, label = "ADA", size = 4, fontface = "bold",colour="grey40") +
    annotate("text", x = 3.125, y = -5, label = "Human", size = 4, fontface = "bold",colour="grey40") 
  
  #save to pdf
  ggsave(glue("bar_r_{chart_cat}.pdf"), plot = barplot, device = "pdf", width = 10, height = 3,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
}

###Prioritisation by treatment utility
AUF_1 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  
  #filter to specimen of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #filter to antimicrobials in the specified list
    filter(Antimicrobial%in%ab_list) %>%
    
    #arrange in descending order of U
    arrange(desc(U)) %>% select(Antimicrobial,U) %>% 
    
    #select ast panel
    mutate(U = round(U,1)) %>% dplyr::slice(1:panel_size) %>% 
    
    #rename columns
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "U")
  
}

###Prioritisation by Intravenous utility
AUF_IV = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  
  #filter to specimen of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #filter to antimicrobials specified in list
    filter(Antimicrobial%in%ab_list) %>% 
    
    #arrange in descending order of iv U
    arrange(desc(Urosepsis_U)) %>% select(Antimicrobial,Urosepsis_U) %>% 
    
    #select ast panel
    mutate(Urosepsis_U = round(Urosepsis_U,1)) %>% dplyr::slice(1:panel_size) %>% 
    
    #rename columns for table
    rename(`Antimicrobial ranking` = "Antimicrobial",`Intravenous Rx Utility` = "Urosepsis_U")
  
}

###Prioritisation by Oral utility
AUF_PO = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  
  #filter to specimen of interest
  df %>% filter(micro_specimen_id==spec_id) %>%
    
    #filter to antimicrobials specified in list
    filter(Antimicrobial%in%ab_list) %>%
    
    #arrange in descending order of oral utility
    arrange(desc(Outpatient_U)) %>% select(Antimicrobial,Outpatient_U) %>% 
    
    #select ast panel
    mutate(Outpatient_U = round(Outpatient_U,1)) %>% dplyr::slice(1:panel_size) %>% 
    
    #amend names for table
    rename(`Antimicrobial ranking` = "Antimicrobial",`Oral Rx Utility` = "Outpatient_U")
  
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
    sapply(function(x) if (x %in% names(map)) map[names(map) == x] %>% unlist() else x)
  
}

###Assigning treatment recommendations
assign_PDRx <- function(df,probab_df,method_used,ab_list1=ab_singles) {
  
  #empty test recommendation df
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  #iterate over specimens
  for (i in 1:nrow(df)) {
    
    #make recommendations based on U descending order
    rec <- probab_df %>% AUF_1(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                               ab_list=ab_list1) %>% 
      select(1)
    
    #bind to recommendation df
    test_recs <- cbind(test_recs,rec)
    
    #status update
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  #transpose recs df and bind specimen ids
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  #make vector of column names for recommendation columns
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  colnames(test_recs) <- testrec_cols
  
  #join recommendations to probability df
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Intravenous treatment recommendations
assign_Intravenous <- function(df,probab_df,method_used,ab_list1=iv_ab_singles) {
  
  #empty test recommendations df
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  #iterate over specimens
  for (i in 1:nrow(df)) {
    
    #specimen recommendations based on descending u
    rec <- probab_df %>% AUF_IV(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    #add to recommendations df
    test_recs <- cbind(test_recs,rec)
    
    #status update
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  #transpose df and bind specimen numbers
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  #make and add column names for recommendations
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  colnames(test_recs) <- testrec_cols
  
  #add recommendation columns to probability df
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Oral treatment recommendations
assign_Oral <- function(df,probab_df,method_used,ab_list1=oral_ab_singles) {
  
  #empty recommendations df
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  #iterate over specimens
  for (i in 1:nrow(df)) {
    
    #make recommendations based on descending order of U
    rec <- probab_df %>% AUF_PO(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    #bind to recommendations df
    test_recs <- cbind(test_recs,rec)
    
    #status update
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  #transpose df and add specimen ids
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  #iterate over antibiotic list
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning standard recommendations
assign_standard_AST <- function(df,choice_1=NA,choice_2=NA,choice_3=NA,
                                choice_4=NA,choice_5=NA,choice_6=NA) {
  
  #manually compose standard panel
  df %>% mutate(
    STANDARD_AST_1 = choice_1,
    STANDARD_AST_2 = choice_2,
    STANDARD_AST_3 = choice_3,
    STANDARD_AST_4 = choice_4,
    STANDARD_AST_5 = choice_5,
    STANDARD_AST_6 = choice_6
  )
  
}
assign_standard_IV <- function(df,choice_1=NA,choice_2=NA,choice_3=NA,
                               choice_4=NA,choice_5=NA,choice_6=NA) {
  
  #manually compose standard iv panel
  df %>% mutate(
    STANDARD_IV_1 = choice_1,
    STANDARD_IV_2 = choice_2,
    STANDARD_IV_3 = choice_3,
    STANDARD_IV_4 = choice_4,
    STANDARD_IV_5 = choice_5,
    STANDARD_IV_6 = choice_6
  )
  
}
assign_standard_oral <- function(df,choice_1=NA,choice_2=NA,choice_3=NA,
                                 choice_4=NA,choice_5=NA,choice_6=NA) {
  
  #manually compose standard oral panel
  df %>% mutate(
    STANDARD_PO_1 = choice_1,
    STANDARD_PO_2 = choice_2,
    STANDARD_PO_3 = choice_3,
    STANDARD_PO_4 = choice_4,
    STANDARD_PO_5 = choice_5,
    STANDARD_PO_6 = choice_6
  )
  
}

###Binding labels for plot dataframes
label_binder <- function(vec,label) {
  
  #bind labels to plot dataframe
  data.frame(vec,Metric=label,Weight=weightseq)
  
}

###Illness severity result checker
result_checker <- function(df,probs_df,acuity_score) {
  
  #filter to acuity score of interest
  df <- df %>% filter(acuity==acuity_score)  
  
  #filter to cases in microsimulation df
  probs_df <- probs_df %>% semi_join(
    df,by="micro_specimen_id"
  )
  
  #make starting df copy for later
  urines_abx <- df
  
  #get number of s/i results for iv, oral, and overall AUF recommendations
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='S'|PDIVRx_1_result=='I'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='S'|PDPORx_1_result=='I'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='S'|PDRx_1_result=='I'))
  iv_resr <- nrow(df %>% filter(PDIVRx_1_result=='R'))
  po_resr <- nrow(df %>% filter(PDPORx_1_result=='R'))
  overall_resr <- nrow(df %>% filter(PDRx_1_result=='R'))
  
  #get percentage susceptibility rates
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  ivr_perc <- iv_resr/nrow(df)*100
  por_perc <- po_resr/nrow(df)*100
  overallr_perc <- overall_resr/nrow(df)*100
  
  #access iv s percentage
  iv_s_access <- (nrow(df %>% filter((PDIVRx_1_result=='S'|PDIVRx_1_result=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  
  #access oral s percentage
  po_s_access <- (nrow(df %>% filter((PDPORx_1_result=='S'|PDPORx_1_result=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  
  #access s percentage overalls
  overall_s_access <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df))*100
  
  #oral s percentage overall
  overall_s_oral <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
  
  #iv s percentage overalls
  overall_s_iv <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df))*100
  
  #access r percentage overall
  overall_r_access <- (nrow(df %>% filter((PDRx_1_result=='R') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df))*100
  
  #oral r percentage overall
  overall_r_oral <- (nrow(df %>% filter((PDRx_1_result=='R') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
  
  #iv r percentage overall
  overall_r_iv <- (nrow(df %>% filter((PDRx_1_result=='R') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df))*100
  
  #key with auf recommendation 1
  urkey <- urines_abx %>% select(micro_specimen_id,PDRx_1) %>% 
    mutate(Antimicrobial=ab_name(PDRx_1) %>% str_replace("/","-")) %>% 
    select(-PDRx_1)
  
  #count of auf recommendations
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  
  #percentages for auf recommendations
  abrx1_percs <- abrx1_df %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  #assign to global environment
  iv_perc <<- iv_perc
  po_perc <<- po_perc
  overall_perc <<- overall_perc
  ivr_perc <<- ivr_perc
  por_perc <<- por_perc
  overallr_perc <<- overallr_perc
  iv_s_access <<- iv_s_access
  po_s_access <<- po_s_access
  overall_s_access <<- overall_s_access
  overall_s_oral <<- overall_s_oral
  overall_s_iv <<- overall_s_iv
  overall_r_access <<- overall_r_access
  overall_r_oral <<- overall_r_oral
  overall_r_iv <<- overall_r_iv
  abrx1_df <<- abrx1_df
  abrx1_percs <<- abrx1_percs
  
  ###Antibiotics prescribed
  overall_res_px_abx <- nrow(df %>% filter(Px_Abx_result=='S'|Px_Abx_result=='I'))
  overall_perc_px_abx <- overall_res_px_abx/nrow(df)*100
  overall_resr_px_abx <- nrow(df %>% filter(Px_Abx_result=='R'))
  overallr_perc_px_abx <- overall_resr_px_abx/nrow(df)*100
  
  overall_s_access_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                                   Px_Abx %in% access_modified_abx_map))/
                                nrow(df))*100
  overall_s_oral_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                                 Px_Abx %in% oral_modified_abx_map))/
                              nrow(df))*100
  overall_s_iv_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                               Px_Abx %in% iv_modified_abx_map))/
                            nrow(df))*100
  
  overall_r_access_px_abx <- (nrow(df %>% filter((Px_Abx_result=='R') &
                                                   Px_Abx %in% access_modified_abx_map))/
                                nrow(df))*100
  overall_r_oral_px_abx <- (nrow(df %>% filter((Px_Abx_result=='R') &
                                                 Px_Abx %in% oral_modified_abx_map))/
                              nrow(df))*100
  overall_r_iv_px_abx <- (nrow(df %>% filter((Px_Abx_result=='R') &
                                               Px_Abx %in% iv_modified_abx_map))/
                            nrow(df))*100
  
  urkey_abx_px <- urines_abx %>% select(micro_specimen_id,Prescribed_abx) %>% 
    mutate(Antimicrobial=Prescribed_abx) %>% 
    select(-Prescribed_abx)
  
  abrx1_df_px_abx <- df %>% count(Px_Abx) %>% arrange(desc(n))
  abrx1_percs_px_abx <- abrx1_df_px_abx %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  overall_perc_px_abx <<- overall_perc_px_abx
  overallr_perc_px_abx <<- overallr_perc_px_abx
  overall_s_access_px_abx <<- overall_s_access_px_abx
  overall_s_oral_px_abx <<- overall_s_oral_px_abx
  overall_s_iv_px_abx <<- overall_s_iv_px_abx
  overall_r_access_px_abx <<- overall_r_access_px_abx
  overall_r_oral_px_abx <<- overall_r_oral_px_abx
  overall_r_iv_px_abx <<- overall_r_iv_px_abx
  abrx1_df_px_abx <<- abrx1_df_px_abx
  abrx1_percs_px_abx <<- abrx1_percs_px_abx
  
}

###Counting the number of S or I results per antimicrobial per panel
number_results_per_panel <- function(df,start_col,end_col,which_abs,result_1,result_2) {
  
  start_col <- enquo(start_col)
  end_col <- enquo(end_col)
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == result_1 |
                    (df %>%
                       select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                       dplyr::slice(i)) == result_2)
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Results per panel across all variants of UF and standard approaches
rpp_ast <- function(df) {
  
  df$PDRx_rpp_ass <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_singles,"S","I")
  print(1/60)
  df$PDRx_rpp_acs <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_combos,"S","I")
  print(2/60)
  df$PDRx_rpp_aas <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_abs,"S","I")
  print(3/60)
  df$PDRx_rpp_iss <- df %>% number_results_per_panel(PDRx_1,PDRx_6,iv_singles,"S","I")
  print(4/60)
  df$PDRx_rpp_ics <- df %>% number_results_per_panel(PDRx_1,PDRx_6,iv_combos,"S","I")
  print(5/60)
  df$PDRx_rpp_ais <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_ivs,"S","I")
  print(6/60)
  df$PDRx_rpp_oss <- df %>% number_results_per_panel(PDRx_1,PDRx_6,oral_singles,"S","I")
  print(7/60)
  df$PDRx_rpp_ocs <- df %>% number_results_per_panel(PDRx_1,PDRx_6,oral_combos,"S","I")
  print(8/60)
  df$PDRx_rpp_aos <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_orals,"S","I")
  print(9/60)
  df$PDRx_rpp_acss <- df %>% number_results_per_panel(PDRx_1,PDRx_6,access_singles,"S","I")
  print(10/60)
  df$PDRx_rpp_accs <- df %>% number_results_per_panel(PDRx_1,PDRx_6,access_combos,"S","I")
  print(11/60)
  df$PDRx_rpp_aacs <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_access,"S","I")
  print(12/60)
  df$PDRx_rpp_wass <- df %>% number_results_per_panel(PDRx_1,PDRx_6,watch_singles,"S","I")
  print(13/60)
  df$PDRx_rpp_wacs <- df %>% number_results_per_panel(PDRx_1,PDRx_6,watch_combos,"S","I")
  print(14/60)
  df$PDRx_rpp_awas <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_watch,"S","I")
  print(15/60)
  df$PDRx_rpp_asr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_singles,"R","NT")
  print(16/60)
  df$PDRx_rpp_acr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_combos,"R","NT")
  print(17/60)
  df$PDRx_rpp_aar <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_abs,"R","NT")
  print(18/60)
  df$PDRx_rpp_isr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,iv_singles,"R","NT")
  print(19/60)
  df$PDRx_rpp_icr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,iv_combos,"R","NT")
  print(20/60)
  df$PDRx_rpp_air <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_ivs,"R","NT")
  print(21/60)
  df$PDRx_rpp_osr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,oral_singles,"R","NT")
  print(22/60)
  df$PDRx_rpp_ocr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,oral_combos,"R","NT")
  print(23/60)
  df$PDRx_rpp_aor <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_orals,"R","NT")
  print(24/60)
  df$PDRx_rpp_acsr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,access_singles,"R","NT")
  print(25/60)
  df$PDRx_rpp_accr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,access_combos,"R","NT")
  print(26/60)
  df$PDRx_rpp_aacr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_access,"R","NT")
  print(27/60)
  df$PDRx_rpp_wasr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,watch_singles,"R","NT")
  print(28/60)
  df$PDRx_rpp_wacr <- df %>% number_results_per_panel(PDRx_1,PDRx_6,watch_combos,"R","NT")
  print(29/60)
  df$PDRx_rpp_awar <- df %>% number_results_per_panel(PDRx_1,PDRx_6,all_watch,"R","NT")
  print(30/60)
  
  ###Standard results per panel
  df$STANDARD_AST_rpp_ass <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_singles,"S","I")
  print(31/60)
  df$STANDARD_AST_rpp_acs <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_combos,"S","I")
  print(32/60)
  df$STANDARD_AST_rpp_aas <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_abs,"S","I")
  print(33/60)
  df$STANDARD_AST_rpp_iss <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_singles,"S","I")
  print(34/60)
  df$STANDARD_AST_rpp_ics <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_combos,"S","I")
  print(35/60)
  df$STANDARD_AST_rpp_ais <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_ivs,"S","I")
  print(36/60)
  df$STANDARD_AST_rpp_oss <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_singles,"S","I")
  print(37/60)
  df$STANDARD_AST_rpp_ocs <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_combos,"S","I")
  print(38/60)
  df$STANDARD_AST_rpp_aos <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_orals,"S","I")
  print(39/60)
  df$STANDARD_AST_rpp_acss <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_singles,"S","I")
  print(40/60)
  df$STANDARD_AST_rpp_accs <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_combos,"S","I")
  print(41/60)
  df$STANDARD_AST_rpp_aacs <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_access,"S","I")
  print(42/60)
  df$STANDARD_AST_rpp_wass <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_singles,"S","I")
  print(43/60)
  df$STANDARD_AST_rpp_wacs <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_combos,"S","I")
  print(44/60)
  df$STANDARD_AST_rpp_awas <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_watch,"S","I")
  print(45/60)
  df$STANDARD_AST_rpp_asr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_singles,"R","NT")
  print(46/60)
  df$STANDARD_AST_rpp_acr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_combos,"R","NT")
  print(47/60)
  df$STANDARD_AST_rpp_aar <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_abs,"R","NT")
  print(48/60)
  df$STANDARD_AST_rpp_isr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_singles,"R","NT")
  print(49/60)
  df$STANDARD_AST_rpp_icr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_combos,"R","NT")
  print(50/60)
  df$STANDARD_AST_rpp_air <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_ivs,"R","NT")
  print(51/60)
  df$STANDARD_AST_rpp_osr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_singles,"R","NT")
  print(52/60)
  df$STANDARD_AST_rpp_ocr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_combos,"R","NT")
  print(53/60)
  df$STANDARD_AST_rpp_aor <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_orals,"R","NT")
  print(54/60)
  df$STANDARD_AST_rpp_acsr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_singles,"R","NT")
  print(55/60)
  df$STANDARD_AST_rpp_accr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_combos,"R","NT")
  print(56/60)
  df$STANDARD_AST_rpp_aacr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_access,"R","NT")
  print(57/60)
  df$STANDARD_AST_rpp_wasr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_singles,"R","NT")
  print(58/60)
  df$STANDARD_AST_rpp_wacr <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_combos,"R","NT")
  print(59/60)
  df$STANDARD_AST_rpp_awar <- df %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_watch,"R","NT")
  print(60/60)
  
  df
  
}

###Counting the number of S or I results per antimicrobial (total)
number_ab_results <- function(df,start_col,end_col,which_abs,result_1,result_2) {
  
  all_si <- c()
  start_col <- enquo(start_col)
  end_col <- enquo(end_col)
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. ==result_1) %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%dplyr::slice(i) %>% unlist(),all_abs))) %>% 
      dplyr::slice(i) %>% t() %>% data.frame() %>% filter(. ==result_2) %>% rownames()
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Assembling dot plot dataframe
create_df <- function(df, col_name, aware_result, panel,type,result) {
  result <- df %>%
    select({{ col_name }}) %>%
    cbind(AWaRe_results = aware_result, Panel = panel,Agents=type,Result=result) %>%
    as.data.frame()
  colnames(result) <- c("n", "AWaRe_results", "Panel","Agents","Result")
  return(result)
}
assemble_dotplot_df <- function(df) {
  
  PDRx_rpp_ass <- df %>% create_df(PDRx_rpp_ass,"PDRx\nsingle S","PDRx","Single","S")
  PDRx_rpp_acs <- df %>% create_df(PDRx_rpp_acs,"PDRx\ncombo S","PDRx","Combo","S")
  PDRx_rpp_aas <- df %>% create_df(PDRx_rpp_aas,"PDRx\nall S","PDRx","All","S")
  PDRx_rpp_iss <- df %>% create_df(PDRx_rpp_iss,"PDRx\nsingle IV S","PDRx","Single IV","S")
  PDRx_rpp_ics <- df %>% create_df(PDRx_rpp_ics,"PDRx\ncombo IV S","PDRx","Combo IV","S")
  PDRx_rpp_ais <- df %>% create_df(PDRx_rpp_ais,"PDRx\nall IV S","PDRx","All IV","S")
  PDRx_rpp_oss <- df %>% create_df(PDRx_rpp_oss,"PDRx\nsingle oral S","PDRx","Single Oral","S")
  PDRx_rpp_ocs <- df %>% create_df(PDRx_rpp_ocs,"PDRx\ncombo oral S","PDRx","Combo Oral","S")
  PDRx_rpp_aos <- df %>% create_df(PDRx_rpp_aos,"PDRx\nall oral S","PDRx","All Oral","S")
  PDRx_rpp_acss <- df %>% create_df(PDRx_rpp_acss,"PDRx\nsingle Access S","PDRx","Single Access","S")
  PDRx_rpp_accs <- df %>% create_df(PDRx_rpp_accs,"PDRx\ncombo Access S","PDRx","Combo Access","S")
  PDRx_rpp_aacs <- df %>% create_df(PDRx_rpp_aacs,"PDRx\nall Access S","PDRx","All Access","S")
  PDRx_rpp_wass <- df %>% create_df(PDRx_rpp_wass,"PDRx\nsingle Watch S","PDRx","Single Watch","S")
  PDRx_rpp_wacs <- df %>% create_df(PDRx_rpp_wacs,"PDRx\ncombo watch S","PDRx","Combo Watch","S")
  PDRx_rpp_awas <- df %>% create_df(PDRx_rpp_awas,"PDRx\nall watch S","PDRx","All Watch","S")
  PDRx_rpp_asr <- df %>% create_df(PDRx_rpp_asr,"PDRx\nsingle R","PDRx","Single","R")
  PDRx_rpp_acr <- df %>% create_df(PDRx_rpp_acr,"PDRx\ncombo R","PDRx","Combo","R")
  PDRx_rpp_aar <- df %>% create_df(PDRx_rpp_aar,"PDRx\nall R","PDRx","All","R")
  PDRx_rpp_isr <- df %>% create_df(PDRx_rpp_isr,"PDRx\nsingle IV R","PDRx","Single IV","R")
  PDRx_rpp_icr <- df %>% create_df(PDRx_rpp_icr,"PDRx\ncombo IV R","PDRx","Combo IV","R")
  PDRx_rpp_air <- df %>% create_df(PDRx_rpp_air,"PDRx\nall IV R","PDRx","All IV","R")
  PDRx_rpp_osr <- df %>% create_df(PDRx_rpp_osr,"PDRx\nsingle oral R","PDRx","Single Oral","R")
  PDRx_rpp_ocr <- df %>% create_df(PDRx_rpp_ocr,"PDRx\ncombo oral R","PDRx","Combo Oral","R")
  PDRx_rpp_aor <- df %>% create_df(PDRx_rpp_aor,"PDRx\nall oral R","PDRx","All Oral","R")
  PDRx_rpp_acsr <- df %>% create_df(PDRx_rpp_acsr,"PDRx\nsingle Access R","PDRx","Single Access","R")
  PDRx_rpp_accr <- df %>% create_df(PDRx_rpp_accr,"PDRx\ncombo Access R","PDRx","Combo Access","R")
  PDRx_rpp_aacr <- df %>% create_df(PDRx_rpp_aacr,"PDRx\nall Access R","PDRx","All Access","R")
  PDRx_rpp_wasr <- df %>% create_df(PDRx_rpp_wasr,"PDRx\nsingle Watch R","PDRx","Single Watch","R")
  PDRx_rpp_wacr <- df %>% create_df(PDRx_rpp_wacr,"PDRx\ncombo Watch R","PDRx","Combo Watch","R")
  PDRx_rpp_awar <- df %>% create_df(PDRx_rpp_awar,"PDRx\nall Watch R","PDRx","All Watch","R")
  STANDARD_AST_rpp_ass <- df %>% create_df(STANDARD_AST_rpp_ass,"Standard\nsingle S","Standard","Single","S")
  STANDARD_AST_rpp_acs <- df %>% create_df(STANDARD_AST_rpp_acs,"Standard\ncombo S","Standard","Combo","S")
  STANDARD_AST_rpp_aas <- df %>% create_df(STANDARD_AST_rpp_aas,"Standard\nall S","Standard","All","S")
  STANDARD_AST_rpp_iss <- df %>% create_df(STANDARD_AST_rpp_iss,"Standard\nsingle IV S","Standard","Single IV","S")
  STANDARD_AST_rpp_ics <- df %>% create_df(STANDARD_AST_rpp_ics,"Standard\ncombo IV S","Standard","Combo IV","S")
  STANDARD_AST_rpp_ais <- df %>% create_df(STANDARD_AST_rpp_ais,"Standard\nall IV S","Standard","All IV","S")
  STANDARD_AST_rpp_oss <- df %>% create_df(STANDARD_AST_rpp_oss,"Standard\nsingle oral S","Standard","Single Oral","S")
  STANDARD_AST_rpp_ocs <- df %>% create_df(STANDARD_AST_rpp_ocs,"Standard\ncombo oral S","Standard","Combo Oral","S")
  STANDARD_AST_rpp_aos <- df %>% create_df(STANDARD_AST_rpp_aos,"Standard\nall oral S","Standard","All Oral","S")
  STANDARD_AST_rpp_acss <- df %>% create_df(STANDARD_AST_rpp_acss,"Standard\nsingle Access S","Standard","Single Access","S")
  STANDARD_AST_rpp_accs <- df %>% create_df(STANDARD_AST_rpp_accs,"Standard\ncombo Access S","Standard","Combo Access","S")
  STANDARD_AST_rpp_aacs <- df %>% create_df(STANDARD_AST_rpp_aacs,"Standard\nall Access S","Standard","All Access","S")
  STANDARD_AST_rpp_wass <- df %>% create_df(STANDARD_AST_rpp_wass,"Standard\nsingle Watch S","Standard","Single Watch","S")
  STANDARD_AST_rpp_wacs <- df %>% create_df(STANDARD_AST_rpp_wacs,"Standard\ncombo Watch S","Standard","Combo Watch","S")
  STANDARD_AST_rpp_awas <- df %>% create_df(STANDARD_AST_rpp_awas,"Standard\nall Watch S","Standard","All Watch","S")
  STANDARD_AST_rpp_asr <- df %>% create_df(STANDARD_AST_rpp_asr,"Standard\nsingle R","Standard","Single","R")
  STANDARD_AST_rpp_acr <- df %>% create_df(STANDARD_AST_rpp_acr,"Standard\ncombo R","Standard","Combo","R")
  STANDARD_AST_rpp_aar <- df %>% create_df(STANDARD_AST_rpp_aar,"Standard\nall R","Standard","All","R")
  STANDARD_AST_rpp_isr <- df %>% create_df(STANDARD_AST_rpp_isr,"Standard\nsingle IV R","Standard","Single IV","R")
  STANDARD_AST_rpp_icr <- df %>% create_df(STANDARD_AST_rpp_icr,"Standard\ncombo IV R","Standard","Combo IV","R")
  STANDARD_AST_rpp_air <- df %>% create_df(STANDARD_AST_rpp_air,"Standard\nall IV R","Standard","All IV","R")
  STANDARD_AST_rpp_osr <- df %>% create_df(STANDARD_AST_rpp_osr,"Standard\nsingle oral R","Standard","Single Oral","R")
  STANDARD_AST_rpp_ocr <- df %>% create_df(STANDARD_AST_rpp_ocr,"Standard\ncombo oral R","Standard","Combo Oral","R")
  STANDARD_AST_rpp_aor <- df %>% create_df(STANDARD_AST_rpp_aor,"Standard\nall oral R","Standard","All Oral","R")
  STANDARD_AST_rpp_acsr <- df %>% create_df(STANDARD_AST_rpp_acsr,"Standard\nsingle Access R","Standard","Single Access","R")
  STANDARD_AST_rpp_accr <- df %>% create_df(STANDARD_AST_rpp_accr,"Standard\ncombo Access R","Standard","Combo Access","R")
  STANDARD_AST_rpp_aacr <- df %>% create_df(STANDARD_AST_rpp_aacr,"Standard\nall Access R","Standard","All Access","R")
  STANDARD_AST_rpp_wasr <- df %>% create_df(STANDARD_AST_rpp_wasr,"Standard\nsingle Watch R","Standard","Single Watch","R")
  STANDARD_AST_rpp_wacr <- df %>% create_df(STANDARD_AST_rpp_wacr,"Standard\ncombo Watch R","Standard","Combo Watch","R")
  STANDARD_AST_rpp_awar <- df %>% create_df(STANDARD_AST_rpp_awar,"Standard\nall Watch R","Standard","All Watch","R")
  
  acs_df <- data.frame(rbind(PDRx_rpp_ass,PDRx_rpp_acs,PDRx_rpp_aas,PDRx_rpp_iss,PDRx_rpp_ics,
                             PDRx_rpp_ais,PDRx_rpp_oss,PDRx_rpp_ocs,PDRx_rpp_aos,PDRx_rpp_acss,
                             PDRx_rpp_accs,PDRx_rpp_aacs,PDRx_rpp_wass,PDRx_rpp_wacs,PDRx_rpp_awas,
                             PDRx_rpp_asr,PDRx_rpp_acr,PDRx_rpp_aar,PDRx_rpp_isr,PDRx_rpp_icr,
                             PDRx_rpp_air,PDRx_rpp_osr,PDRx_rpp_ocr,PDRx_rpp_aor,PDRx_rpp_acsr,
                             PDRx_rpp_accr,PDRx_rpp_aacr,PDRx_rpp_wasr,PDRx_rpp_wacr,PDRx_rpp_awar,
                             STANDARD_AST_rpp_ass,STANDARD_AST_rpp_acs,STANDARD_AST_rpp_aas,STANDARD_AST_rpp_iss,
                             STANDARD_AST_rpp_ics,STANDARD_AST_rpp_ais,STANDARD_AST_rpp_oss,STANDARD_AST_rpp_ocs,
                             STANDARD_AST_rpp_aos,STANDARD_AST_rpp_acss,STANDARD_AST_rpp_accs,STANDARD_AST_rpp_aacs,
                             STANDARD_AST_rpp_wass ,STANDARD_AST_rpp_wacs,STANDARD_AST_rpp_awas,STANDARD_AST_rpp_asr,
                             STANDARD_AST_rpp_acr,STANDARD_AST_rpp_aar,STANDARD_AST_rpp_isr,STANDARD_AST_rpp_icr,
                             STANDARD_AST_rpp_air,STANDARD_AST_rpp_osr,STANDARD_AST_rpp_ocr,STANDARD_AST_rpp_aor,
                             STANDARD_AST_rpp_acsr,STANDARD_AST_rpp_accr,STANDARD_AST_rpp_aacr,STANDARD_AST_rpp_wasr,
                             STANDARD_AST_rpp_wacr,STANDARD_AST_rpp_awar))
  
  acs_df <- acs_df %>% group_by(AWaRe_results) %>% mutate(iqr_min=quantile(n)[2],
                                                          iqr_max=quantile(n)[4]) %>% 
    ungroup() %>% 
    mutate(iqr_min=case_when(iqr_min<0 ~ 0,TRUE~iqr_min))
  acs_df <- acs_df %>% rename(Approach="Panel")
  
  summdf <- acs_df %>% group_by(AWaRe_results) %>% mutate(median_n = median(n)) %>%
    ungroup() %>% left_join(acs_df %>% group_by(AWaRe_results) %>% 
                              count(n,name="n_points"),by=c("AWaRe_results","n")) %>% distinct()
  
  acs_df
  
}

###Main dot plot of all S results and WHO Access S results
main_dotplotter <- function(df,PDRx_1,standard_1,PDRx_2,standard_2,
                            left_label,right_label,sens_addendum="") {
  
  plot_df <- df %>% filter(grepl(PDRx_1,AWaRe_results) |
                             grepl(standard_1,AWaRe_results) |
                             grepl(PDRx_2,AWaRe_results) |
                             grepl(standard_2,AWaRe_results) )
  
  plot_df$AWaRe_results <- factor(plot_df$AWaRe_results,levels = c(PDRx_1,
                                                                   standard_1,
                                                                   PDRx_2,
                                                                   standard_2))
  
  plot_df <- plot_df %>% mutate(Approach=case_when(
    Approach=="PDRx"~"ADA",TRUE~Approach
  ))
  
  df_plot <- ggplot(plot_df,aes(x=AWaRe_results,y=n,color=Approach)) +
    geom_jitter(colour="black", alpha=0.01, width=0.1,height=0.15) +
    stat_summary(geom="point",fun="median",size=4)+
    geom_errorbar(aes(ymin=iqr_min,ymax=iqr_max,width=0,color=Approach),show.legend = F)+
    ylab("")+
    ggtitle(glue("Microsimulation study:\nNumber of susceptible results provided per specimen\n{(sens_addendum)}")) +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust=3),
          plot.margin = unit(c(0.1,0.1,0.1,0.25),"inches"),
          plot.title = element_text(hjust=0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank())+
    geom_vline(xintercept = 2.5,linetype="dashed",color="grey") +
    scale_y_continuous(limits = c(-0.15,7),breaks=c(0:6)) +
    scale_color_manual(values=c("#00BFC4","#F8766D"))+
    geom_text(x=1.5,y=6.75,label=glue("{left_label}"),color="#3C3C3C",size=4) +
    geom_text(x=3.5,y=6.75,label=glue("{right_label}"),color="#3C3C3C",size=4)
  
  ggsave(glue("ur_synd_UF_{left_label}_{right_label}_plot.pdf"), plot = df_plot, device = "pdf", width = 6, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  df_plot
  
}

###Determining total number of S or I results provided by personalised panel
number_SorI_pdast <- function(df,which_abs) {
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(PDRx_1:PDRx_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(PDRx_1:PDRx_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                       dplyr::slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Determining total number of S or I results provided by standard panel
number_SorI_standard <- function(df,which_abs) {
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(STANDARD_AST_1,STANDARD_AST_2,STANDARD_AST_3,
                                                           STANDARD_AST_4,STANDARD_AST_5,STANDARD_AST_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                     dplyr::slice(i)) == "S" |
                    (df %>%
                       select(all_of(intersect(df %>% select(STANDARD_AST_1,STANDARD_AST_2,STANDARD_AST_3,
                                                             STANDARD_AST_4,STANDARD_AST_5,STANDARD_AST_6) %>%dplyr::slice(i) %>% unlist(),which_abs))) %>% 
                       dplyr::slice(i)) == "I")
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Cleaning and joining observations information
obs_clean <- function(microsimdf,eddf,triagedf,probsdf,edpxdf,
                      prdfname,mdfname){
  
  staykey <- eddf %>% select(stay_id,intime) %>% distinct(stay_id,.keep_all = T) %>% 
    rename(charttime="intime")
  triagedf <- triagedf %>% left_join(staykey) %>% relocate(charttime,.before = "temperature") %>% 
    mutate(chartdate=as.Date(charttime)) %>% select(subject_id,chartdate,acuity) %>% 
    filter(!is.na(acuity))
  microsimdf <- microsimdf %>% select(-acuity)
  microsimdf <- microsimdf %>% semi_join(triagedf,by=c("subject_id","chartdate"))
  probsdf <- probsdf %>% semi_join(microsimdf,by="micro_specimen_id")
  microsimdf <- microsimdf %>% left_join(triagedf) %>%
    distinct(micro_specimen_id,.keep_all = T)
  acuitykey <- microsimdf %>% select(micro_specimen_id,acuity)
  probsdf <- probsdf %>% left_join(acuitykey)
  
  edpxdf <- MIMER::clean_antibiotics(edpxdf,drug_col=drug)
  
  edpxdf <- edpxdf %>% filter(is_abx) %>% 
    mutate(chartdate=as.Date(charttime))
  
  edpxdf <- edpxdf %>% rename(ab_name="abx_name")
  
  edpxdf <- edpxdf %>% mutate(ab_name=case_when(grepl("Piperacillin",ab_name)~
                                                  "Piperacillin/tazobactam",
                                                TRUE~ab_name))
  edpxdf <- edpxdf %>% semi_join(microsimdf,by=c("chartdate","subject_id"))
  
  edpxdf <- edpxdf %>%
    mutate(order = match(ab_name, all_singles %>% ab_name())) %>%
    distinct(subject_id,chartdate,ab_name,.keep_all = T) %>% 
    arrange(subject_id,chartdate, order) %>%
    group_by(subject_id,chartdate) %>%
    summarize(ab_name = str_c(ab_name, collapse = "_"), .groups = "drop") %>% 
    ungroup() %>% 
    mutate(ab_name=str_replace_all(ab_name,"/","-"))
  
  microsimdf <- microsimdf %>% semi_join(edpxdf,by=c("subject_id","chartdate"))
  edpxdf_key <- edpxdf %>% select(ab_name,subject_id) %>% rename(
    Prescribed_abx = "ab_name"
  )
  microsimdf <- microsimdf %>% left_join(edpxdf_key)
  modified_abx_map <- combined_antimicrobial_map
  names(modified_abx_map) <- str_replace_all(names(modified_abx_map),
                                             " & ", "_")
  microsimdf <- microsimdf %>% filter(Prescribed_abx%in%names(modified_abx_map))
  microsimdf <- microsimdf %>% mutate(Px_Abx = abcombo_replace(Prescribed_abx,modified_abx_map))
  
  probsdf <- probsdf %>% semi_join(microsimdf,by="micro_specimen_id")
  
  assign(prdfname,probsdf,envir=.GlobalEnv)
  assign(mdfname,microsimdf,envir=.GlobalEnv)
  
}

###Line plot of antimicrobial recommendations
ab_lineplot <- function(plotdf,pltitle) {
  
  ggplot(plotdf,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
    
    ###Line plot
    geom_line()+
    
    ###Theme and title
    theme_minimal()+
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )+
    ggtitle(glue("{pltitle}"))+
    
    ###Illness severity on x axis
    scale_x_discrete(breaks = seq(0,3))+
    
    ###Manual colour scale to align antimicrobials
    scale_color_manual(values=collist)+
    
    ###Percentage on y axis
    ylim(0,100)
  
}

##Data upload

util_probs_df <- read_csv("dce_util_probs_df.csv")
ur_util <- read_csv("interim_ur_util.csv")
train_abx <- read_csv("train_abx.csv")
urine_df <- read_csv("pos_urines.csv")
abx_ref <- read_csv("abx_ref_all_combos.csv")
vitalsign <- read_csv("vitalsign.csv")
triage <- read_csv("triage.csv")
edstays <- read_csv("edstays.csv")
pyxis <- read_csv("pyxis.csv") %>% rename(drug="name")
scores <- read_csv("scores_df.csv")

##Antimicrobial mapping lists

###Main antimicrobial lists
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

###Additional mapping lists for specific circumstances
singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%all_singles]
iv_singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%iv_singles]
oral_singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%oral_singles]
modified_abx_map <- combined_antimicrobial_map
names(modified_abx_map) <- str_replace_all(names(modified_abx_map),
                                           " & ", "_")
modified_abx_map_brief <- modified_abx_map[intersect(names(modified_abx_map),(util_probs_df$Antimicrobial))]
iv_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_ivs]
oral_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_orals]
access_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_access]
iv_mod_brief <- intersect(modified_abx_map_brief,iv_combos)
oral_mod_brief <- intersect(modified_abx_map_brief,oral_combos)
access_mod_brief <- intersect(modified_abx_map_brief,all_access)
all_abs <- all_singles
long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
ab_singles <- names(singles_map)
iv_ab_singles <- names(iv_singles_map)
oral_ab_singles <- names(oral_singles_map)

##Observations information

###Cleaning and adding to microsim and prob dfs
ur_util %>% obs_clean(edstays,triage,util_probs_df,pyxis,
                      "util_probs_df","ur_util")

###Writing updated microsim and probability dataframes to CSV
write_csv(util_probs_df,"ur_synd_util_probs_df_pre_util.csv")
write_csv(ur_util,"ur_synd_ur_util_pre_recs.csv")

###Filtering to cases where there were urinary symptoms and/or diagnoses
staykey <- edstays %>% select(stay_id,intime) %>% distinct(stay_id,.keep_all = T) %>% 
  rename(charttime="intime")
triagekey <- triage %>% left_join(staykey) %>% relocate(charttime,.before = "temperature") %>% 
  mutate(chartdate=as.Date(charttime)) %>% select(subject_id,chartdate,acuity,stay_id) %>% 
  filter(!is.na(acuity))
compkey <- triage %>% select(stay_id,chiefcomplaint) %>% 
  group_by(stay_id) %>%
  summarise(chiefcomplaint = paste(chiefcomplaint, collapse = ", "), .groups = "drop") %>% 
  ungroup()
diag_key <- ed_diagnosis %>% select(stay_id,icd_title) %>% 
  group_by(stay_id) %>%
  summarise(icd_title = paste(icd_title, collapse = ", "), .groups = "drop") %>% 
  ungroup()

ed_util <- ur_util %>% select(-c(acuity,chiefcomplaint)) %>% semi_join(triagekey,by=c("subject_id","chartdate")) %>% 
  left_join(triagekey,by=c("subject_id","chartdate"))
ur_util <- ed_util %>% left_join(compkey) %>% left_join(diag_key)

ur_icds <- ur_util %>% filter(grepl("(pyelonephritis|urinary tract infection|tract|cystitis|dysuria|painful micturition|nocturia|urgency|incomplete bladder|bladder pain|frequency of micturition|hematuria|pyuria|abnormal findings in urine|hesitancy|retention of urine|flank pain)",icd_title,
                                    ignore.case = T))

ur_comps <- ur_util %>% filter(grepl("(pyelonephritis|urinary tract infection|tract|cystitis|dysuria|painful micturition|nocturia|urgency|emptying|frequency|hematuria|pyuria|urine|urinary|hesitancy|retention|flank pain)",chiefcomplaint,
                                     ignore.case = T))
ur_icdcomps <- rbind(ur_icds,ur_comps)

ur_util_icdcom <- ur_util %>% semi_join(ur_icds,by="subject_id")
ur_util <- ur_util_icdcom
util_probs_df <- util_probs_df %>% semi_join(ur_util,by="subject_id")
##Utility analysis

###Convert acuity to 0-3 scale
ur_util <- ur_util %>% mutate(acuity=5-acuity) %>% 
  mutate(acuity=acuity-1)
util_probs_df <- util_probs_df %>% mutate(acuity=5-acuity) %>% 
  mutate(acuity=acuity-1)

##Calculate U
util_probs_df <- util_probs_df %>% 
  U_calculation(b = 1)

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
util_probs_df <- util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

###Table of utility values
util_probs_df %>% group_by(Antimicrobial) %>% filter(Antimicrobial%in%names(singles_map)) %>% 
  summarise(Median_util=median(U)) %>% 
  arrange(desc(Median_util))

###Formulary agents (for plot colour coding)
formulary_agents <- c()

###Lists of single-agent, combination, intravenous, and oral variants
uplot1 <- c(rep("U",3),rep("Urosepsis_U",2),rep("Outpatient_U",3))
uplot2 <- c(rep("Treatment",3),rep("Intravenous treatment",2),rep("Oral treatment",3))
uplot3 <- c("","(single agent)"," (combinations)"," (single agent)"," (combinations)",
            ""," (single agent)"," (combinations)")
length(uplot3)

###Iterate over and plot above variants
for (i in seq_along(uplot1)){
  
  util_probs_df %>% utility_plot(!!sym(uplot1[i]),uplot2[i],uplot3[i])
  
}

##Individual treatment recommendations

##Overall
ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_") %>% 
  mutate(across(starts_with("PDRx_"), ~ abcombo_reverse(., combined_antimicrobial_map)))

##IV
ur_util <- ur_util %>% assign_Intravenous(util_probs_df,"Intravenous_") %>% 
  mutate(across(starts_with("Intravenous_"), ~ abcombo_reverse(., combined_antimicrobial_map)))

###Oral
ur_util <- ur_util %>% assign_Oral(util_probs_df,"Oral_") %>% 
  mutate(across(starts_with("Oral_"), ~ abcombo_reverse(., combined_antimicrobial_map)))

##Combinations overall
ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_comb_",names(modified_abx_map_brief)) %>% 
  mutate(across(starts_with("PDRx_comb_"), ~ abcombo_reverse(., modified_abx_map_brief)))

##IV
ur_util <- ur_util %>% assign_Intravenous(util_probs_df,"Intravenous_comb_",names(iv_mod_brief)) %>% 
  mutate(across(starts_with("Intravenous_comb_"), ~ abcombo_reverse(., iv_mod_brief)))

###Oral
ur_util <- ur_util %>% assign_Oral(util_probs_df,"Oral_comb_",names(oral_mod_brief)) %>% 
  mutate(across(starts_with("Oral_comb_"), ~ abcombo_reverse(., oral_mod_brief)))


###Standard panel treatment & AST recommendations
ur_util <- ur_util %>% 
  assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO") %>% 
  mutate(STANDARD_IV_1="TZP",
         STANDARD_IV_2="MEM",
         STANDARD_PO_1="NIT",
         STANDARD_PO_2="SXT") %>% 
  mutate(Px_Abx=abcombo_reverse(Px_Abx,combined_antimicrobial_map))

###Get susceptibility results for recommendations
ur_util <- ur_util %>%
  rowwise() %>%
  mutate(PDIVRx_1_result = get(Intravenous_1),
         PDPORx_1_result = get(Oral_1),
         PDRx_1_result = get(PDRx_1),
         STIVRx_1_result = get(STANDARD_IV_1),
         STIVRx_2_result = get(STANDARD_IV_2),
         STPORx_1_result = get(STANDARD_PO_1),
         STPORx_2_result = get(STANDARD_PO_2),
         Px_Abx_result = get(Px_Abx)) %>%
  ungroup()

###Write dfs with attached recommendations to CSVs
write_csv(util_probs_df,"ur_synd_util_probs_df_final.csv")
write_csv(ur_util,"ur_synd_ur_util_final.csv")

##Chi-squared testing for proportion of isolates susceptible

###N isolates
n_all <- c(nrow(ur_util), nrow(ur_util))

###All antibiotics
x_all <- c(nrow(ur_util %>% filter((PDRx_1_result=="S"|PDRx_1_result=="I"))), nrow(ur_util %>% filter((Px_Abx_result=="S"|Px_Abx_result=="I"))))
perc_all1 <- round(prop.test(x_all, n_all)$estimate[1], 3)
perc_all2 <- round(prop.test(x_all, n_all)$estimate[2], 3)
p2_all <- round(prop.test(x_all, n_all)$p.value, 3)
ci1_all <- round(prop.test(x_all, n_all)$conf.int[1], 3)
ci2_all <- round(prop.test(x_all, n_all)$conf.int[2], 3)
glue("Automated recommendation proportion overall: {perc_all1}
        Actual prescriptions proportion overall: {perc_all2}
        P value: {p2_all}
        95% CI of difference: {ci1_all} to {ci2_all}")

##Access antibiotics
x_acs <- c(nrow(ur_util %>% filter((PDRx_1_result=="S"|PDRx_1_result=="I")&PDRx_1%in%access_singles)), nrow(ur_util %>% filter((Px_Abx_result=="S"|Px_Abx_result=="I")&Px_Abx%in%access_modified_abx_map)))
perc_acs1 <- round(prop.test(x_acs, n_all)$estimate[1], 3)
perc_acs2 <- round(prop.test(x_acs, n_all)$estimate[2], 3)
p2_acs <- round(prop.test(x_acs, n_all)$p.value, 3)
ci1_acs <- round(prop.test(x_acs, n_all)$conf.int[1], 3)
ci2_acs <- round(prop.test(x_acs, n_all)$conf.int[2], 3)
glue("Automated recommendation proportion Access: {perc_acs1}
        Actual prescriptions proportion Access: {perc_acs2}
        P value: {p2_acs}
        95% CI of difference: {ci1_acs} to {ci2_acs}")

###Oral antibiotics
x_orals <- c(nrow(ur_util %>% filter((PDRx_1_result=="S"|PDRx_1_result=="I")&PDRx_1%in%oral_singles)), nrow(ur_util %>% filter((Px_Abx_result=="S"|Px_Abx_result=="I")&Px_Abx%in%oral_modified_abx_map)))
perc_orals1 <- round(prop.test(x_orals, n_all)$estimate[1], 3)
perc_orals2 <- round(prop.test(x_orals, n_all)$estimate[2], 3)
p2_orals <- round(prop.test(x_orals, n_all)$p.value, 3)
ci1_orals <- round(prop.test(x_orals, n_all)$conf.int[1], 3)
ci2_orals <- round(prop.test(x_orals, n_all)$conf.int[2], 3)
glue("Automated recommendation proportion orals: {perc_orals1}
        Actual prescriptions proportion orals: {perc_orals2}
        P value: {p2_orals}
        95% CI of difference: {ci1_orals} to {ci2_orals}")

###IV antibiotics
x_ivs <- c(nrow(ur_util %>% filter((PDRx_1_result=="S"|PDRx_1_result=="I")&PDRx_1%in%iv_singles)), nrow(ur_util %>% filter((Px_Abx_result=="S"|Px_Abx_result=="I")&Px_Abx%in%iv_modified_abx_map)))
perc_ivs1 <- round(prop.test(x_ivs, n_all)$estimate[1], 3)
perc_ivs2 <- round(prop.test(x_ivs, n_all)$estimate[2], 3)
p2_ivs <- round(prop.test(x_ivs, n_all)$p.value, 3)
ci1_ivs <- round(prop.test(x_ivs, n_all)$conf.int[1], 3)
ci2_ivs <- round(prop.test(x_ivs, n_all)$conf.int[2], 3)
glue("Automated recommendation proportion IVs: {perc_ivs1}
        Actual prescriptions proportion IVs: {perc_ivs2}
        P value: {p2_ivs}
        95% CI of difference: {ci1_ivs} to {ci2_ivs}")


##Illness severity analysis

###Set empty vectors
iv_perclist <- c()
po_perclist <- c()
overall_perclist <- c()
ivr_perclist <- c()
por_perclist <- c()
overallr_perclist <- c()
ivac_list <- c()
poac_list <- c()
ovac_list <- c()
ovor_list <- c()
oviv_list <- c()
ovacr_list <- c()
ovorr_list <- c()
ovivr_list <- c()
pdrx1_list <- data.frame(PDRx_1=as.ab(ab_singles))
percs_list <- data.frame(PDRx_1=as.ab(ab_singles))
cam <- combined_antimicrobial_map %>% unlist()
names(cam) <- NULL
overall_perclist_px_abx <- c()
overallr_perclist_px_abx <- c()
ovac_list_px_abx <- c()
ovor_list_px_abx <- c()
oviv_list_px_abx <- c()
ovacr_list_px_abx <- c()
ovorr_list_px_abx <- c()
ovivr_list_px_abx <- c()
pdrx1_list_px_abx <- data.frame(Px_Abx=cam)
percs_list_px_abx <- data.frame(Px_Abx=cam)
weightseq <- seq(0,3)

###Iterate along illness severity scores
for(weight in seq_along(weightseq)) {
  
  ###Check results at illness severity iteration value
  result_checker(ur_util,util_probs_df,weightseq[weight])
  
  ###Populate vectors
  iv_perclist <- c(iv_perclist,iv_perc)
  po_perclist <- c(po_perclist,po_perc)
  overall_perclist <- c(overall_perclist,overall_perc)
  ivr_perclist <- c(ivr_perclist,ivr_perc)
  por_perclist <- c(por_perclist,por_perc)
  overallr_perclist <- c(overallr_perclist,overallr_perc)
  ivac_list <- c(ivac_list,iv_s_access)
  poac_list <- c(poac_list,po_s_access)
  ovac_list <- c(ovac_list,overall_s_access)
  ovor_list <- c(ovor_list,overall_s_oral)
  oviv_list <- c(oviv_list,overall_s_iv)
  ovacr_list <- c(ovacr_list,overall_r_access)
  ovorr_list <- c(ovorr_list,overall_r_oral)
  ovivr_list <- c(ovivr_list,overall_r_iv)
  pdrx1_list <- pdrx1_list %>% left_join(abrx1_df,by="PDRx_1")
  colnames(pdrx1_list)[weightseq[weight]+2] <- weightseq[weight]
  percs_list <- percs_list %>% left_join(abrx1_percs,by="PDRx_1")
  colnames(percs_list)[weightseq[weight]+2] <- weightseq[weight]
  
  overall_perclist_px_abx <- c(overall_perclist_px_abx,overall_perc_px_abx)
  overallr_perclist_px_abx <- c(overallr_perclist_px_abx,overallr_perc_px_abx)
  ovac_list_px_abx <- c(ovac_list_px_abx,overall_s_access_px_abx)
  ovor_list_px_abx <- c(ovor_list_px_abx,overall_s_oral_px_abx)
  oviv_list_px_abx <- c(oviv_list_px_abx,overall_s_iv_px_abx)
  ovacr_list_px_abx <- c(ovacr_list_px_abx,overall_r_access_px_abx)
  ovorr_list_px_abx <- c(ovorr_list_px_abx,overall_r_oral_px_abx)
  ovivr_list_px_abx <- c(ovivr_list_px_abx,overall_r_iv_px_abx)
  pdrx1_list_px_abx <- pdrx1_list_px_abx %>% left_join(abrx1_df_px_abx,by="Px_Abx")
  colnames(pdrx1_list_px_abx)[weightseq[weight]+2] <- weightseq[weight]
  percs_list_px_abx <- percs_list_px_abx %>% left_join(abrx1_percs_px_abx,by="Px_Abx")
  colnames(percs_list_px_abx)[weightseq[weight]+2] <- weightseq[weight]
  
}

###Bind treatment type labels for AUF recommendations
iv_perclist <- iv_perclist %>% label_binder("All agents")
po_perclist <- po_perclist %>% label_binder("All agents")
overall_perclist <- overall_perclist %>% label_binder("All agents")
ivr_perclist <- ivr_perclist %>% label_binder("All agents")
por_perclist <- por_perclist %>% label_binder("All agents")
overallr_perclist <- overallr_perclist %>% label_binder("All agents")
ivac_list <- ivac_list %>% label_binder("Access agents")
poac_list <- poac_list %>% label_binder("Access agents")
ovac_list <- ovac_list %>% label_binder("Access agents")
ovor_list <- ovor_list %>% label_binder("Oral agents")
oviv_list <- oviv_list %>% label_binder("IV agents")
ovacr_list <- ovacr_list %>% label_binder("Access agents")
ovorr_list <- ovorr_list %>% label_binder("Oral agents")
ovivr_list <- ovivr_list %>% label_binder("IV agents")

###Rename percentage columns for AUF recommendations
iv_perclist <- iv_perclist %>% rename(Percentage = "vec")
po_perclist <- po_perclist %>% rename(Percentage = "vec")
overall_perclist <- overall_perclist %>% rename(Percentage = "vec")
ivr_perclist <- ivr_perclist %>% rename(Percentage = "vec")
por_perclist <- por_perclist %>% rename(Percentage = "vec")
overallr_perclist <- overallr_perclist %>% rename(Percentage = "vec")
ivac_list <- ivac_list %>% rename(Percentage = "vec")
poac_list <- poac_list %>% rename(Percentage = "vec")
ovac_list <- ovac_list %>% rename(Percentage = "vec")
ovor_list <- ovor_list %>% rename(Percentage = "vec")
oviv_list <- oviv_list %>% rename(Percentage = "vec")
ovacr_list <- ovacr_list %>% rename(Percentage = "vec")
ovorr_list <- ovorr_list %>% rename(Percentage = "vec")
ovivr_list <- ovivr_list %>% rename(Percentage = "vec")

###Bind treatment type labels for antimicrobials actually prescribed
overall_perclist_px_abx <- overall_perclist_px_abx %>% label_binder("All agents")
overallr_perclist_px_abx <- overallr_perclist_px_abx %>% label_binder("All agents")
ovac_list_px_abx <- ovac_list_px_abx %>% label_binder("Access agents")
ovor_list_px_abx <- ovor_list_px_abx %>% label_binder("Oral agents")
oviv_list_px_abx <- oviv_list_px_abx %>% label_binder("IV agents")
ovacr_list_px_abx <- ovacr_list_px_abx %>% label_binder("Access agents")
ovorr_list_px_abx <- ovorr_list_px_abx %>% label_binder("Oral agents")
ovivr_list_px_abx <- ovivr_list_px_abx %>% label_binder("IV agents")

###Rename percentage columns for antimicrobials actually prescribed
overall_perclist_px_abx <- overall_perclist_px_abx %>% rename(Percentage = "vec")
overallr_perclist_px_abx <- overallr_perclist_px_abx %>% rename(Percentage = "vec")
ovac_list_px_abx <- ovac_list_px_abx %>% rename(Percentage = "vec")
ovor_list_px_abx <- ovor_list_px_abx %>% rename(Percentage = "vec")
oviv_list_px_abx <- oviv_list_px_abx %>% rename(Percentage = "vec")
ovacr_list_px_abx <- ovacr_list_px_abx %>% rename(Percentage = "vec")
ovorr_list_px_abx <- ovorr_list_px_abx %>% rename(Percentage = "vec")
ovivr_list_px_abx <- ovivr_list_px_abx %>% rename(Percentage = "vec")

###Compile dataframes from vectors
iv_xg_plot_df <- data.frame(rbind(
  iv_perclist,ivac_list
))
po_xg_plot_df <- data.frame(rbind(
  po_perclist,poac_list
))
overall_xg_plot_df <- data.frame(rbind(
  overall_perclist,ovac_list,ovor_list,oviv_list
))
overall_xg_plot_df_px_abx <- data.frame(rbind(
  overall_perclist_px_abx,ovac_list_px_abx,ovor_list_px_abx,oviv_list_px_abx
))
overallr_xg_plot_df <- data.frame(rbind(
  overallr_perclist,ovacr_list,ovorr_list,ovivr_list
))
overallr_xg_plot_df_px_abx <- data.frame(rbind(
  overallr_perclist_px_abx,ovacr_list_px_abx,ovorr_list_px_abx,ovivr_list_px_abx
))

##Plots for illness severity analysis  

###Replace NAs in AUF df with zero i.e., no cases of recommendation
pdrx1_list[is.na(pdrx1_list)] <- 0
percs_list[is.na(percs_list)] <- 0

###Remove antimicrobials that were never recommended
pdrx1_df <- percs_list %>% filter(rowSums(select(.,2:ncol(percs_list)))!=0)

###Long form df for plotting
pdrx1_df <- melt(pdrx1_df)
colnames(pdrx1_df) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df <- pdrx1_df %>% mutate(Antimicrobial = ab_name(Antimicrobial))

###Write plot df to csv
write_csv(pdrx1_df,"ur_synd_abplot_df.csv")

###Replace NAs in human prescription df with 0
pdrx1_list_px_abx[is.na(pdrx1_list_px_abx)] <- 0
percs_list_px_abx[is.na(percs_list_px_abx)] <- 0

###Reclassify all antimicrobials with <15 prescriptions across all severities
pdrx1_df_px_abx <- percs_list_px_abx %>% mutate(Px_Abx=case_when((`0`+`1`+`2`+`3`)<15~"Other",
                                                                 TRUE~Px_Abx))

###Remove antimicrobials with no prescriptions
pdrx1_df_px_abx <- pdrx1_df_px_abx %>% filter(rowSums(select(.,2:ncol(pdrx1_df_px_abx)))!=0)

###Long-form df for plot
pdrx1_df_px_abx <- melt(pdrx1_df_px_abx)
colnames(pdrx1_df_px_abx) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")

###Clean up antimicrobial names for plot
pdrx1_df_px_abx <- pdrx1_df_px_abx %>% mutate(Antimicrobial=
                                                abcombo_replace(Antimicrobial,
                                                                combined_antimicrobial_map)) %>% 
  mutate(Antimicrobial=str_replace_all(Antimicrobial,"-","/")) %>% 
  mutate(Antimicrobial=str_replace_all(Antimicrobial,"_"," & "))

###Get percentage of recommendations per antimicrobial
pdrx1_df_px_abx <- pdrx1_df_px_abx %>% group_by(Antimicrobial,`Illness severity`) %>% 
  summarise(`Percentage of first-line recommendations`=sum(`Percentage of first-line recommendations`)) %>% 
  ungroup()

###Write plot df to csv
write_csv(pdrx1_df_px_abx,"ur_synd_abplot_df_px_abx.csv")

###Standardise antimicrobial colours for plots
abcollist <- rbind(pdrx1_df %>% filter(`Illness severity`==1) %>%
                     arrange(desc(`Percentage of first-line recommendations`)) %>% 
                     distinct(Antimicrobial),
                   pdrx1_df_px_abx %>% filter(`Illness severity`==1) %>%
                     arrange(desc(`Percentage of first-line recommendations`))
                   %>% distinct(Antimicrobial)) %>%
  distinct(Antimicrobial) %>% unlist()
collist <-  hue_pal()(length(abcollist))
names(collist) <- abcollist

###Factorise to order by frequency at severity 0
pdrx1_df$Antimicrobial <- factor(pdrx1_df$Antimicrobial,
                                 levels=pdrx1_df %>% filter(`Illness severity`==0) %>%
                                   arrange(desc(`Percentage of first-line recommendations`)) %>% 
                                   pull(Antimicrobial))
pdrx1_df_px_abx$Antimicrobial <- factor(pdrx1_df_px_abx$Antimicrobial,
                                        levels=pdrx1_df_px_abx %>% filter(`Illness severity`==0) %>%
                                          arrange(desc(`Percentage of first-line recommendations`)) %>% 
                                          pull(Antimicrobial))

###Line plot of AUF antimicrobial recommendations
abplot <- pdrx1_df %>% ab_lineplot(
  "ADA antibiotic choices according to illness severity\n(coded UTI diagnoses only)")
abplot

###Save plot to pdf
ggsave(glue("ur_synd_illness_abplot.pdf"), plot = abplot, device = "pdf", width = 10, height = 3,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Line plot of human antimicrobial prescriptions
abplot_px_abx <- pdrx1_df_px_abx %>% ab_lineplot(
  "Human antibiotic choices in ED according to illness severity\n(coded UTI diagnoses only)")
abplot_px_abx

###Save to pdf
ggsave(glue("ur_synd_illness_abplot_px_abx.pdf"), plot = abplot_px_abx, device = "pdf", width = 10, height = 3,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Proportion plot dataframes
write_csv(iv_xg_plot_df,"ur_synd_iv_xg_plot_df.csv")
write_csv(po_xg_plot_df,"ur_synd_po_xg_plot_df.csv")
write_csv(overall_xg_plot_df,"ur_synd_overall_xg_plot_df.csv")

###Barplots of urinary pathogen coverage
overall_xg_plot_df$Type <- "AUF"
overall_xg_plot_df_px_abx$Type <- "Human"
overall_bar_plot(overall_xg_plot_df,overall_xg_plot_df_px_abx,"(oral|iv)",
                 "Access agents","Watch agents",4,
                 "AWaRe category (in patients with coded UTI diagnoses)")
overall_bar_plot(overall_xg_plot_df,overall_xg_plot_df_px_abx,"(Access|iv)",
                 "Oral agents","Non-oral agents",4,
                 "Oral administrability (in patients with coded UTI diagnoses)")
overall_bar_plot(overall_xg_plot_df,overall_xg_plot_df_px_abx,"(Access|oral)",
                 "IV agents","Non-IV agents",4,
                 "IV administrability (in patients with coded UTI diagnoses)")

###Barplots of urinary pathogen resistance
overallr_xg_plot_df$Type <- "AUF"
overallr_xg_plot_df_px_abx$Type <- "Human"
overallr_bar_plot(overallr_xg_plot_df,overallr_xg_plot_df_px_abx,"(oral|iv)",
                  "Access agents","Watch agents",4,
                  "AWaRe category (in patients with coded UTI diagnoses)")
overallr_bar_plot(overallr_xg_plot_df,overallr_xg_plot_df_px_abx,"(Access|iv)",
                  "Oral agents","Non-oral agents",4,
                  "Oral administrability (in patients with coded UTI diagnoses)")
overallr_bar_plot(overallr_xg_plot_df,overallr_xg_plot_df_px_abx,"(Access|oral)",
                  "IV agents","Non-IV agents",4,
                  "IV administrability (in patients with coded UTI diagnoses)")

##AST performance analysis - number of results per panel

###Number of AST results per panel
ur_util <- ur_util %>% rpp_ast()

###Assemble data frame for dot plot data visualisation
acs_df <- ur_util %>% assemble_dotplot_df()
acs_df %>% count(AWaRe_results)
write_csv(acs_df,"ur_synd_uf_ast_sourcedata_aware_dotplot.csv")

###Dot plot of number of all S results and Access S results per panel
main_aware_plot <- acs_df %>% main_dotplotter("PDRx\nsingle S","Standard\nsingle S","PDRx\nsingle Access S","Standard\nsingle Access S",
                                              "All agents","WHO access agents","(in patients with coded UTI diagnoses)")

###Dot plot of number of oral S results and IV S results per panel
main_ivo_plot <- acs_df %>% main_dotplotter("PDRx\nsingle IV S","Standard\nsingle IV S","PDRx\nsingle oral S","Standard\nsingle oral S",
                                            "IV agents","Oral agents","(in patients with coded UTI diagnoses)")

###Total number of S or I results provided by personalised panel
ur_util$n_allS_PDRx6 <- ur_util %>% number_SorI_pdast(all_abs)

###Number of Access category S or I results provided by personalised panel
access_abs <- all_access
ur_util$n_acS_PDRx6 <- ur_util %>% number_SorI_pdast(access_abs)

###Number of oral S or I results provided by personalised panel
oral_antibs <- c("AMP","SAM","CIP",
                 "SXT","NIT")
ur_util$n_poS_PDRx6 <- ur_util %>% number_SorI_pdast(oral_antibs)

###Number of IV S or I results provided by personalised panel
iv_antibs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
               "GEN","SXT","VAN")
ur_util$n_ivS_PDRx6 <- ur_util %>% number_SorI_pdast(iv_antibs)

###Total number of S or I results provided by standard panel
ur_util$n_allS_standard6 <- ur_util %>% number_SorI_standard(all_abs)

###Number of Access category S or I results provided by standard panel
ur_util$n_acS_standard6 <- ur_util %>% number_SorI_standard(access_abs)

###Number of oral S or I results provided by standard panel
ur_util$n_poS_standard6 <- ur_util %>% number_SorI_standard(oral_antibs)

###Number of IV S or I results provided by standard panel
ur_util$n_ivS_standard6 <- ur_util %>% number_SorI_standard(iv_antibs)

###Wilxocon signed ranks test of AST results per panel
p_acs <- round(wilcox.test(ur_util %>% pull(n_acS_PDRx6), ur_util %>% pull(n_acS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_all <- round(wilcox.test(ur_util %>% pull(n_allS_PDRx6), ur_util %>% pull(n_allS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_pos <- round(wilcox.test(ur_util %>% pull(n_poS_PDRx6), ur_util %>% pull(n_poS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_ivs <- round(wilcox.test(ur_util %>% pull(n_ivS_PDRx6), ur_util %>% pull(n_ivS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)

cil_acs <- round(wilcox.test(ur_util %>% pull(n_acS_PDRx6), ur_util %>% pull(n_acS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[1], 3)
ciu_acs <- round(wilcox.test(ur_util %>% pull(n_acS_PDRx6), ur_util %>% pull(n_acS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[2], 3)
cil_all <- round(wilcox.test(ur_util %>% pull(n_allS_PDRx6), ur_util %>% pull(n_allS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[1], 3)
ciu_all <- round(wilcox.test(ur_util %>% pull(n_allS_PDRx6), ur_util %>% pull(n_allS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[2], 3)
cil_pos <- round(wilcox.test(ur_util %>% pull(n_poS_PDRx6), ur_util %>% pull(n_poS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[1], 3)
ciu_pos <- round(wilcox.test(ur_util %>% pull(n_poS_PDRx6), ur_util %>% pull(n_poS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[2], 3)
cil_ivs <- round(wilcox.test(ur_util %>% pull(n_ivS_PDRx6), ur_util %>% pull(n_ivS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[1], 3)
ciu_ivs <- round(wilcox.test(ur_util %>% pull(n_ivS_PDRx6), ur_util %>% pull(n_ivS_standard6), paired = TRUE, conf.int = TRUE)$conf.int[2], 3)

z_access <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(ur_util$n_acS_PDRx6), " ~ ", quo_text(ur_util$n_acS_standard6))), data = ur_util), "standardized")[1]
acs_effectsize <- round(z_access / sqrt(nrow(ur_util)), 3)

z_all <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(ur_util$n_allS_PDRx6), " ~ ", quo_text(ur_util$n_allS_standard6))), data = ur_util), "standardized")[1]
all_effectsize <- round(z_all / sqrt(nrow(ur_util)), 3)

z_po <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(ur_util$n_poS_PDRx6), " ~ ", quo_text(ur_util$n_poS_standard6))), data = ur_util), "standardized")[1]
po_effectsize <- round(z_po / sqrt(nrow(ur_util)), 3)

z_iv <- statistic(wilcoxsign_test(as.formula(paste0(quo_text(ur_util$n_ivS_PDRx6), " ~ ", quo_text(ur_util$n_ivS_standard6))), data = ur_util), "standardized")[1]
iv_effectsize <- round(z_iv / sqrt(nrow(ur_util)), 3)

glue("All results:
        Median {median(ur_util$n_allS_PDRx6)} overall S/I results per specimen (IQR {quantile(ur_util$n_allS_PDRx6)[2]} to {quantile(ur_util$n_allS_PDRx6)[4]})
        Median {median(ur_util$n_allS_standard6)} overall S/I results per specimen (IQR {quantile(ur_util$n_allS_standard6)[2]} to {quantile(ur_util$n_allS_standard6)[4]})
        p {ifelse(p_all<0.001,'<0.001',p_all)}, effect size {all_effectsize}
        
        Access results:
        Median {median(ur_util$n_acS_PDRx6)} Access S/I results per specimen (IQR {quantile(ur_util$n_acS_PDRx6)[2]} to {quantile(ur_util$n_acS_PDRx6)[4]})
        Median {median(ur_util$n_acS_standard6)} Access S/I results per specimen (IQR {quantile(ur_util$n_acS_standard6)[2]} to {quantile(ur_util$n_acS_standard6)[4]})
        p {ifelse(p_acs<0.001,'<0.001',p_acs)}, effect size {acs_effectsize}
       
       Oral results:
        Median {median(ur_util$n_poS_PDRx6)} Oral S/I results per specimen (IQR {quantile(ur_util$n_poS_PDRx6)[2]} to {quantile(ur_util$n_poS_PDRx6)[4]})
        Median {median(ur_util$n_poS_standard6)} Oral S/I results per specimen (IQR {quantile(ur_util$n_poS_standard6)[2]} to {quantile(ur_util$n_poS_standard6)[4]})
        p {ifelse(p_pos<0.001,'<0.001',p_pos)}, effect size {po_effectsize}
       
       IV results:
        Median {median(ur_util$n_ivS_PDRx6)} IV S/I results per specimen (IQR {quantile(ur_util$n_ivS_PDRx6)[2]} to {quantile(ur_util$n_ivS_PDRx6)[4]})
        Median {median(ur_util$n_ivS_standard6)} IV S/I results per specimen (IQR {quantile(ur_util$n_ivS_standard6)[2]} to {quantile(ur_util$n_ivS_standard6)[4]})
        p {ifelse(p_ivs<0.001,'<0.001',p_ivs)}, effect size {iv_effectsize}")

###Write final manuscript dfs to csv
write_csv(util_probs_df,"ur_synd_util_probs_df_manuscript.csv")
write_csv(ur_util,"ur_synd_ur_util_manuscript.csv")
