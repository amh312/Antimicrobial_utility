#SPECIALTY SUBANALYSIS

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
}

###Utility score calculator
abs_calc <- function(val,prob) {
  
  ifelse(val>=0,val*prob,abs(val)*(1-prob))
  
}

###Utility data visualisation
utility_plot <- function(df, variable,application,modification="") {
  
  variable <- enquo(variable)
  formulary_agents <- c()
  
  axiscols <- if_else(
    df %>% pull(Antimicrobial) %in% formulary_agents,
    "seagreen", "black")
  
  df <- df %>% mutate(Antimicrobial = str_replace_all(Antimicrobial,"_"," & "))
  
  if (application=="Intravenous treatment") {
    
    df <- df %>% filter(Urosepsis_U != min(Urosepsis_U))
    
  } else if (application=="Oral treatment") {
    
    df <- df %>% filter(Outpatient_U != min(Outpatient_U))
    
  } else if (application=="AST") {
    
    df <- df %>% filter(AST_utility != min(AST_utility) & single_agent)
    
  }
  
  if (grepl("single",modification,ignore.case=T)) {
    
    df <- df %>% filter(single_agent)
    
  } else if (grepl("combination",modification,ignore.case=T)) {
    
    df <- df %>% filter(!single_agent)
    
  }
  
  thisplot <- ggplot(df %>% mutate(Antimicrobial = 
                                     factor(Antimicrobial,
                                            levels = df %>% group_by(Antimicrobial) %>% 
                                              summarise(Median_util=median(!!variable)) %>% 
                                              arrange(Median_util) %>% select(Antimicrobial) %>% unlist())), 
                     aes(x=!!variable,y=Antimicrobial,fill=Antimicrobial)) +
    geom_boxplot(outlier.color = NULL,
                 outlier.alpha = 0.3) +
    theme_minimal() +
    ggtitle(glue("{application} utility distribution on\nmicrosimulation dataset{modification}")) +
    theme(legend.position = "None",axis.text.y = element_text(
      colour = axiscols))+
    xlab(glue("{application} utility"))+
    theme()
  
  ggsave(glue("utility_{application}_{modification}.pdf"), plot = thisplot, device = "pdf", width = 10, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(thisplot)
  
}

###Separated-specialty sensitivity analysis
spec_sens_analysis <- function(df,rankings_df,specialty) {
  
  ###Engineer scores dataframe
  scores <- rankings_df %>% dce_cleanup(characteristics)
  
  ###Mlogit object for ranked logit
  dce_logdata <- mlogit.data(scores,
                             choice = "Rank", 
                             shape = "long", 
                             chid.var = "id", 
                             alt.var = "Antibiotic", 
                             ranked = TRUE)
  
  ###Ranked logit formula
  dce_formula <- Rank ~ CDI_highrisk + Toxicity_highrisk +  
    Oral_option + UTI_specific + IV_option + High_cost + Access + Reserve | 0
  
  ###Fit model
  dce_model <- mlogit(dce_formula, data = dce_logdata, method = "bfgs")
  
  ###Extract coefficients from model
  coefficients <- coef(dce_model)
  
  ###Add coefficients to scores dataframe
  scores <- data.frame(
    Coefficient = names(coefficients),
    Value = coefficients,
    OR = exp(coefficients),
    OR_dif = exp(coefficients) - 1
  ) %>% 
    mutate(stan_OR = (OR - 1) / max(abs(OR - 1))) %>% 
    mutate(colour = case_when(
      Value > 0 ~ "B", 
      Value < 0 ~ "A"
    )) %>% 
    mutate(Coefficient = c("High CDI risk","High toxicity risk","Oral option",
                           "UTI-specific","IV option","High cost","Access category",
                           "Reserve category")) %>% 
    arrange(Value)
  score_rownames <- rownames(scores)
  
  ###Make empty dataframe for ceofficient bootstrap
  bstr_dcecoefs <- matrix(NA, nrow = 1000, ncol = length(coef(mlogit(dce_formula, data = dce_logdata, method = "bfgs"))))
  
  ###Iterate over 1000 seeds
  for (i in 1:1000) {
    set.seed(i)
    
    ###Index to sample with replacement
    bstr_ind <- sample(unique(dce_logdata$id), replace = TRUE)
    
    ###Subset by index
    bstr_ind_data <- dce_logdata %>% subset(id %in% bstr_ind)
    
    ###Fit model
    bstr_mod <- mlogit(dce_formula, data = bstr_ind_data, method = "bfgs")
    
    ###Get coefficients
    bstr_dcecoefs[i, ] <- coef(bstr_mod)
  }
  
  ###Make dataframe with blank iq ranges
  iqs <- data.frame(Coefficient=factor(names(coef(bstr_mod))),Q5=NA,Q95=NA)
  
  ###Add labels to iq dataframe
  iqs <- iqs %>% mutate(Coefficient=case_when(
    grepl("CDI",Coefficient)~ "High CDI risk",
    grepl("Toxicity",Coefficient)~ "High toxicity risk",
    grepl("UTI",Coefficient)~ "UTI-specific",
    grepl("Access",Coefficient)~ "Access category",
    grepl("Reserve",Coefficient)~ "Reserve category",
    TRUE~Coefficient
  )) %>% 
    mutate(Coefficient=str_replace_all(Coefficient,"_"," "))
  
  ###Add IQ ranges of bootstrapped coefficients to dataframe
  for (i in 1:ncol(bstr_dcecoefs)) {
    
    iqs$Q5[i] <- quantile(bstr_dcecoefs[,i],probs=seq(0,1,0.05))['5%']
    iqs$Q95[i] <- quantile(bstr_dcecoefs[,i],probs=seq(0,1,0.05))['95%']
    
  }
  
  ###Join IQ ranges to DCE dataframe
  scores <- scores %>% left_join(iqs)
  rownames(scores) <- score_rownames
  
  ###Write cleaned DCE df to csv
  write_csv(scores,"scores_df.csv")
  
  ###Factorise rankings for plot
  scores$Coefficient <- factor(scores$Coefficient, levels=
                                 scores %>% arrange(Value) %>% 
                                 select(Coefficient) %>% unlist())
  
  ###Plot DCE rankings
  sens_features <- ggplot(scores,aes(x=Value,y=Coefficient,fill=colour)) +
    
    ###Bar chart
    geom_col() +
    
    ###Cross line
    geom_hline(aes(yintercept=0)) +
    
    ###Titles and labels
    ylab("Drug property") +
    xlab("Coefficient value for drug selection probability") +
    ggtitle(glue("The effect of different antimicrobial drug properties on\n{specialty} clinician prescribing preference in the UTI\nscenario discrete choice experiment"))+
    
    ###Central line
    geom_vline(xintercept = 0,colour="grey40")+
    
    ###Confidence intervals
    geom_errorbar(aes(xmin = Q5, xmax = Q95), width = 0.1,colour="grey40")+
    
    ###Theme
    theme_minimal()+
    theme(panel.grid = element_blank(),
          legend.position = "None")
  
  ggsave(glue("{specialty}_importances.pdf"), plot = sens_features, device = "pdf", width = 8, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  print(sens_features)
  
  
  ###Antimicrobial dummy variables in probability prediction dataframe
  df$abx_name_ <- as.factor(df$Antimicrobial)
  df <- df %>% mutate(
    abx_name_ = str_replace_all(abx_name_,"-",".")
  )
  dummy_vars <- model.matrix(~ abx_name_ - 1, data = df)
  df <- cbind(df, dummy_vars) %>% tibble() %>% 
    select(-abx_name_)
  
  ##Utility score calculation
  
  df <- df %>%
    utility_preprocess(scores,cost_list,ur_util)
  
  ###Calculate overall utility score
  df <- df %>% U_calculation(b = 1)
  
  ##Utility analysis
  
  ###Utility data visualisation
  df %>% group_by(Antimicrobial) %>% 
    summarise(Median_util=median(U)) %>% 
    arrange(desc(Median_util))
  formulary_agents <- c()
  
  df
  
}

###Preprocessing probability/utility dataframe
utility_preprocess <- function(probsdf,scoredf,costdf,urref_df){
  
  #cdi risk utility
  cdi_value <- scoredf[rownames(scoredf)=="CDI_highrisk",] %>% 
    select(Value) %>% unlist()
  
  #toxicity risk utility
  tox_value <- scoredf[rownames(scoredf)=="Toxicity_highrisk",] %>% 
    select(Value) %>% unlist()
  
  #nitrofurantoin as only uti-specific agent
  uti_specifics <- c("Nitrofurantoin")
  
  #uti-specific utility
  uti_value <- scoredf[rownames(scoredf)=="UTI_specific",] %>% 
    select(Value) %>% unlist()
  
  #list of access agents
  access_abs <- c("AMP","SAM","CZO",
                  "GEN","SXT","NIT") %>% ab_name() %>% 
    str_replace("/","-")
  
  #access agent combinations
  access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  
  #put singles and combos together
  access_abs <- c(access_abs, access_combos)
  
  #access agent utility
  access_value <- scoredf[rownames(scoredf)=="Access",] %>% 
    select(Value) %>% unlist()
  
  #vector of oral options
  oral_abs <- c("AMP","SAM","CIP",
                "SXT","NIT") %>% ab_name() %>% 
    str_replace("/","-")
  
  #oral combinations
  oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  
  #put singles and combos together
  oral_abs <- c(oral_abs, oral_combos)
  
  #oral agent utility
  oral_value <- scoredf[rownames(scoredf)=="Oral_option",] %>% 
    select(Value) %>% unlist()
  
  #vector of iv agents
  iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
              "GEN","SXT","VAN") %>% ab_name() %>% 
    str_replace("/","-") 
  
  #iv combinations
  iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  
  #put singles and combos together
  iv_abs <- c(iv_abs, iv_combos)
  
  #iv agent utility
  iv_value <- scoredf[rownames(scoredf)=="IV_option",] %>% 
    select(Value) %>% unlist()
  
  #no reserve category antibiotics in this study
  reserve_abs <- c()
  
  #reserve class utility
  reserve_value <- scoredf[rownames(scoredf)=="Reserve",] %>% 
    select(Value) %>% unlist()
  
  #highcost utility
  cost_value <- scoredf[rownames(scoredf)=="High_cost",] %>% 
    select(Value) %>% unlist()
  
  #get drugs from cost list
  drug_order <- costdf %>% distinct(`Generic Name`) %>% unlist()
  
  #get costs
  costdf <- costdf %>% group_by(`Generic Name`) %>% summarise(orig_cost=min(Cost)) %>% ungroup() %>% 
    mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
  
  #get agent combinations
  comb <- combn(costdf$`Generic Name`, 2, simplify = FALSE)
  
  #sum costs for combinations
  cost_comb <- combn(costdf$orig_cost, 2, function(x) sum(x))
  
  #combinations into df
  costdf2 <- data.frame(
    Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
    orig_cost = cost_comb
  )  
  
  #binding singles and combos together
  costdf <- costdf %>% rename(Antimicrobial="Generic Name")
  costdf <- data.frame(rbind(costdf,costdf2))
  
  #normalise minimum costs
  costdf <- costdf %>% mutate(min_cost = orig_cost/max(orig_cost),
                              Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
  
  #join costs to probability df
  probsdf <- probsdf %>% left_join(costdf)
  
  #check medians and iqrs of cost
  costdf %>% dplyr::slice(1:13) %>% summarise(
    MED = Median(orig_cost),
    iq1 = quantile(orig_cost)[2],
    iq3 = quantile(orig_cost)[4]
  )
  
  #attach individual utilities to dataframe
  probsdf %>% 
    
    #cdi
    mutate(value_CDI = cdi_value,
           util_CDI = abs_calc(value_CDI,prob_CDI),
           
           #toxicity
           value_tox = tox_value,
           util_tox = abs_calc(value_tox,prob_tox),
           
           #uti-specificity
           UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
           value_UTI = uti_value,
           util_uti = abs_calc(value_UTI,UTI_specific),
           
           #access agent
           Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
           value_access = access_value,
           util_access = abs_calc(value_access,Access_agent),
           
           #oral agent
           Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
           value_oral = oral_value,
           util_oral = abs_calc(value_oral,Oral_agent),
           
           #iv agent
           IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
           value_iv = iv_value,
           util_iv = abs_calc(value_iv,IV_agent),
           
           #reserve agent
           Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
           value_reserve = reserve_value,
           util_reserve = abs_calc(value_reserve,Reserve_agent),
           
           #cost
           Highcost_agent = min_cost,
           value_highcost = cost_value,
           util_highcost = abs_calc(value_highcost,Highcost_agent),
           
           #single agent
           single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE)
    )
  
}

###Utility calculation
U_calculation <- function(df,b=1) {
    
    df %>%
      
      #calculate U
      mutate(U=S*(util_tox+util_CDI+util_uti+util_access+util_oral+
                    util_reserve+util_highcost+(util_iv*exp(acuity)*b))) %>% 
      
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

###Utility plot
utilplot_iterate <- function(df,filterterm){
  
  rxvar <- c("U","Urosepsis_U","Outpatient_U")
  rxroute <- c("Treatment","Intravenous treatment","Oral treatment")
  spectypes <- c("infection specialties","intensive care","medical specialties",
                 "surgical specialties","general practice")
  rxtype <- c(","," (single agent,"," (combinations,")
  suffix <- c("",")",")")
  specroutelist <- c()
  
  for (i in seq_along(spectypes)){
    
    for (j in seq_along(rxtype)){
      
      specroutelist <- c(specroutelist,
                         glue("{rxtype[j]} {spectypes[i]}{suffix[j]}"))
      
    }
    
  }
  
  specroutelist <- specroutelist[grepl(filterterm,specroutelist)]
  
  for (i in seq_along(rxroute)){
    
    for (j in seq_along(specroutelist)){
      
      df %>% utility_plot(!!sym(rxvar[i]),
                          rxroute[i],
                          specroutelist[j])
      
    }
    
  }
  
}

##Read-in and filter
rankings <- read_csv("labelled_DCE_results.csv")
util_probs_df_2 <- read_csv("probs_df_overall.csv")
ref_util<- read_csv("util_probs_df_final.csv")

###Join treatment to dataframe
acuities <- ref_util %>% select(micro_specimen_id,acuity) %>% distinct(
  micro_specimen_id,.keep_all = T)

util_probs_df_2 <- util_probs_df_2 %>% semi_join(ref_util,by="micro_specimen_id") %>% 
  left_join(acuities,by="micro_specimen_id")

###Read-in and filter questionnaire dataframes
rankings_infection <- rankings %>% filter(Specialty=="Infection") %>% select(-Specialty)
rankings_itu <- rankings %>% filter(Specialty=="Intensive care") %>% select(-Specialty)
rankings_medicine <- rankings %>% filter(Specialty=="Medicine") %>% select(-Specialty)
rankings_surgery <- rankings %>% filter(Specialty=="Surgery") %>% select(-Specialty)
rankings_gp <- rankings %>% filter(Specialty=="General Practice") %>% select(-Specialty)

###Examine coefficients and apply to utility function
util_infection <- util_probs_df_2 %>% spec_sens_analysis(rankings_infection,"Infection")
util_itu <- util_probs_df_2 %>% spec_sens_analysis(rankings_itu,"Intensive care")
util_medicine <- util_probs_df_2 %>% spec_sens_analysis(rankings_medicine,"Medicine")
util_surgery <- util_probs_df_2 %>% spec_sens_analysis(rankings_surgery,"Surgery")
util_gp <- util_probs_df_2 %>% spec_sens_analysis(rankings_gp,"General Practice")

###Examine utility function output for each group
util_infection %>% utilplot_iterate("infection")
util_itu %>% utilplot_iterate("intensive care")
util_medicine %>% utilplot_iterate("medical")
util_surgery %>% utilplot_iterate("surgical")
util_gp %>% utilplot_iterate("general practice")
