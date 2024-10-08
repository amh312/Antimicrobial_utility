#UTILITY CALCULATIONS

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                 overall_tox = factor(overall_tox),
                 sepsis_ae=factor(sepsis_ae))
}

###Utility score calculator
calculate_utilities <- function(df,formulary_list=NULL) {
  
  df <- df %>% mutate(overall_util = util_uti + util_access +
                                       util_oral + util_iv +
                                       util_reserve + util_highcost,
                      overall_oral_util = util_uti + util_access +
                        util_oral +
                        util_reserve + util_highcost,
                      overall_iv_util = util_uti + util_access +
                          util_iv +
                          util_reserve + util_highcost,
                      R_penalty=R*prob_sepsisae,
                S_utility = S*overall_util,
                S_PO_utility = S*overall_oral_util,
                S_IV_utility=S*overall_iv_util,
                ent_S_utility = ent_S*overall_util,
                AMPR_utility = case_when(
                  Antimicrobial=="Ampicillin" ~
                    AMP_R_value*R, TRUE~0
                ),
                SAMR_utility = case_when(
                  Antimicrobial=="Ampicillin-sulbactam" ~
                    SAM_R_value*R, TRUE~0
                ),
                TZPR_utility = case_when(
                  Antimicrobial=="Piperacillin-tazobactam" ~
                    TZP_R_value*R, TRUE~0
                ),
                CZOR_utility = case_when(
                  Antimicrobial=="Cefazolin" ~
                    CZO_R_value*R, TRUE~0
                ),
                CROR_utility = case_when(
                  Antimicrobial=="Ceftriaxone" ~
                    CRO_R_value*R, TRUE~0
                ),
                CAZR_utility = case_when(
                  Antimicrobial=="Ceftazidime" ~
                    CAZ_R_value*R, TRUE~0
                ),
                FEPR_utility = case_when(
                  Antimicrobial=="Cefepime" ~
                    FEP_R_value*R, TRUE~0
                ),
                MEMR_utility = case_when(
                  Antimicrobial=="Meropenem" ~
                    MEM_R_value*R, TRUE~0
                ),
                CIPR_utility = case_when(
                  Antimicrobial=="Ciprofloxacin" ~
                    CIP_R_value*R, TRUE~0
                ),
                GENR_utility = case_when(
                  Antimicrobial=="Gentamicin" ~
                    GEN_R_value*R, TRUE~0
                ),
                SXTR_utility = case_when(
                  Antimicrobial=="Trimethoprim-sulfamethoxazole" ~
                    SXT_R_value*R, TRUE~0
                ),
                NITR_utility = case_when(
                  Antimicrobial=="Nitrofurantoin" ~
                    NIT_R_value*R, TRUE~0
                ),
                VANR_utility = case_when(
                  Antimicrobial=="Vancomycin" ~
                    VAN_R_value*R, TRUE~0
                ),
                Formulary_agent = case_when(
                  Antimicrobial%in%formulary_list ~
                    TRUE, TRUE~FALSE
                ),
                Formulary_utility = 
                  Formulary_agent*R,
                AST_utility = (S_utility+
                  AMPR_utility +
                  SAMR_utility +
                  TZPR_utility +
                  CZOR_utility +
                  CROR_utility +
                  CAZR_utility +
                  FEPR_utility +
                  MEMR_utility +
                  CIPR_utility +
                  GENR_utility +
                  SXTR_utility +
                  NITR_utility +
                  VANR_utility +
                  Formulary_utility)*single_agent,
                Rx_utility = S_utility -
                  R_penalty,
                ent_AMPR_utility = case_when(
                  Antimicrobial=="Ampicillin" ~
                    AMP_R_value*ent_R, TRUE~0
                ),
                ent_SAMR_utility = case_when(
                  Antimicrobial=="Ampicillin-sulbactam" ~
                    SAM_R_value*ent_R, TRUE~0
                ),
                ent_TZPR_utility = case_when(
                  Antimicrobial=="Piperacillin-tazobactam" ~
                    TZP_R_value*ent_R, TRUE~0
                ),
                ent_CZOR_utility = case_when(
                  Antimicrobial=="Cefazolin" ~
                    CZO_R_value*ent_R, TRUE~0
                ),
                ent_CROR_utility = case_when(
                  Antimicrobial=="Ceftriaxone" ~
                    CRO_R_value*ent_R, TRUE~0
                ),
                ent_CAZR_utility = case_when(
                  Antimicrobial=="Ceftazidime" ~
                    CAZ_R_value*ent_R, TRUE~0
                ),
                ent_FEPR_utility = case_when(
                  Antimicrobial=="Cefepime" ~
                    FEP_R_value*ent_R, TRUE~0
                ),
                ent_MEMR_utility = case_when(
                  Antimicrobial=="Meropenem" ~
                    MEM_R_value*ent_R, TRUE~0
                ),
                ent_CIPR_utility = case_when(
                  Antimicrobial=="Ciprofloxacin" ~
                    CIP_R_value*ent_R, TRUE~0
                ),
                ent_GENR_utility = case_when(
                  Antimicrobial=="Gentamicin" ~
                    GEN_R_value*ent_R, TRUE~0
                ),
                ent_SXTR_utility = case_when(
                  Antimicrobial=="Trimethoprim-sulfamethoxazole" ~
                    SXT_R_value*ent_R, TRUE~0
                ),
                ent_NITR_utility = case_when(
                  Antimicrobial=="Nitrofurantoin" ~
                    NIT_R_value*ent_R, TRUE~0
                ),
                ent_VANR_utility = case_when(
                  Antimicrobial=="Vancomycin" ~
                    VAN_R_value*ent_R, TRUE~0
                ),
                ent_Formulary_agent = case_when(
                  Antimicrobial%in%formulary_list ~
                    TRUE, TRUE~FALSE
                ),
                ent_Formulary_utility = 
                  Formulary_agent*ent_R,
                ent_AST_utility = (ent_S_utility+
                                 ent_AMPR_utility +
                                 ent_SAMR_utility +
                                 ent_TZPR_utility +
                                 ent_CZOR_utility +
                                 ent_CROR_utility +
                                 ent_CAZR_utility +
                                 ent_FEPR_utility +
                                 ent_MEMR_utility +
                                 ent_CIPR_utility +
                                 ent_GENR_utility +
                                 ent_SXTR_utility +
                                 ent_NITR_utility +
                                 ent_VANR_utility +
                                 ent_Formulary_utility)*single_agent,
                ent_Rx_utility = (overall_util * ent_S) -
                  (ent_R*(prob_sepsisae - util_CDI -
                     util_tox)))
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~S_IV_utility - R_penalty),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~S_PO_utility -  R_penalty),
      ent_Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~(overall_iv_util*ent_S)-(ent_R*(prob_sepsisae - util_CDI -
                                                           util_tox))),
      ent_Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~(overall_oral_util*ent_S)-(ent_R*(prob_sepsisae - util_CDI -
                                                            util_tox)))
    )
  
}

###Utility data visualisation
utility_plot <- function(df, variable,application,modification="") {
  
  variable <- enquo(variable)
  
  axiscols <- if_else(
    df %>% pull(Antimicrobial) %in% formulary_agents,
    "seagreen", "black")
  
  df <- df %>% mutate(Antimicrobial = str_replace_all(Antimicrobial,"_"," & "))
  
  if (application=="Intravenous treatment") {
    
    df <- df %>% filter(Urosepsis_Rx_utility != min(Urosepsis_Rx_utility))
    
  } else if (application=="Oral treatment") {
    
    df <- df %>% filter(Outpatient_Rx_utility != min(Outpatient_Rx_utility))
    
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
  
  ggsave(glue("utility_{application}.pdf"), plot = thisplot, device = "pdf", width = 10, height = 8,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(thisplot)
  
}

###Sensitivity analysis probability distribution check
dens_check <- function(df,abx_agent,result,alpha,beta) {
  
  result <- enquo(result)
  
  test_dist <- tibble(Probability=rbeta(nrow(ur_util),alpha,beta),Distribution="Test")
  actual_dist <- tibble(
    Probability=df %>% filter(Antimicrobial==abx_agent) %>% pull(!!result),
    Distribution="Actual")
  dist_df <- tibble(rbind(test_dist,actual_dist))
  
  thisplot <- ggplot(dist_df,
                     aes(x=Probability,color=Distribution)) + 
    geom_density()+
    ggtitle(glue("{abx_agent} resistance probability"))+
    xlim(0,1)+
    ylab("") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  print(thisplot)
  
}

###Replacing resistance probability distribution with simulated distribution
dist_replace <- function(df,df2,abx_agent,result,result2,alpha,beta) {
  
  replacement <- rbeta(nrow(df2),alpha,beta)
  
  dfrows <- which(df$Antimicrobial== abx_agent)
  dfcols <- which(colnames(df)==result)
  
  dfrows2 <- which(df$Antimicrobial== abx_agent)
  dfcols2 <- which(colnames(df)==result2)
  
  df[dfrows,dfcols] <- replacement
  df[dfrows2,dfcols2] <- 1-replacement
  
  df
  
}

###Sensitivity analysis varying resistance probabilities
dens_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  sens_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(sens_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","r_prob","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- 0.5
    b <- 18.5
    
    probs_df_2 <- probs_df
    probs_df_3 <- probs_df
    
    for(i in 1:9) {
      
      densy <- probs_df_2 %>% dens_check(iterabs[j],R,a,b)
      print(densy)
      
      sens_df <- probs_df_3 %>% dist_replace(df,iterabs[j],"R","S",a,b) %>%
        calculate_utilities()
      
      overall_median <- round(median(sens_df %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- sens_df %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      sens_row <- sens_df %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  r_prob=median(R),
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      sens_cum <- data.frame(rbind(sens_cum,sens_row))
      
      a <- a+i
      b <- b-(i/2)
      
    }
    
  }
  
  sens_cum
  
}

###Data visualisation of resistance probability sensitivity analysis
dens_sens_plot <- function(df,measure,uf) {
  
  uf <- enquo(uf)
  
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  
  for(i in seq_along(iterabs)) {
  
  df_spec_plot <- df %>% filter(Antimicrobial==iterabs[i])
  med_plot <- df_spec_plot %>% mutate(Antimicrobial="All antimicrobials",
                                      med_util=overall_med,
                                      lower_iqr=overall_med,
                                      upper_iqr=overall_med)
  best_plot <- df_spec_plot %>% mutate(Antimicrobial="Best antimicrobial",
                                      med_util=best_ut,
                                      lower_iqr=best_ut,
                                      upper_iqr=best_ut)
  df_spec_plot <- data.frame(rbind(df_spec_plot,med_plot,best_plot))
  
  df_spec_plot$Antimicrobial <- factor(df_spec_plot$Antimicrobial,
                                       levels=c("All antimicrobials","Best antimicrobial",
                                                iterabs[i]))
  
  df_plot <- ggplot(df_spec_plot, aes(x = r_prob)) +
    geom_line(aes(y = as.numeric(med_util), group = Antimicrobial, color = Antimicrobial)) +
    geom_ribbon(aes(y = as.numeric(med_util),
                    ymin = as.numeric(lower_iqr),
                    ymax = as.numeric(upper_iqr),
                    group = Antimicrobial, fill = Antimicrobial), alpha = 0.3) +
    ylim(-1.5,1.5) +
    xlim(0,1) +
    ggtitle(glue("Effect of varying population resistance probability on {iterabs[i]} {measure} utility")) +
    xlab(glue("Population median probability of {iterabs[i]} resistance")) +
    ylab(glue("{measure} utility")) +
    theme_minimal()
  
  ggsave(glue("dens_sens_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 10, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(df_plot)
  
  }
  
}

###Sensitivity analysis varying characteristic importance
uti_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = a,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}
access_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = a,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}
oral_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = a,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}
iv_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = a,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}
cdi_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = a,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}
tox_util_sens <- function(df,probs_df,uf) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- -1
    probs_df_2 <- probs_df
    
    for(i in 1:9) {
      
      cdi_util_probs <- predict(underlying_cdi, probs_df_2,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_2,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_2,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = a,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_2 <- probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_2 <- probs_df_2 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_2 <- form_probs_df_2 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_2 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_2 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=a,
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
      a <- a+0.25
      
    }
    
  }
  
  util_cum
  
}

###Data visualisation of importance sensitivity analysis
util_sens_plot <- function(df,measure,value) {
  
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  
  for(i in seq_along(iterabs)) {
    
    df_spec_plot <- df %>% filter(Antimicrobial==iterabs[i])
    med_plot <- df_spec_plot %>% mutate(Antimicrobial="All antimicrobials",
                                        med_util=overall_med,
                                        lower_iqr=overall_med,
                                        upper_iqr=overall_med)
    best_plot <- df_spec_plot %>% mutate(Antimicrobial="Best antimicrobial",
                                         med_util=best_ut,
                                         lower_iqr=best_ut,
                                         upper_iqr=best_ut)
    df_spec_plot <- data.frame(rbind(df_spec_plot,med_plot,best_plot))
    
    df_spec_plot$Antimicrobial <- factor(df_spec_plot$Antimicrobial,
                                         levels=c("All antimicrobials","Best antimicrobial",
                                                  iterabs[i]))
    
    df_plot <- ggplot(df_spec_plot, aes(x = spec_value)) +
      geom_line(aes(y = as.numeric(med_util), group = Antimicrobial, color = Antimicrobial)) +
      geom_ribbon(aes(y = as.numeric(med_util),
                      ymin = as.numeric(lower_iqr),
                      ymax = as.numeric(upper_iqr),
                      group = Antimicrobial, fill = Antimicrobial), alpha = 0.3) +
      ylim(-1.5,1.5) +
      xlim(-1,1) +
      ggtitle(glue("Effect of varying {value} importance on {iterabs[i]} {measure} utility")) +
      xlab(glue("Importance of {value}")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()
    
    ggsave(glue("{value}_sens_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 10, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
}

###Sensitivity analysis probability distribution check
dens_check_2 <- function(df,abx_agent,result,alpha,beta,characteristic) {
  
  result <- enquo(result)
  
  test_dist <- tibble(Probability=rbeta(nrow(ur_util),alpha,beta),Distribution="Test")
  actual_dist <- tibble(
    Probability=df %>% filter(Antimicrobial==abx_agent) %>% pull(!!result),
    Distribution="Actual")
  dist_df <- tibble(rbind(test_dist,actual_dist))
  
  thisplot <- ggplot(dist_df,
                     aes(x=Probability,color=Distribution)) + 
    geom_density()+
    ggtitle(glue("{abx_agent} {characteristic} probability"))+
    xlim(0,1)+
    ylab("") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  print(thisplot)
  
}

###Replacing other probability distributions with simulated distribution
dist_replace_2 <- function(df,df2,abx_agent,parameter,alpha,beta) {
  
  replacement <- rbeta(nrow(df2),alpha,beta)
  
  dfrows <- which(df$Antimicrobial== abx_agent)
  dfcols <- which(colnames(df)==parameter)
  
  df[dfrows,dfcols] <- replacement
  
  df
  
}

###Sensitivity analysis varying characteristic probabilities
cdi_prob_sens <- function(df,probs_df,uf,characteristic,characteristic_col,char_col,char_col2,characteristic_name) {
  
  uf <- enquo(uf)
  char_col <- enquo(char_col)
  char_col2 <- enquo(char_col2)
  
  sens_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(sens_cum) <- c("med_util","best_util","lower_iqr","upper_iqr",characteristic,"overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- 50
    b <- 1850
    
    probs_df_2 <- probs_df
    probs_df_3 <- probs_df
    
    for(i in 1:9) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_util_probs <- predict(underlying_cdi, probs_df_3,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_3,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_3,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      prob_vector <- probs_df_3 %>% dist_replace_2(df,iterabs[j],characteristic_col,a,b) %>% 
        pull(!!char_col2)
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(prob_CDI = prob_vector,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_3 <- probs_df_3 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_3 <- probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_3 <- probs_df_3 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_3 <- form_probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_3 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_3 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      sens_row <- probs_df_3 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  !!char_col:=median(!!char_col2),
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      sens_cum <- data.frame(rbind(sens_cum,sens_row))
      
      a <- a+(i*100)
      b <- b-((i*100)/2)
      
    }
    
  }
  
  sens_cum
  
}
tox_prob_sens <- function(df,probs_df,uf,characteristic,characteristic_col,char_col,char_col2,characteristic_name) {
  
  uf <- enquo(uf)
  char_col <- enquo(char_col)
  char_col2 <- enquo(char_col2)
  
  sens_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(sens_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr",characteristic,"overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- 50
    b <- 1850
    
    probs_df_2 <- probs_df
    probs_df_3 <- probs_df
    
    for(i in 1:9) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_util_probs <- predict(underlying_cdi, probs_df_3,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_3,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_3,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      prob_vector <- probs_df_3 %>% dist_replace_2(df,iterabs[j],characteristic_col,a,b) %>% 
        pull(!!char_col2)
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = prob_vector,
               value_tox = tox_value,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = sepsis_util_probs,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_3 <- probs_df_3 %>% calculate_utilities()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_3 <- probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_3 <- probs_df_3 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_3 <- form_probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_3 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_3 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      sens_row <- probs_df_3 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  !!char_col:=median(!!char_col2),
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      sens_cum <- data.frame(rbind(sens_cum,sens_row))
      
      a <- a+(i*100)
      b <- b-((i*100)/2)
      
    }
    
  }
  
  sens_cum
  
}
sepsisae_prob_sens <- function(df,probs_df,uf,characteristic,characteristic_col,char_col,char_col2,characteristic_name) {
  
  uf <- enquo(uf)
  char_col <- enquo(char_col)
  char_col2 <- enquo(char_col2)
  
  sens_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(sens_cum) <- c("med_util","best_ut","lower_iqr","upper_iqr",characteristic,"overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    a <- 50
    b <- 1850
    
    probs_df_2 <- probs_df
    probs_df_3 <- probs_df
    
    for(i in 1:9) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_util_probs <- predict(underlying_cdi, probs_df_3,type="response")
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_util_probs <- predict(underlying_tox, probs_df_3,type="response")
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ####Sepsis adverse outcome risk utility
      sepsis_util_probs <- predict(underlying_sepsis, probs_df_3,type="response")
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores[rownames(scores)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores[rownames(scores)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores[rownames(scores)=="High_cost",] %>% 
        select(Value) %>% unlist()
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Enterococcus removal sensitivity analysis key
      ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
        rename(ent_R = "R", ent_S = "S")
      
      prob_vector <- probs_df_3 %>% dist_replace_2(df,iterabs[j],characteristic_col,a,b) %>% 
        pull(!!char_col2)
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(prob_CDI = cdi_util_probs,
               value_CDI = cdi_value,
               util_CDI = prob_CDI * value_CDI,
               prob_tox = tox_util_probs,
               value_tox = tox_util_probs,
               util_tox = prob_tox * value_tox,
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = UTI_specific * value_UTI,
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = Access_agent * value_access,
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = Oral_agent * value_oral,
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = IV_agent * value_iv,
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = Reserve_agent * value_reserve,
               Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
               value_highcost = cost_value,
               util_highcost = Highcost_agent * value_highcost,
               prob_sepsisae = prob_vector,
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_3 <- probs_df_3 %>% calculate_utilities()
      
      print(iterabs[j])
      probs_df_3 %>% 
        group_by(Antimicrobial) %>% summarise(meds=median(Rx_utility)) %>% 
        print()
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
        str_replace_all("/","-")
      probs_df_3 <- probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      ###Dataframe for formulary agent sensitivity analysis
      form_probs_df_3 <- probs_df_3 %>% calculate_utilities(
        formulary_list = c("Ceftriaxone","Ciprofloxacin"))
      form_probs_df_3 <- form_probs_df_3 %>% filter(Antimicrobial %in% abx_in_train)
      
      overall_median <- round(median(probs_df_3 %>% select(!!uf) %>% unlist() %>%  as.numeric()),3)
      
      best_util <- probs_df_3 %>% group_by(Antimicrobial) %>% summarise(
        abmeds = median(!!uf)) %>% ungroup() %>% 
        summarise(bestab = max(abmeds)) %>% unlist()
      
      sens_row <- probs_df_3 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  !!char_col:=median(!!char_col2),
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      sens_cum <- data.frame(rbind(sens_cum,sens_row))
      
      a <- a+(i*100)
      b <- b-((i*100)/2)
      
    }
    
  }
  
  sens_cum
  
}

###Data visualisation of characteristic probability sensitivity analysis
dens_sens_plot_2 <- function(df,characteristic,measure,char_col) {
  
  char_col <- enquo(char_col)
  
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  
  for(i in seq_along(iterabs)) {
    
    df_spec_plot <- df %>% filter(Antimicrobial==iterabs[i])
    med_plot <- df_spec_plot %>% mutate(Antimicrobial="All antimicrobials",
                                        med_util=overall_med,
                                        lower_iqr=overall_med,
                                        upper_iqr=overall_med)
    best_plot <- df_spec_plot %>% mutate(Antimicrobial="Best antimicrobial",
                                         med_util=best_ut,
                                         lower_iqr=best_ut,
                                         upper_iqr=best_ut)
    df_spec_plot <- data.frame(rbind(df_spec_plot,med_plot,best_plot))
    
    df_spec_plot$Antimicrobial <- factor(df_spec_plot$Antimicrobial,
                                         levels=c("All antimicrobials","Best antimicrobial",
                                                  iterabs[i]))
    
    df_plot <- ggplot(df_spec_plot, aes(x = !!char_col)) +
      geom_line(aes(y = as.numeric(med_util), group = Antimicrobial, color = Antimicrobial)) +
      geom_ribbon(aes(y = as.numeric(med_util),
                      ymin = as.numeric(lower_iqr),
                      ymax = as.numeric(upper_iqr),
                      group = Antimicrobial, fill = Antimicrobial), alpha = 0.3) +
      ylim(-1.5,1.5) +
      xlim(0,1) +
      ggtitle(glue("Effect of varying {iterabs[i]} {characteristic} risk on {iterabs[i]} {measure} utility")) +
      xlab(glue("Population median probability of {characteristic} following {iterabs[i]}")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()
    
    ggsave(glue("dens_sens_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 10, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
}

###Separated-specialty sensitivity analysis
spec_sens_analysis <- function(df,results_df,specialty) {
###Engineer scores dataframe
rankings <- results_df %>% select(10:ncol(results_df)) %>% slice(2:nrow(results_df))
colnames(rankings) <- c("Abelfenide","Acetemran","Adenomadin","Adrevenac",
                        "Amrodine","Choriotroban","Cormide","Decloxone",
                        "Dexaset","Endoxolol","Olanzasys","Pansolid",
                        "Protestryl")
Characteristic <- c("AWaRe","CDI","Toxicity","UTI","Oral","IV","Cost")
Abelfenide <- c("Reserve","Low","Low","No","Yes","Yes","High")
Acetemran <- c("Watch","Low","Low","Yes","Yes","No","Low")
Adenomadin <- c("Access","High","Low","No","Yes","Yes","Low")
Adrevenac <- c("Watch","High","Low","No","No","Yes","High")
Amrodine <- c("Reserve","High","Low","No","No","Yes","High")
Choriotroban <- c("Access","Low","High","No","No","Yes","Low")
Cormide <- c("Watch","Low","Low","No","No","Yes","High")
Decloxone <- c("Watch","Low","High","No","No","Yes","High")
Dexaset <- c("Access","Low","Low","No","Yes","Yes","Low")
Endoxolol <- c("Access","Low","Low","Yes","Yes","No","Low")
Olanzasys <- c("Watch","High","Low","No","No","Yes","Low")
Pansolid <- c("Access","Low","Low","No","Yes","No","Low")
Protestryl <- c("Watch","High","Low","No","Yes","Yes","Low")
characteristics <- data.frame(Abelfenide,Acetemran,Adenomadin,
                              Adrevenac,Amrodine,Choriotroban,
                              Cormide,Decloxone,Dexaset,
                              Endoxolol,Olanzasys,Pansolid,Protestryl) %>% 
  t() %>% data.frame()

colnames(characteristics) <- Characteristic
characteristics <- characteristics %>% 
  mutate(Antibiotic = rownames(characteristics)) %>% relocate(
    Antibiotic,.before = AWaRe
  )
rownames(characteristics) <- NULL
scores <- rankings %>% pivot_longer(Abelfenide:Protestryl)
colnames(scores) <- c("Antibiotic","Rank")
scores <- scores %>% left_join(characteristics,by="Antibiotic")

scores <- scores %>% mutate(Access = case_when(AWaRe=="Access"~1,TRUE~0),
                            Watch = case_when(AWaRe=="Watch"~1,TRUE~0),
                            Reserve = case_when(AWaRe=="Reserve"~1,TRUE~0))

scores <- scores %>% rename(CDI_highrisk = "CDI",
                            Toxicity_highrisk = "Toxicity",
                            UTI_specific = "UTI",
                            Oral_option = "Oral",
                            IV_option = "IV",
                            High_cost = "Cost")
scores[scores=="High"] <- "1"
scores[scores=="Low"] <- "0"
scores[scores=="Yes"] <- "1"
scores[scores=="No"] <- "0"

scores

scores <- scores %>% select(-AWaRe) %>%
  mutate(across(c(Rank, CDI_highrisk, Toxicity_highrisk, UTI_specific,Access,
                  Watch,Reserve), as.numeric)) 

repeat_val <- nrow(scores)/13
seq_rep <- seq(1,nrow(scores)/13)
seq_rep1 <- c()

for (i in seq_rep) {
  
  seq_rep2 <- rep(seq_rep[i],13)
  seq_rep1 <- append(seq_rep1,seq_rep2)
  
  
}

scores <- scores %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                            id = seq_rep1) %>% 
  mutate(across(Oral_option:High_cost,as.numeric))

##Extracting characteristic importance weights

##Apply ridge regression
mlogit_data <- mlogit.data(scores, choice = "Rank", shape = "long", 
                           chid.var = "id", alt.var = "Antibiotic", 
                           ranked = TRUE)

formula_no_int <- Rank ~ CDI_highrisk + Toxicity_highrisk +  
  Oral_option + UTI_specific + IV_option + High_cost + Access + Reserve

X <- model.matrix(formula_no_int, data = mlogit_data)
y <- mlogit_data$Rank
fit <- cv.glmnet(X, y, family = "multinomial", alpha = 0)

lambda_min <- fit$lambda.min
lambda_1se <- fit$lambda.1se
cat("Lambda.min:", lambda_min, "\n")
cat("Lambda.1se:", lambda_1se, "\n")

tmp_coeffs <- coef(fit, s = "lambda.min")
scores <- tmp_coeffs$`TRUE` %>% as.matrix() %>% data.frame()
colnames(scores) <- "Value"
scores <- scores %>% mutate(Coefficient = rownames(scores),
                            OR = exp(Value),
                            OR_dif = OR-1)
scores <- scores %>% mutate(colour = case_when(
  Value > 0 ~ "B", Value < 0 ~ "A"
)) %>% slice(-1) %>% slice(-1)

scores <- scores %>% mutate(stan_OR = (OR-1)/max(abs(OR-1)))
scores$Coefficient <- c("High CDI risk","High toxicity risk","Oral option",
                        "UTI-specific","IV option","High cost","Access category",
                        "Reserve category")

scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())

###Visualise odds ratios for antimicrobial characteristics
sens_features <- ggplot(scores,aes(x=OR_dif,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Odds ratio for drug selection") +
  ggtitle(glue("The effect of different antimicrobial drug properties\non {specialty} clinician prescribing preference in UTI scenario"))+
  scale_x_continuous(labels = function(x) x+1)+
  geom_vline(xintercept = 0)

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

##CDI prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
abx_columns <- grep("^abx_name_", names(train_abx), value = TRUE)
service_columns <- grep("^curr_service_", names(train_abx), value = TRUE)
predictors <- c("pCDI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                abx_columns, service_columns, "pICU", "pSEPSIS")
predictors <- predictors[predictors != "abx_name_Ampicillin"]
formula <- reformulate(predictors, response = "CDI")
cdi_fit <- log_reg_spec %>%
  fit(formula,data = train_abx)

underlying_cdi <- extract_fit_engine(cdi_fit)
coef(underlying_cdi) %>% round(2)

###Evaluate on test data
cdi_test_probs <- predict(underlying_cdi, test_abx,type="response")

###Attach to probability prediction dataframe
cdi_util_key <- ur_util %>% select(micro_specimen_id,pHADM,MALE,pICU,
                                   CDI:pSEPSIS) %>% 
  select(-AKI)
df <- df %>% 
  left_join(cdi_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

##Toxicity prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
predictors <- c("prAKI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                abx_columns, service_columns, "pICU", "pSEPSIS","admission_sepsis",
                "SIRS","highCRP")
predictors <- predictors[predictors != "abx_name_Ampicillin"]
formula <- reformulate(predictors, response = "overall_tox")
tox_fit <- log_reg_spec %>%
  fit(formula,
      data = train_abx)

underlying_tox <- extract_fit_engine(tox_fit)
coef(underlying_tox) %>% round(2)

###Evaluate model on test data
tox_test_probs <- predict(underlying_tox, test_abx,type="response")

###Attach to probability prediction dataframe
tox_util_key <- ur_util %>% select(micro_specimen_id,overall_tox,
                                   admission_sepsis,
                                   SIRS,highCRP)
df <- df %>% 
  left_join(tox_util_key,by="micro_specimen_id")

##Sepsis adverse outcomes prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
predictors <- c("prAKI","pCDI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                abx_columns, service_columns, "pICU", "pSEPSIS")
predictors <- predictors[predictors != "abx_name_Ampicillin"]
formula <- reformulate(predictors, response = "sepsis_ae")
sepsis_fit <- log_reg_spec %>%
  fit(formula,
      data = train_abx)
underlying_sepsis <- extract_fit_engine(sepsis_fit)
coef(underlying_sepsis) %>% round(2)

###Evaluate model on test data
sepsis_test_probs <- predict(underlying_sepsis, test_abx,type="response")

###Attach to probability predictions dataframe
sepsis_util_key <- ur_util %>% select(micro_specimen_id,sepsis_ae)
df <- df %>% 
  left_join(sepsis_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

##Sync columns
curr_service_columns_1 <- grep("^curr_service_", names(train_abx), value = TRUE)
curr_service_columns_2 <- grep("^curr_service_", names(df), value = TRUE)
columns_to_add <- setdiff(curr_service_columns_1, curr_service_columns_2)
for (col in columns_to_add) {
  df[[col]] <- 0
}

##Utility score calculation

###CDI risk utility
cdi_util_probs <- predict(underlying_cdi, df,type="response")
cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
  select(Value) %>% unlist()

###Toxicity risk utility
tox_util_probs <- predict(underlying_tox, df,type="response")
tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
  select(Value) %>% unlist()

####Sepsis adverse outcome risk utility
sepsis_util_probs <- predict(underlying_sepsis, df,type="response")

###UTI-specific utility
uti_specifics <- c("Nitrofurantoin")
uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
  select(Value) %>% unlist()

###Access category utility
access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
access_abs <- c(access_abs, access_combos)

access_value <- scores[rownames(scores)=="Access",] %>% 
  select(Value) %>% unlist()

###Oral option utilituy
oral_abs <- c("AMP","SAM","CIP",
              "SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
oral_abs <- c(oral_abs, oral_combos)

oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
  select(Value) %>% unlist()

###IV option utility
iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-") 
iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
iv_abs <- c(iv_abs, iv_combos)

iv_value <- scores[rownames(scores)=="IV_option",] %>% 
  select(Value) %>% unlist()

###Reserve category utility
reserve_abs <- c()
reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
  select(Value) %>% unlist()

###High-cost agent utility
highcost_abs <- c()
cost_value <- scores[rownames(scores)=="High_cost",] %>% 
  select(Value) %>% unlist()

###AST R result utility
Rval_key <- ur_util %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)

###Enterococcus removal sensitivity analysis key
ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
  rename(ent_R = "R", ent_S = "S")

###Attach individual utilities to dataframe
df <- df %>% 
  mutate(prob_CDI = cdi_util_probs,
         value_CDI = cdi_value,
         util_CDI = prob_CDI * value_CDI,
         prob_tox = tox_util_probs,
         value_tox = tox_value,
         util_tox = prob_tox * value_tox,
         UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
         value_UTI = uti_value,
         util_uti = UTI_specific * value_UTI,
         Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
         value_access = access_value,
         util_access = Access_agent * value_access,
         Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
         value_oral = oral_value,
         util_oral = Oral_agent * value_oral,
         IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
         value_iv = iv_value,
         util_iv = IV_agent * value_iv,
         Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
         value_reserve = reserve_value,
         util_reserve = Reserve_agent * value_reserve,
         Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
         value_highcost = cost_value,
         util_highcost = Highcost_agent * value_highcost,
         prob_sepsisae = sepsis_util_probs,
         single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE)
  ) %>% left_join(Rval_key,by="micro_specimen_id") %>% left_join(ent_sens_key,by=c("micro_specimen_id","Antimicrobial")) %>% 
  mutate(ent_R = case_when(is.na(ent_R)~1,TRUE~ent_R),
         ent_S = case_when(is.na(ent_S)~0,TRUE~ent_S))

###Calculate overall utility score
df <- df %>% calculate_utilities()

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
df <- df %>% filter(Antimicrobial %in% abx_in_train)

###Dataframe for formulary agent sensitivity analysis
form_df <- df %>% calculate_utilities(
  formulary_list = c("Ceftriaxone","Ciprofloxacin"))
form_df <- form_df %>% filter(Antimicrobial %in% abx_in_train)

##Utility analysis

###Utility data visualisation
df %>% group_by(Antimicrobial) %>% 
  summarise(Median_util=median(Rx_utility)) %>% 
  arrange(desc(Median_util))
formulary_agents <- c()

df

}

##Read_in and final preprocessing

###Read in
abx <- read_csv("interim_abx.csv")
train_abx <- read_csv("train_abx.csv")
test_abx <- read_csv("test_abx.csv")
util_probs_df <- read_csv("probs_df_overall.csv")
ent_probs_df <- read_csv("ent_probs_df.csv")
ur_util <- read_csv("interim_ur_util.csv")
micro <- read_csv("micro_clean2.csv")
mic_ref <- micro %>% anti_join(ur_util,by="subject_id")
results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")

###Re-factorising outcome variables on abx dataframes after read-in
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()

##Survey results

###Engineer scores dataframe
rankings <- results %>% select(10:ncol(results)) %>% slice(2:nrow(results))
colnames(rankings) <- c("Abelfenide","Acetemran","Adenomadin","Adrevenac",
                        "Amrodine","Choriotroban","Cormide","Decloxone",
                        "Dexaset","Endoxolol","Olanzasys","Pansolid",
                        "Protestryl")
Characteristic <- c("AWaRe","CDI","Toxicity","UTI","Oral","IV","Cost")
Abelfenide <- c("Reserve","Low","Low","No","Yes","Yes","High")
Acetemran <- c("Watch","Low","Low","Yes","Yes","No","Low")
Adenomadin <- c("Access","High","Low","No","Yes","Yes","Low")
Adrevenac <- c("Watch","High","Low","No","No","Yes","High")
Amrodine <- c("Reserve","High","Low","No","No","Yes","High")
Choriotroban <- c("Access","Low","High","No","No","Yes","Low")
Cormide <- c("Watch","Low","Low","No","No","Yes","High")
Decloxone <- c("Watch","Low","High","No","No","Yes","High")
Dexaset <- c("Access","Low","Low","No","Yes","Yes","Low")
Endoxolol <- c("Access","Low","Low","Yes","Yes","No","Low")
Olanzasys <- c("Watch","High","Low","No","No","Yes","Low")
Pansolid <- c("Access","Low","Low","No","Yes","No","Low")
Protestryl <- c("Watch","High","Low","No","Yes","Yes","Low")
characteristics <- data.frame(Abelfenide,Acetemran,Adenomadin,
                              Adrevenac,Amrodine,Choriotroban,
                              Cormide,Decloxone,Dexaset,
                              Endoxolol,Olanzasys,Pansolid,Protestryl) %>% 
  t() %>% data.frame()

colnames(characteristics) <- Characteristic
characteristics <- characteristics %>% 
  mutate(Antibiotic = rownames(characteristics)) %>% relocate(
    Antibiotic,.before = AWaRe
  )
rownames(characteristics) <- NULL
scores <- rankings %>% pivot_longer(Abelfenide:Protestryl)
colnames(scores) <- c("Antibiotic","Rank")
scores <- scores %>% left_join(characteristics,by="Antibiotic")

scores <- scores %>% mutate(Access = case_when(AWaRe=="Access"~1,TRUE~0),
                            Watch = case_when(AWaRe=="Watch"~1,TRUE~0),
                            Reserve = case_when(AWaRe=="Reserve"~1,TRUE~0))

scores <- scores %>% rename(CDI_highrisk = "CDI",
                            Toxicity_highrisk = "Toxicity",
                            UTI_specific = "UTI",
                            Oral_option = "Oral",
                            IV_option = "IV",
                            High_cost = "Cost")
scores[scores=="High"] <- "1"
scores[scores=="Low"] <- "0"
scores[scores=="Yes"] <- "1"
scores[scores=="No"] <- "0"

scores

scores <- scores %>% select(-AWaRe) %>%
  mutate(across(c(Rank, CDI_highrisk, Toxicity_highrisk, UTI_specific,Access,
                  Watch,Reserve), as.numeric)) 

repeat_val <- nrow(scores)/13
seq_rep <- seq(1,nrow(scores)/13)
seq_rep1 <- c()

for (i in seq_rep) {
  
  seq_rep2 <- rep(seq_rep[i],13)
  seq_rep1 <- append(seq_rep1,seq_rep2)
  
  
}

scores <- scores %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                            id = seq_rep1) %>% 
  mutate(across(Oral_option:High_cost,as.numeric))

##Extracting characteristic importance weights

##Apply ridge regression
mlogit_data <- mlogit.data(scores, choice = "Rank", shape = "long", 
                           chid.var = "id", alt.var = "Antibiotic", 
                           ranked = TRUE)

formula_no_int <- Rank ~ CDI_highrisk + Toxicity_highrisk +  
  Oral_option + UTI_specific + IV_option + High_cost + Access + Reserve

X <- model.matrix(formula_no_int, data = mlogit_data)
y <- mlogit_data$Rank
fit <- cv.glmnet(X, y, family = "multinomial", alpha = 0)

lambda_min <- fit$lambda.min
lambda_1se <- fit$lambda.1se
cat("Lambda.min:", lambda_min, "\n")
cat("Lambda.1se:", lambda_1se, "\n")

tmp_coeffs <- coef(fit, s = "lambda.min")
scores <- tmp_coeffs$`TRUE` %>% as.matrix() %>% data.frame()
colnames(scores) <- "Value"
scores <- scores %>% mutate(Coefficient = rownames(scores),
                            OR = exp(Value),
                            OR_dif = OR-1)
scores <- scores %>% mutate(colour = case_when(
  Value > 0 ~ "B", Value < 0 ~ "A"
)) %>% slice(-1) %>% slice(-1)

scores <- scores %>% mutate(stan_OR = (OR-1)/max(abs(OR-1)))
scores$Coefficient <- c("High CDI risk","High toxicity risk","Oral option",
                        "UTI-specific","IV option","High cost","Access category",
                        "Reserve category")

scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())

###Visualise odds ratios for antimicrobial characteristics
ORplot <- ggplot(scores,aes(x=OR_dif,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Odds ratio for drug selection") +
  ggtitle("The effect of different antimicrobial drug properties\non clinician prescribing preference in UTI scenario")+
  scale_x_continuous(labels = function(x) x+1)+
  geom_vline(xintercept = 0)

ggsave(glue("ORplot.pdf"), plot = ORplot, device = "pdf", width = 10, height = 8,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Antimicrobial dummy variables in probability prediction dataframe
util_probs_df$abx_name_ <- as.factor(util_probs_df$Antimicrobial)
util_probs_df <- util_probs_df %>% mutate(
  abx_name_ = str_replace_all(abx_name_,"-",".")
)
dummy_vars <- model.matrix(~ abx_name_ - 1, data = util_probs_df)
util_probs_df <- cbind(util_probs_df, dummy_vars) %>% tibble() %>% 
  select(-abx_name_)

##CDI prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
abx_columns <- grep("^abx_name_", names(train_abx), value = TRUE)
service_columns <- grep("^curr_service_", names(train_abx), value = TRUE)
predictors <- c("pCDI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                abx_columns, service_columns, "pICU", "pSEPSIS")
predictors <- predictors[predictors != "abx_name_Ampicillin"]
formula <- reformulate(predictors, response = "CDI")
cdi_fit <- log_reg_spec %>%
  fit(formula,data = train_abx)

underlying_cdi <- extract_fit_engine(cdi_fit)
coef(underlying_cdi) %>% round(2)

###Evaluate on test data
cdi_test_probs <- predict(underlying_cdi, test_abx,type="response")
roc_curve <- roc(test_abx$CDI, cdi_test_probs)
plot(roc_curve, col = "darkgreen", lwd = 2,main="CDI prediction ROC")
auc_value <- pROC::auc(roc_curve)
print(auc_value)

###Attach to probability prediction dataframe
cdi_util_key <- ur_util %>% select(micro_specimen_id,pHADM,MALE,pICU,
                                   CDI:pSEPSIS) %>% 
  select(-AKI)
util_probs_df <- util_probs_df %>% 
  left_join(cdi_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

##Toxicity prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
predictors <- c("prAKI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                abx_columns, service_columns, "pICU", "pSEPSIS","admission_sepsis",
                "SIRS","highCRP")
predictors <- predictors[predictors != "abx_name_Ampicillin"]
formula <- reformulate(predictors, response = "overall_tox")
tox_fit <- log_reg_spec %>%
  fit(formula,
      data = train_abx)

underlying_tox <- extract_fit_engine(tox_fit)
coef(underlying_tox) %>% round(2)

###Evaluate model on test data
tox_test_probs <- predict(underlying_tox, test_abx,type="response")
roc_curve <- roc(test_abx$overall_tox, tox_test_probs)
plot(roc_curve, col = "blue", lwd = 2,main="Toxicity prediction ROC")
auc_value <- pROC::auc(roc_curve)
print(auc_value)

###Attach to probability prediction dataframe
tox_util_key <- ur_util %>% select(micro_specimen_id,overall_tox,
                                   admission_sepsis,
                                   SIRS,highCRP)
util_probs_df <- util_probs_df %>% 
  left_join(tox_util_key,by="micro_specimen_id")

##Sepsis adverse outcomes prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
predictors <- c("prAKI","pCDI", "pHADM", "age65", "pCKD", "pDIAB", "pLIVER", "pCARD", "pCVA", "pCA", "MALE",
                 service_columns, "pICU", "pSEPSIS")
formula <- reformulate(predictors, response = "sepsis_ae")
sepsis_fit <- log_reg_spec %>%
  fit(formula,
      data = train_abx)
underlying_sepsis <- extract_fit_engine(sepsis_fit)
coef(underlying_sepsis) %>% round(2)

###Evaluate model on test data
sepsis_test_probs <- predict(underlying_sepsis, test_abx,type="response")
roc_curve <- roc(test_abx$sepsis_ae, sepsis_test_probs)
plot(roc_curve, col = "red", lwd = 2,main="Sepsis adverse events prediction ROC")
auc_value <- pROC::auc(roc_curve)
print(auc_value)

###Attach to probability predictions dataframe
sepsis_util_key <- ur_util %>% select(micro_specimen_id,sepsis_ae)
util_probs_df <- util_probs_df %>% 
  left_join(sepsis_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

##Sync columns
curr_service_columns_1 <- grep("^curr_service_", names(train_abx), value = TRUE)
curr_service_columns_2 <- grep("^curr_service_", names(util_probs_df), value = TRUE)
columns_to_add <- setdiff(curr_service_columns_1, curr_service_columns_2)
for (col in columns_to_add) {
  util_probs_df[[col]] <- 0
}

##Utility score calculation

###CDI risk utility
cdi_util_probs <- predict(underlying_cdi, util_probs_df,type="response")
cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
  select(Value) %>% unlist()

###Toxicity risk utility
tox_util_probs <- predict(underlying_tox, util_probs_df,type="response")
tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
  select(Value) %>% unlist()

####Sepsis adverse outcome risk utility
sepsis_util_probs <- predict(underlying_sepsis, util_probs_df,type="response")

###UTI-specific utility
uti_specifics <- c("Nitrofurantoin")
uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
  select(Value) %>% unlist()

###Access category utility
access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
access_abs <- c(access_abs, access_combos)

access_value <- scores[rownames(scores)=="Access",] %>% 
  select(Value) %>% unlist()

###Oral option utilituy
oral_abs <- c("AMP","SAM","CIP",
              "SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
oral_abs <- c(oral_abs, oral_combos)

oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
  select(Value) %>% unlist()

###IV option utility
iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-") 
iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
iv_abs <- c(iv_abs, iv_combos)

iv_value <- scores[rownames(scores)=="IV_option",] %>% 
  select(Value) %>% unlist()

###Reserve category utility
reserve_abs <- c()
reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
  select(Value) %>% unlist()

###High-cost agent utility
highcost_abs <- c()
cost_value <- scores[rownames(scores)=="High_cost",] %>% 
  select(Value) %>% unlist()

###AST R result utility
Rval_key <- ur_util %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)

###Enterococcus removal sensitivity analysis key
ent_sens_key <- ent_probs_df %>% select(micro_specimen_id,Antimicrobial,S,R) %>% 
  rename(ent_R = "R", ent_S = "S")

###Attach individual utilities to dataframe
util_probs_df <- util_probs_df %>% 
  mutate(prob_CDI = cdi_util_probs,
         value_CDI = cdi_value,
         util_CDI = prob_CDI * value_CDI,
         prob_tox = tox_util_probs,
         value_tox = tox_value,
         util_tox = prob_tox * value_tox,
         UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
         value_UTI = uti_value,
         util_uti = UTI_specific * value_UTI,
         Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
         value_access = access_value,
         util_access = Access_agent * value_access,
         Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
         value_oral = oral_value,
         util_oral = Oral_agent * value_oral,
         IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
         value_iv = iv_value,
         util_iv = IV_agent * value_iv,
         Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
         value_reserve = reserve_value,
         util_reserve = Reserve_agent * value_reserve,
         Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 0.25, TRUE~0),
         value_highcost = cost_value,
         util_highcost = Highcost_agent * value_highcost,
         prob_sepsisae = sepsis_util_probs,
         single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE)
  ) %>% left_join(Rval_key,by="micro_specimen_id") %>% left_join(ent_sens_key,by=c("micro_specimen_id","Antimicrobial")) %>% 
  mutate(ent_R = case_when(is.na(ent_R)~1,TRUE~ent_R),
         ent_S = case_when(is.na(ent_S)~0,TRUE~ent_S))

write_csv(util_probs_df,"pre_util_probs_df.csv")

###Calculate overall utility score
util_probs_df <- util_probs_df %>% calculate_utilities()

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
util_probs_df <- util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

###Dataframe for formulary agent sensitivity analysis
form_util_probs_df <- util_probs_df %>% calculate_utilities(
  formulary_list = c("Ceftriaxone","Ciprofloxacin"))
form_util_probs_df <- form_util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

##Utility analysis

###Utility data visualisation
util_probs_df %>% group_by(Antimicrobial) %>% 
  summarise(Median_util=median(Rx_utility)) %>% 
  arrange(desc(Median_util))
formulary_agents <- c()
util_probs_df %>% utility_plot(Rx_utility,"Treatment")
util_probs_df %>% utility_plot(Rx_utility,"Treatment", " (single agent)")
util_probs_df %>% utility_plot(Rx_utility,"Treatment", " (combinations)")
util_probs_df %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent)")
util_probs_df %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations)")
util_probs_df %>% utility_plot(Outpatient_Rx_utility,"Oral treatment")
util_probs_df %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent)")
util_probs_df %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations)")
util_probs_df %>% utility_plot(AST_utility,"AST")

##Sensitivity analysis

###Probability density check across all antimicrobials
all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN")
long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_abs, all_combos)
long_allcombos <- combn(long_allabs, 2, FUN = function(x) paste(x, collapse = "_"))
long_allabs <- c(long_allabs, long_allcombos)

for (i in seq_along(long_allabs)) {
  
  util_probs_df %>% dens_check(long_allabs[i],S,1,4) 
  
}

###Sensitivity analysis varying resistance probabilities
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN") %>% ab_name() %>% 
  str_replace("/","-")
rx_dens_sens <- ur_util %>% dens_sens(util_probs_df,Rx_utility)
ast_dens_sens <- ur_util %>% dens_sens(util_probs_df,AST_utility)

rx_dens_sens %>% dens_sens_plot("Treatment")
ast_dens_sens %>% dens_sens_plot("Testing")

###Sensitivity analysis varying weighting values of different factors
uti_sens_df <- ur_util %>% uti_util_sens(util_probs_df,Rx_utility)
access_sens_df <- ur_util %>% access_util_sens(util_probs_df,Rx_utility)
oral_sens_df <- ur_util %>% oral_util_sens(util_probs_df,Rx_utility)
iv_sens_df <- ur_util %>% iv_util_sens(util_probs_df,Rx_utility)
cdi_sens_df <- ur_util %>% cdi_util_sens(util_probs_df,Rx_utility)
tox_sens_df <- ur_util %>% tox_util_sens(util_probs_df,Rx_utility)

uti_sens_df %>% util_sens_plot("Treatment","UTI-specificity")
access_sens_df %>% util_sens_plot("Treatment","Access category")
oral_sens_df %>% util_sens_plot("Treatment","Oral option")
iv_sens_df %>% util_sens_plot("Treatment","IV option")
cdi_sens_df %>% util_sens_plot("Treatment","CDI")
tox_sens_df %>% util_sens_plot("Treatment","Toxicity")

cdi_sens_df_ast <- ur_util %>% cdi_util_sens(util_probs_df,AST_utility)
uti_sens_df_ast <- ur_util %>% uti_util_sens(util_probs_df,AST_utility)
access_sens_df_ast <- ur_util %>% access_util_sens(util_probs_df,AST_utility)
oral_sens_df_ast <- ur_util %>% oral_util_sens(util_probs_df,AST_utility)
iv_sens_df_ast <- ur_util %>% iv_util_sens(util_probs_df,AST_utility)
cdi_sens_df_ast <- ur_util %>% cdi_util_sens(util_probs_df,AST_utility)
tox_sens_df_ast <- ur_util %>% tox_util_sens(util_probs_df,AST_utility)

uti_sens_df_ast %>% util_sens_plot("Testing","UTI-specificity")
access_sens_df_ast %>% util_sens_plot("Testing","Access category")
oral_sens_df_ast %>% util_sens_plot("Testing","Oral option")
iv_sens_df_ast %>% util_sens_plot("Testing","IV option")
cdi_sens_df_ast %>% util_sens_plot("Testing","CDI")
tox_sens_df_ast %>% util_sens_plot("Testing","Toxicity")

###Sensitivity analysis with variation of CDI, toxicity and sepsis AE risk
cdi_prob_df <- ur_util %>% cdi_prob_sens(util_probs_df,Rx_utility,"cdi_prob","prob_CDI",cdi_prob,prob_CDI,"CDI")
tox_prob_df <- ur_util %>% tox_prob_sens(util_probs_df,Rx_utility,"tox_prob","prob_tox",tox_prob,prob_tox,"toxicity")
sepsisae_prob_df <- ur_util %>% sepsisae_prob_sens(util_probs_df,Rx_utility,"sepsisae_prob","prob_sepsisae",sepsisae_prob,prob_sepsisae,"sepsis adverse events")

cdi_prob_df %>% dens_sens_plot_2("CDI","Treatment",cdi_prob)
tox_prob_df %>% dens_sens_plot_2("toxicity","Treatment",tox_prob)
sepsisae_prob_df %>% dens_sens_plot_2("sepsis adverse events","Treatment",sepsisae_prob)

cdi_prob_df_ast <- ur_util %>% cdi_prob_sens(util_probs_df,AST_utility,"cdi_prob","prob_CDI",cdi_prob,prob_CDI,"CDI")
tox_prob_df_ast <- ur_util %>% tox_prob_sens(util_probs_df,AST_utility,"tox_prob","prob_tox",tox_prob,prob_tox,"toxicity")
sepsisae_prob_df_ast <- ur_util %>% sepsisae_prob_sens(util_probs_df,AST_utility,"sepsisae_prob","prob_sepsisae",sepsisae_prob,prob_sepsisae,"sepsis adverse events")

cdi_prob_df_ast %>% dens_sens_plot_2("CDI","AST",cdi_prob)
tox_prob_df_ast %>% dens_sens_plot_2("toxicity","AST",tox_prob)
sepsisae_prob_df_ast %>% dens_sens_plot_2("sepsis adverse events","AST",sepsisae_prob)

###Sensitivity analysis with Enterococci not used in probability predictions
util_probs_df %>% utility_plot(ent_Rx_utility,"Treatment", " (Enterococcus removed)")
util_probs_df %>% utility_plot(ent_Rx_utility,"Treatment", " (Enterococcus removed, single agent)")
util_probs_df %>% utility_plot(ent_Rx_utility,"Treatment", " (Enterococcus removed, combinations)")
util_probs_df %>% utility_plot(ent_Urosepsis_Rx_utility,"Intravenous treatment", " (Enterococcus removed, single agent)")
util_probs_df %>% utility_plot(ent_Urosepsis_Rx_utility,"Intravenous treatment", " (Enterococcus removed, combinations)")
util_probs_df %>% utility_plot(ent_Outpatient_Rx_utility,"Oral treatment", " (Enterococcus removed)")
util_probs_df %>% utility_plot(ent_Outpatient_Rx_utility,"Oral treatment", " (Enterococcus removed, single agent)")
util_probs_df %>% utility_plot(ent_Outpatient_Rx_utility,"Oral treatment", " (Enterococcus removed, combinations)")
util_probs_df %>% utility_plot(ent_AST_utility,"AST", " (Enterococcus removed)")

###Write utility function dataframes
write_csv(util_probs_df,"utility_dataframe.csv")
write_csv(form_util_probs_df,"form_utility_dataframe.csv")

##Sensitivity analysis - separated specialties for questionnaire

###Read-in and filter questionnaire dataframes
results <- read_csv("labelled_ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
util_probs_df_2 <- read_csv("probs_df_overall.csv")
results_infection <- results %>% filter(Specialty=="Infection") %>% select(-Specialty)
results_itu <- results %>% filter(Specialty=="Intensive care") %>% select(-Specialty)
results_medicine <- results %>% filter(Specialty=="Medicine") %>% select(-Specialty)
results_surgery <- results %>% filter(Specialty=="Surgery") %>% select(-Specialty)
results_gp <- results %>% filter(Specialty=="General Practice") %>% select(-Specialty)

###Examine coefficients and apply to utility function
util_infection <- util_probs_df_2 %>% spec_sens_analysis(results_infection,"Infection")
util_itu <- util_probs_df_2 %>% spec_sens_analysis(results_itu,"Intensive care")
util_medicine <- util_probs_df_2 %>% spec_sens_analysis(results_medicine,"Medicine")
util_surgery <- util_probs_df_2 %>% spec_sens_analysis(results_surgery,"Surgery")
util_gp <- util_probs_df_2 %>% spec_sens_analysis(results_gp,"General Practice")

###Examine utility function output for each group
util_infection %>% utility_plot(Rx_utility,"Treatment",", infection specialties")
util_infection %>% utility_plot(Rx_utility,"Treatment", " (single agent, infection specialties)")
util_infection %>% utility_plot(Rx_utility,"Treatment", " (combinations, infection specialties)")
util_infection %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, infection specialties)")
util_infection %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, infection specialties)")
util_infection %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", infection specialties")
util_infection %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, infection specialties)")
util_infection %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, infection specialties)")
util_infection %>% utility_plot(AST_utility,"AST",", infection specialties")

util_itu %>% utility_plot(Rx_utility,"Treatment",", intensive care")
util_itu %>% utility_plot(Rx_utility,"Treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Rx_utility,"Treatment", " (combinations)")
util_itu %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, intensive care)")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", intensive care")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, intensive care)")
util_itu %>% utility_plot(AST_utility,"AST",", intensive care")

util_medicine %>% utility_plot(Rx_utility,"Treatment",", medical specialties")
util_medicine %>% utility_plot(Rx_utility,"Treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Rx_utility,"Treatment", " (combinations, medical specialties)")
util_medicine %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, medical specialties)")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", medical specialties")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, medical specialties)")
util_medicine %>% utility_plot(AST_utility,"AST",", medical specialties")

util_surgery %>% utility_plot(Rx_utility,"Treatment",", surgical specialties")
util_surgery %>% utility_plot(Rx_utility,"Treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Rx_utility,"Treatment", " (combinations, surgical specialties)")
util_surgery %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, surgical specialties)")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", surgical specialties")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, surgical specialties)")
util_surgery %>% utility_plot(AST_utility,"AST",", surgical specialties")

util_gp %>% utility_plot(Rx_utility,"Treatment",", general practice")
util_gp %>% utility_plot(Rx_utility,"Treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Rx_utility,"Treatment", " (combinations, general practice)")
util_gp %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, general practice)")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", general practice")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, general practice)")
util_gp %>% utility_plot(AST_utility,"AST",", general practice")




