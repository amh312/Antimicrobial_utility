#UTILITY CALCULATIONS

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
calculate_utilities <- function(df,formulary_list=c(),R_weight=1,MEWS=0) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox*tox_AUC + util_CDI*CDI_AUC + util_iv*MEWS*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util*AUC,
                      AMPR_utility = case_when(
                        Antimicrobial=="Ampicillin" ~
                          abs(R_weight)*AMP_R_value*R, TRUE~0
                      ),
                      SAMR_utility = case_when(
                        Antimicrobial=="Ampicillin-sulbactam" ~
                          abs(R_weight)*SAM_R_value*R, TRUE~0
                      ),
                      TZPR_utility = case_when(
                        Antimicrobial=="Piperacillin-tazobactam" ~
                          abs(R_weight)*TZP_R_value*R, TRUE~0
                      ),
                      CZOR_utility = case_when(
                        Antimicrobial=="Cefazolin" ~
                          abs(R_weight)*CZO_R_value*R, TRUE~0
                      ),
                      CROR_utility = case_when(
                        Antimicrobial=="Ceftriaxone" ~
                          abs(R_weight)*CRO_R_value*R, TRUE~0
                      ),
                      CAZR_utility = case_when(
                        Antimicrobial=="Ceftazidime" ~
                          abs(R_weight)*CAZ_R_value*R, TRUE~0
                      ),
                      FEPR_utility = case_when(
                        Antimicrobial=="Cefepime" ~
                          abs(R_weight)*FEP_R_value*R, TRUE~0
                      ),
                      MEMR_utility = case_when(
                        Antimicrobial=="Meropenem" ~
                          abs(R_weight)*MEM_R_value*R, TRUE~0
                      ),
                      CIPR_utility = case_when(
                        Antimicrobial=="Ciprofloxacin" ~
                          abs(R_weight)*CIP_R_value*R, TRUE~0
                      ),
                      GENR_utility = case_when(
                        Antimicrobial=="Gentamicin" ~
                          abs(R_weight)*GEN_R_value*R, TRUE~0
                      ),
                      SXTR_utility = case_when(
                        Antimicrobial=="Trimethoprim-sulfamethoxazole" ~
                          abs(R_weight)*SXT_R_value*R, TRUE~0
                      ),
                      NITR_utility = case_when(
                        Antimicrobial=="Nitrofurantoin" ~
                          abs(R_weight)*NIT_R_value*R, TRUE~0
                      ),
                      VANR_utility = case_when(
                        Antimicrobial=="Vancomycin" ~
                          abs(R_weight)*VAN_R_value*R, TRUE~0
                      ),
                      Formulary_agent = case_when(
                        Antimicrobial%in%formulary_list ~
                          TRUE, TRUE~FALSE
                      ),
                      Formulary_utility = 
                        abs(R_weight)*Formulary_agent*R,
                      AST_utility = (short_util+
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
                      Rx_utility = S_utility)
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~0,
        TRUE~Rx_utility),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~0,
        TRUE~Rx_utility)
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
  
  ggsave(glue("utility_{application}_{modification}.pdf"), plot = thisplot, device = "pdf", width = 10, height = 8,
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
    theme_minimal()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
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
    
    for(i in 0:10) {
      
      if (i==0) {
        
        sens_df <- probs_df_3 %>% mutate(
          R=case_when(Antimicrobial==iterabs[j]~0,TRUE~R),
          S=case_when(Antimicrobial==iterabs[j]~1,TRUE~S))
        
      } else if (i==10) {
        
        sens_df <- probs_df_3 %>% mutate(
          R=case_when(Antimicrobial==iterabs[j]~1,TRUE~R),
          S=case_when(Antimicrobial==iterabs[j]~0,TRUE~S)
        )
        
      } else {
      
      densy <- probs_df_2 %>% dens_check(iterabs[j],R,a,b)
      print(densy)
      
      sens_df <- probs_df_3 %>% dist_replace(df,iterabs[j],"R","S",a,b) }
      
      sens_df <- sens_df %>% calculate_utilities()
      
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
    xlim(0,1) +
    ggtitle(glue("Effect of varying {iterabs[i]} resistance\nprobability on {iterabs[i]} {measure} utility")) +
    xlab(glue("Population median probability of\n{iterabs[i]} resistance")) +
    ylab(glue("{measure} utility")) +
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(glue("dens_sens_1_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(df_plot)
  
  }
  
}

###Sensitivity analysis varying characteristic Desirability
R_util_sens <- function(df,probs_df,uf,min_val,max_val) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    probs_df_2 <- probs_df
      
      results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
      
      ##Survey results
      
      ###Engineer scores dataframe
      rankings <- results %>% select(10:ncol(results)) %>% dplyr::slice(2:nrow(results))
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
      
      for (p in seq_rep) {
        
        seq_rep2 <- rep(seq_rep[p],13)
        seq_rep1 <- append(seq_rep1,seq_rep2)
        
        
      }
      
      scores <- scores %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                                  id = seq_rep1) %>% 
        mutate(across(Oral_option:High_cost,as.numeric))
      
      ##Extracting characteristic Desirability weights
      
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
      )) %>% dplyr::slice(-1) %>% dplyr::slice(-1)
      
      scores <- scores %>% mutate(stan_OR = (OR-1)/max(abs(OR-1)))
      scores$Coefficient <- c("High CDI risk","High toxicity risk","Oral option",
                              "UTI-specific","IV option","High cost","Access category",
                              "Reserve category")
      
      scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())
      scores$Value <- scores$OR
      weight_sq <- seq(min_val,max_val,length.out=11)
      
      for(i in seq_along(weight_sq)) {
      
      scores2 <- scores  
      R_wt_value <- weight_sq[i]
      
      cdi_value <- scores2[rownames(scores2)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores2[rownames(scores2)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores2[rownames(scores2)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores2[rownames(scores2)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utilituy
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores2[rownames(scores2)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores2[rownames(scores2)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores2[rownames(scores2)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores2[rownames(scores2)=="High_cost",] %>% 
        select(Value) %>% unlist()
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_2 <- probs_df_2 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities(MEWS=R_wt_value,R_weight = 1)
      
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
      
      print(probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
              summarise(med_util=median(!!uf)))
      
      util_row <- probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% 
        summarise(med_util=median(!!uf),
                  best_ut=best_util,
                  lower_iqr=quantile(!!uf)[2],
                  upper_iqr=quantile(!!uf)[4],
                  spec_value=weight_sq[i],
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      util_cum <- data.frame(rbind(util_cum,util_row))
    
    }
    
  }
  
  util_cum
  
}
util_sens <- function(df,probs_df,uf,variable_criterion,min_val,max_val) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    probs_df_2 <- probs_df
    
    results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
    
    ##Survey results
    
    ###Engineer scores dataframe
    rankings <- results %>% select(10:ncol(results)) %>% dplyr::slice(2:nrow(results))
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
    
    for (p in seq_rep) {
      
      seq_rep2 <- rep(seq_rep[p],13)
      seq_rep1 <- append(seq_rep1,seq_rep2)
      
      
    }
    
    scores <- scores %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                                id = seq_rep1) %>% 
      mutate(across(Oral_option:High_cost,as.numeric))
    
    ##Extracting characteristic Desirability weights
    
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
    )) %>% dplyr::slice(-1) %>% dplyr::slice(-1)
    
    scores <- scores %>% mutate(stan_OR = (OR-1)/max(abs(OR-1)))
    scores$Coefficient <- c("High CDI risk","High toxicity risk","Oral option",
                            "UTI-specific","IV option","High cost","Access category",
                            "Reserve category")
    
    scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())
    scores$Value <- scores$OR
    weight_sq <- seq(min_val,max_val,length.out=11)
    
    for(i in seq_along(weight_sq)) {
      
      scores2 <- scores  
      R_wt_value <- 1
      scores2 <- scores2 %>% mutate(Value=case_when(rownames(scores2)==variable_criterion~weight_sq[i],
                                                    TRUE~Value))
      
      cdi_value <- scores2[rownames(scores2)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores2[rownames(scores2)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores2[rownames(scores2)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores2[rownames(scores2)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utility
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores2[rownames(scores2)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores2[rownames(scores2)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores2[rownames(scores2)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores2[rownames(scores2)=="High_cost",] %>% 
        select(Value) %>% unlist()
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_2 <- probs_df_2 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities(R_weight=R_wt_value)
      print(probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% summarise(medac=median(Rx_utility)))
      
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
                  spec_value=weight_sq[i],
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
    }
    
  }
  
  util_cum
  
}

###Data visualisation of Desirability sensitivity analysis
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
      ggtitle(glue("Effect of varying {value} weight on\n{iterabs[i]} {measure} utility")) +
      xlab(glue("{value}")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())+
      geom_vline(xintercept = 0, color = "lightgrey", size = 0.5)
    
    ggsave(glue("util_sens_{value}_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
}
R_util_sens_plot <- function(df,measure,value) {
  
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
      ggtitle(glue("Effect of varying {value} on\n{iterabs[i]} {measure} utility")) +
      xlab(glue("{value}")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_x_continuous(breaks = seq(min(df$spec_value), max(df$spec_value), by = 1))
    
    ggsave(glue("dens_sens_{value}_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
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
    
    for(i in 0:10) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      
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
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_3 <- probs_df_3 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      if (i==0){
        
        prob_vector <- probs_df_3 %>% mutate(
          !!char_col2:=case_when(Antimicrobial==iterabs[j]~0,
          TRUE~!!char_col2)) %>% 
          pull(!!char_col2) 
        
      } else if (i==10) {
        
        prob_vector <- probs_df_3 %>% mutate(
          !!char_col2:=case_when(Antimicrobial==iterabs[j]~1,
          TRUE~!!char_col2)) %>% 
          pull(!!char_col2) 
        
      } else {
      
      prob_vector <- probs_df_3 %>% dist_replace_2(df,iterabs[j],characteristic_col,a,b) %>% 
        pull(!!char_col2) }
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###calculate overall utility score
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
    
    for(i in 0:10) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      
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
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_3 <- probs_df_3 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      if (i==0){
        
        prob_vector <- probs_df_3 %>% mutate(
          !!char_col2:=case_when(Antimicrobial==iterabs[j]~0,
                                 TRUE~!!char_col2)) %>% 
          pull(!!char_col2) 
        
      } else if (i==10) {
        
        prob_vector <- probs_df_3 %>% mutate(
          !!char_col2:=case_when(Antimicrobial==iterabs[j]~1,
                                 TRUE~!!char_col2)) %>% 
          pull(!!char_col2) 
        
      } else {
        
        prob_vector <- probs_df_3 %>% dist_replace_2(df,iterabs[j],characteristic_col,a,b) %>% 
          pull(!!char_col2) }
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###calculate overall utility score
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
    
    for(i in 0:10) {
      
      densy <- probs_df_2 %>% dens_check_2(iterabs[j],!!char_col2,a,b,characteristic_name)
      print(densy)
      
      ##Utility score calculation
      
      ###CDI risk utility
      cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
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
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_3 <- probs_df_3 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      if (i==0){
        
        prob_vector <- rep(0,nrow(probs_df_3))
        
      } else if (i==10) {
        
        prob_vector <- rep(1,nrow(probs_df_3))
        
      } else {
        
        prob_vector <- rbeta(nrow(probs_df_3),a,b) }
      
      
      ###Attach individual utilities to dataframe
      probs_df_3 <- probs_df_3 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      ###calculate overall utility score
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

###Sensitivity analysis varying AUCs
util_sens_auc <- function(df,probs_df,uf,min_val,max_val) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
    probs_df_2 <- probs_df
    
    results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
    
    ##Survey results
    
    ###Engineer scores dataframe
    rankings <- results %>% select(10:ncol(results)) %>% dplyr::slice(2:nrow(results))
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
    
    for (p in seq_rep) {
      
      seq_rep2 <- rep(seq_rep[p],13)
      seq_rep1 <- append(seq_rep1,seq_rep2)
      
      
    }
    
    scores <- scores %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                                id = seq_rep1) %>% 
      mutate(across(Oral_option:High_cost,as.numeric))
    
    ##Extracting characteristic Desirability weights
    
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
    )) %>% dplyr::slice(-1) %>% dplyr::slice(-1)
    
    scores <- scores %>% mutate(stan_OR = (OR-1)/max(abs(OR-1)))
    scores$Coefficient <- c("High CDI risk","High toxicity risk","Oral option",
                            "UTI-specific","IV option","High cost","Access category",
                            "Reserve category")
    
    scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())
    scores$Value <- scores$OR
    weight_sq <- seq(min_val,max_val,length.out=11)
    
    for(i in seq_along(weight_sq)) {
      
      scores2 <- scores  
      R_wt_value <- 1
      
      cdi_value <- scores2[rownames(scores2)=="CDI_highrisk",] %>% 
        select(Value) %>% unlist()
      
      ###Toxicity risk utility
      tox_value <- scores2[rownames(scores2)=="Toxicity_highrisk",] %>% 
        select(Value) %>% unlist()
      
      
      ###UTI-specific utility
      uti_specifics <- c("Nitrofurantoin")
      uti_value <- scores2[rownames(scores2)=="UTI_specific",] %>% 
        select(Value) %>% unlist()
      
      ###Access category utility
      access_abs <- c("AMP","SAM","CZO",
                      "GEN","SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      access_abs <- c(access_abs, access_combos)
      
      access_value <- scores2[rownames(scores2)=="Access",] %>% 
        select(Value) %>% unlist()
      
      ###Oral option utility
      oral_abs <- c("AMP","SAM","CIP",
                    "SXT","NIT") %>% ab_name() %>% 
        str_replace("/","-")
      oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      oral_abs <- c(oral_abs, oral_combos)
      
      oral_value <- scores2[rownames(scores2)=="Oral_option",] %>% 
        select(Value) %>% unlist()
      
      ###IV option utility
      iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN") %>% ab_name() %>% 
        str_replace("/","-") 
      iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
      iv_abs <- c(iv_abs, iv_combos)
      
      iv_value <- scores2[rownames(scores2)=="IV_option",] %>% 
        select(Value) %>% unlist()
      
      ###Reserve category utility
      reserve_abs <- c()
      reserve_value <- scores2[rownames(scores2)=="Reserve",] %>% 
        select(Value) %>% unlist()
      
      ###High-cost agent utility
      highcost_abs <- c()
      cost_value <- scores2[rownames(scores2)=="High_cost",] %>% 
        select(Value) %>% unlist()
      cost_list <- read_csv("us_drug_cost_list.csv")
      drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
      cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
        mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
      comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
      cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
      cost_list2 <- data.frame(
        Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
        min_cost = cost_comb
      )  
      cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
      cost_list <- data.frame(rbind(cost_list,cost_list2))
      cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                        Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
      probs_df_2 <- probs_df_2 %>% left_join(cost_list)
      
      ###AST R result utility
      Rval_key <- df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
      
      ###Attach individual utilities to dataframe
      probs_df_2 <- probs_df_2 %>% 
        mutate(value_CDI = cdi_value,
               util_CDI = abs_calc(value_CDI,prob_CDI),
               value_tox = tox_value,
               util_tox = abs_calc(value_tox,prob_tox),
               UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
               value_UTI = uti_value,
               util_uti = abs_calc(value_UTI,UTI_specific),
               Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
               value_access = access_value,
               util_access = abs_calc(value_access,Access_agent),
               Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
               value_oral = oral_value,
               util_oral = abs_calc(value_oral,Oral_agent),
               IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
               value_iv = iv_value,
               util_iv = abs_calc(value_iv,IV_agent),
               Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
               value_reserve = reserve_value,
               util_reserve = abs_calc(value_reserve,Reserve_agent),
               Highcost_agent = min_cost,
               value_highcost = cost_value,
               util_highcost = abs_calc(value_highcost,Highcost_agent),
               single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE))
      
      probs_df_2 <- probs_df_2 %>% mutate(AUC=case_when(Antimicrobial==iterabs[j]~weight_sq[i],
                                                        TRUE~AUC))
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities(R_weight=R_wt_value)
      print(probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% summarise(medac=median(Rx_utility)))
      
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
                  spec_value=weight_sq[i],
                  overall_med=overall_median,
                  Antimicrobial=iterabs[j])
      
      
      util_cum <- data.frame(rbind(util_cum,util_row))
      
    }
    
  }
  
  util_cum
  
}

###Data visualisation of analysis varying AUCs
auc_sens_plot <- function(df,measure,uf) {
  
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
    
    df_plot <- ggplot(df_spec_plot, aes(x = spec_value)) +
      geom_line(aes(y = as.numeric(med_util), group = Antimicrobial, color = Antimicrobial)) +
      geom_ribbon(aes(y = as.numeric(med_util),
                      ymin = as.numeric(lower_iqr),
                      ymax = as.numeric(upper_iqr),
                      group = Antimicrobial, fill = Antimicrobial), alpha = 0.3) +
      xlim(0,1) +
      ggtitle(glue("Effect of varying {iterabs[i]} prediction\nAUROC on {iterabs[i]} {measure} utility")) +
      xlab(glue("Population median probability of\n{iterabs[i]} resistance")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    ggsave(glue("auc_sensplot_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
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
      ggtitle(glue("Effect of {characteristic} on\n{iterabs[i]} {measure} utility")) +
      xlab(glue("{characteristic}")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    ggsave(glue("dens_sens_2_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
}

###Separated-specialty sensitivity analysis
spec_sens_analysis <- function(df,results_df,specialty) {
###Engineer scores dataframe
rankings <- results_df %>% select(10:ncol(results_df)) %>% dplyr::slice(2:nrow(results_df))
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

###Apply ranked logit regression
mlogit_data <- mlogit.data(scores, 
                           choice = "Rank", 
                           shape = "long", 
                           chid.var = "id", 
                           alt.var = "Antibiotic", 
                           ranked = TRUE)

formula_no_int <- Rank ~ CDI_highrisk + Toxicity_highrisk +  
  Oral_option + UTI_specific + IV_option + High_cost + Access + Reserve | 0
rol_model <- mlogit(formula_no_int, data = mlogit_data, method = "bfgs")
coefficients <- coef(rol_model)
scores <- data.frame(
  Coefficient = names(coefficients),
  Value = coefficients,
  OR = exp(coefficients),
  OR_dif = exp(coefficients) - 1
)

###Add coefficients to scores dataframe
scores <- scores %>% 
  mutate(stan_OR = (OR - 1) / max(abs(OR - 1))) %>% 
  mutate(colour = case_when(
    Value > 0 ~ "B", 
    Value < 0 ~ "A"
  )) %>% 
  mutate(Coefficient = c("High CDI risk","High toxicity risk","Oral option",
                         "UTI-specific","IV option","High cost","Access category",
                         "Reserve category")) %>% 
  arrange(Value)

###Export unstandardised weights for reference
unstan_vals <- scores %>% select(Value) %>% 
  mutate(Name = c("Wc", "Wt", "Wo", "Wu", "Wi", "Wh", "Wa", "Wr"))
write_csv(unstan_vals, "unstand_weights.csv")

###Set default R weight value
R_wt_value <- 0

###Visualise odds ratios for antimicrobial characteristics
scores$Coefficient <- factor(scores$Coefficient, levels=
                               scores %>% arrange(Value) %>% 
                               select(Coefficient) %>% unlist())

sens_features <- ggplot(scores,aes(x=Value,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle(glue("The effect of different antimicrobial drug properties\non {specialty} clinician prescribing preference in UTI scenario"))+
  geom_vline(xintercept = 0)

ggsave(glue("{specialty}_importances.pdf"), plot = sens_features, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")
print(sens_features)





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

###CDI risk utility
cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
  select(Value) %>% unlist()

###Toxicity risk utility
tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
  select(Value) %>% unlist()

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
cost_list <- read_csv("us_drug_cost_list.csv")
drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(min_cost=min(Cost)) %>% ungroup() %>% 
  mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
cost_comb <- combn(cost_list$min_cost, 2, function(x) sum(x))
cost_list2 <- data.frame(
  Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
  min_cost = cost_comb
)  
cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
cost_list <- data.frame(rbind(cost_list,cost_list2))
cost_list <- cost_list %>% mutate(min_cost = min_cost/max(min_cost),
                                  Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
df <- df %>% left_join(cost_list)

###AST R result utility
Rval_key <- ur_util %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)

###Attach individual utilities to dataframe
df <- df %>% 
  mutate(value_CDI = cdi_value,
         util_CDI = abs_calc(value_CDI,prob_CDI),
         value_tox = tox_value,
         util_tox = abs_calc(value_tox,prob_tox),
         UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
         value_UTI = uti_value,
         util_uti = abs_calc(value_UTI,UTI_specific),
         Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
         value_access = access_value,
         util_access = abs_calc(value_access,Access_agent),
         Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
         value_oral = oral_value,
         util_oral = abs_calc(value_oral,Oral_agent),
         IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
         value_iv = iv_value,
         util_iv = abs_calc(value_iv,IV_agent),
         Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
         value_reserve = reserve_value,
         util_reserve = abs_calc(value_reserve,Reserve_agent),
         Highcost_agent = min_cost,
         value_highcost = cost_value,
         util_highcost = abs_calc(value_highcost,Highcost_agent),
         single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE)
  ) %>% left_join(Rval_key,by="micro_specimen_id")

df <- df %>% left_join(auc_df,by="Antimicrobial")
df$CDI_AUC <- auc_df %>% filter(Antimicrobial=="CDI") %>%
  select(AUC) %>% unlist
df$tox_AUC <- auc_df %>% filter(Antimicrobial=="toxicity") %>%
  select(AUC) %>% unlist
df <- df %>% left_join(recall_df,by="Antimicrobial")
df$CDI_recall <- recall_df %>% filter(Antimicrobial=="CDI") %>%
  select(recall) %>% unlist
df$tox_recall <- recall_df %>% filter(Antimicrobial=="toxicity") %>%
  select(recall) %>% unlist
df <- df %>% left_join(precision_df,by="Antimicrobial")
df$CDI_precision <- precision_df %>% filter(Antimicrobial=="CDI") %>%
  select(precision) %>% unlist
df$tox_precision <- precision_df %>% filter(Antimicrobial=="toxicity") %>%
  select(precision) %>% unlist
df <- df %>% left_join(accuracy_df,by="Antimicrobial")
df$CDI_accuracy <- accuracy_df %>% filter(Antimicrobial=="CDI") %>%
  select(accuracy) %>% unlist
df$tox_accuracy <- accuracy_df %>% filter(Antimicrobial=="toxicity") %>%
  select(accuracy) %>% unlist
df <- df %>% left_join(f1_score_df,by="Antimicrobial")
df$CDI_f1_score <- f1_score_df %>% filter(Antimicrobial=="CDI") %>%
  select(f1_score) %>% unlist
df$tox_f1_score <- f1_score_df %>% filter(Antimicrobial=="toxicity") %>%
  select(f1_score) %>% unlist

###Calculate overall utility score
df <- df %>% calculate_utilities(R_weight = R_wt_value)

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
ur_util <- read_csv("interim_ur_util.csv")
micro <- read_csv("micro_clean2.csv")
mic_ref <- micro %>% anti_join(ur_util,by="subject_id")
results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")

###Make AUC list
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
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
metrics_ablist <- c(combined_antimicrobial_map,"CDI","toxicity")
auc_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
auc_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  auc_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("AUC"))
  
}
colnames(auc_df) <- c("Antimicrobial","AUC")
auc_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
auc_df[38,1] <- "CDI"
auc_df[39,1] <- "toxicity"
auc_df$Antimicrobial <- as.character(auc_df$Antimicrobial)
write_csv(auc_df,"auc_df.csv")
recall_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
recall_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  recall_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Recall"))
  
}
colnames(recall_df) <- c("Antimicrobial","recall")
recall_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
recall_df[38,1] <- "CDI"
recall_df[39,1] <- "toxicity"
recall_df$Antimicrobial <- as.character(recall_df$Antimicrobial)
write_csv(recall_df,"recall_df.csv")
precision_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
precision_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  precision_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Precision"))
  
}
colnames(precision_df) <- c("Antimicrobial","precision")
precision_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
precision_df[38,1] <- "CDI"
precision_df[39,1] <- "toxicity"
precision_df$Antimicrobial <- as.character(precision_df$Antimicrobial)
write_csv(precision_df,"precision_df.csv")
accuracy_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
accuracy_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  accuracy_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("Accuracy"))
  
}
colnames(accuracy_df) <- c("Antimicrobial","accuracy")
accuracy_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
accuracy_df[38,1] <- "CDI"
accuracy_df[39,1] <- "toxicity"
accuracy_df$Antimicrobial <- as.character(accuracy_df$Antimicrobial)
write_csv(accuracy_df,"accuracy_df.csv")
f1_score_df <- data.frame(matrix(nrow=length(metrics_ablist),ncol=0))
f1_score_df$Antimicrobial <- metrics_ablist
for (i in 1:length(metrics_ablist)) {
  
  f1_score_df[i,2] <- read_csv(glue("metrics_{metrics_ablist[i]}.csv")) %>% select(contains("F1"))
  
}
colnames(f1_score_df) <- c("Antimicrobial","f1_score")
f1_score_df$Antimicrobial[1:37] <- names(metrics_ablist)[1:37]
f1_score_df[38,1] <- "CDI"
f1_score_df[39,1] <- "toxicity"
f1_score_df$Antimicrobial <- as.character(f1_score_df$Antimicrobial)
write_csv(f1_score_df,"f1_score_df.csv")

###Re-factorising outcome variables on abx dataframes after read-in
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()

##Survey results

###Engineer scores dataframe
rankings <- results %>% select(10:ncol(results)) %>% dplyr::slice(2:nrow(results))
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

##Extract characteristic importance weights

###Apply ranked logit regression
mlogit_data <- mlogit.data(scores, 
                           choice = "Rank", 
                           shape = "long", 
                           chid.var = "id", 
                           alt.var = "Antibiotic", 
                           ranked = TRUE)

formula_no_int <- Rank ~ CDI_highrisk + Toxicity_highrisk +  
  Oral_option + UTI_specific + IV_option + High_cost + Access + Reserve | 0
rol_model <- mlogit(formula_no_int, data = mlogit_data, method = "bfgs")
coefficients <- coef(rol_model)
scores <- data.frame(
  Coefficient = names(coefficients),
  Value = coefficients,
  OR = exp(coefficients),
  OR_dif = exp(coefficients) - 1
)

###Add coefficients to scores dataframe
scores <- scores %>% 
  mutate(stan_OR = (OR - 1) / max(abs(OR - 1))) %>% 
  mutate(colour = case_when(
    Value > 0 ~ "B", 
    Value < 0 ~ "A"
  )) %>% 
  mutate(Coefficient = c("High CDI risk","High toxicity risk","Oral option",
                         "UTI-specific","IV option","High cost","Access category",
                         "Reserve category")) %>% 
  arrange(Value)

###Export unstandardised weights for reference
unstan_vals <- scores %>% select(Value) %>% 
  mutate(Name = c("Wc", "Wt", "Wo", "Wu", "Wi", "Wh", "Wa", "Wr"))
write_csv(unstan_vals, "unstand_weights.csv")

###Set default R weight value
R_wt_value <- 0

###Visualise odds ratios for antimicrobial characteristics
scores$Coefficient <- factor(scores$Coefficient, levels=
                         scores %>% arrange(Value) %>% 
                         select(Coefficient) %>% unlist())

ORplot <- ggplot(scores,aes(x=Value,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle("The effect of different antimicrobial drug properties\non clinician prescribing preference in UTI scenario")+
  geom_vline(xintercept = 0)

ggsave(glue("ORplot.pdf"), plot = ORplot, device = "pdf", width = 10, height = 8,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")
print(ORplot)

###Antimicrobial dummy variables in probability prediction dataframe
util_probs_df$abx_name_ <- as.factor(util_probs_df$Antimicrobial)
util_probs_df <- util_probs_df %>% mutate(
  abx_name_ = str_replace_all(abx_name_,"-",".")
)
dummy_vars <- model.matrix(~ abx_name_ - 1, data = util_probs_df)
util_probs_df <- cbind(util_probs_df, dummy_vars) %>% tibble() %>% 
  select(-abx_name_)

##Utility score calculation

###CDI risk utility
cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
  select(Value) %>% unlist()

###Toxicity risk utility
tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
  select(Value) %>% unlist()

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
cost_list <- read_csv("us_drug_cost_list.csv")
drug_order <- cost_list %>% distinct(`Generic Name`) %>% unlist()
cost_list <- cost_list %>% group_by(`Generic Name`) %>% summarise(orig_cost=min(Cost)) %>% ungroup() %>% 
  mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
comb <- combn(cost_list$`Generic Name`, 2, simplify = FALSE)
cost_comb <- combn(cost_list$orig_cost, 2, function(x) sum(x))
cost_list2 <- data.frame(
  Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
  orig_cost = cost_comb
)  
cost_list <- cost_list %>% rename(Antimicrobial="Generic Name")
cost_list <- data.frame(rbind(cost_list,cost_list2))
cost_list <- cost_list %>% mutate(min_cost = orig_cost/max(orig_cost),
         Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
util_probs_df <- util_probs_df %>% left_join(cost_list)

###AST R result utility
Rval_key <- ur_util %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)

###Attach individual utilities to dataframe
util_probs_df <- util_probs_df %>% 
  mutate(value_CDI = cdi_value,
         util_CDI = abs_calc(value_CDI,prob_CDI),
         value_tox = tox_value,
         util_tox = abs_calc(value_tox,prob_tox),
         UTI_specific = case_when(Antimicrobial %in% uti_specifics ~ 1, TRUE~0),
         value_UTI = uti_value,
         util_uti = abs_calc(value_UTI,UTI_specific),
         Access_agent = case_when(Antimicrobial %in% access_abs ~ 1, TRUE~0),
         value_access = access_value,
         util_access = abs_calc(value_access,Access_agent),
         Oral_agent = case_when(Antimicrobial %in% oral_abs ~ 1, TRUE~0),
         value_oral = oral_value,
         util_oral = abs_calc(value_oral,Oral_agent),
         IV_agent = case_when(Antimicrobial %in% iv_abs ~ 1, TRUE~0),
         value_iv = iv_value,
         util_iv = abs_calc(value_iv,IV_agent),
         Reserve_agent = case_when(Antimicrobial %in% reserve_abs ~ 1, TRUE~0),
         value_reserve = reserve_value,
         util_reserve = abs_calc(value_reserve,Reserve_agent),
         Highcost_agent = min_cost,
         value_highcost = cost_value,
         util_highcost = abs_calc(value_highcost,Highcost_agent),
         single_agent = case_when(!grepl("_",Antimicrobial) ~ TRUE, TRUE~FALSE)
  ) %>% left_join(Rval_key,by="micro_specimen_id")

util_probs_df <- util_probs_df %>% left_join(auc_df,by="Antimicrobial")
util_probs_df$CDI_AUC <- auc_df %>% filter(Antimicrobial=="CDI") %>%
  select(AUC) %>% unlist
util_probs_df$tox_AUC <- auc_df %>% filter(Antimicrobial=="toxicity") %>%
  select(AUC) %>% unlist
util_probs_df <- util_probs_df %>% left_join(recall_df,by="Antimicrobial")
util_probs_df$CDI_recall <- recall_df %>% filter(Antimicrobial=="CDI") %>%
  select(recall) %>% unlist
util_probs_df$tox_recall <- recall_df %>% filter(Antimicrobial=="toxicity") %>%
  select(recall) %>% unlist
util_probs_df <- util_probs_df %>% left_join(precision_df,by="Antimicrobial")
util_probs_df$CDI_precision <- precision_df %>% filter(Antimicrobial=="CDI") %>%
  select(precision) %>% unlist
util_probs_df$tox_precision <- precision_df %>% filter(Antimicrobial=="toxicity") %>%
  select(precision) %>% unlist
util_probs_df <- util_probs_df %>% left_join(accuracy_df,by="Antimicrobial")
util_probs_df$CDI_accuracy <- accuracy_df %>% filter(Antimicrobial=="CDI") %>%
  select(accuracy) %>% unlist
util_probs_df$tox_accuracy <- accuracy_df %>% filter(Antimicrobial=="toxicity") %>%
  select(accuracy) %>% unlist
util_probs_df <- util_probs_df %>% left_join(f1_score_df,by="Antimicrobial")
util_probs_df$CDI_f1_score <- f1_score_df %>% filter(Antimicrobial=="CDI") %>%
  select(f1_score) %>% unlist
util_probs_df$tox_f1_score <- f1_score_df %>% filter(Antimicrobial=="toxicity") %>%
  select(f1_score) %>% unlist

write_csv(util_probs_df,"pre_util_probs_df.csv")

###Calculate overall utility score
util_probs_df <- util_probs_df %>% calculate_utilities(R_weight = R_wt_value)

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

##Sensitivity analysis

###Stem reference lists
all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN")
long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_abs, all_combos)
long_allcombos <- combn(long_allabs, 2, FUN = function(x) paste(x, collapse = "_"))
long_allabs <- c(long_allabs, long_allcombos)

###Sensitivity analysis varying resistance probabilities
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN") %>% ab_name() %>% 
  str_replace("/","-")
rx_dens_sens <- ur_util %>% dens_sens(util_probs_df,Rx_utility)
rx_dens_sens %>% dens_sens_plot("Treatment")

###Sensitivity analysis varying weighting values of different factors
uti_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"UTI_specific",-2,2)
access_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"Access",-2,2)
oral_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"Oral_option",-2,2)
iv_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"IV_option",-2,2)
cdi_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"CDI_highrisk",-2,2)
tox_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"Toxicity_highrisk",-2,2)
highcost_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"High_cost",-2,2)
reserve_sens_df <- ur_util %>% util_sens(util_probs_df,Rx_utility,"Reserve",-2,2)
R_sens_df <- ur_util %>% R_util_sens(util_probs_df,Rx_utility,0,20)

uti_sens_df %>% util_sens_plot("Treatment","UTI-specificity")
access_sens_df %>% util_sens_plot("Treatment","Access category")
oral_sens_df %>% util_sens_plot("Treatment","Oral option")
iv_sens_df %>% util_sens_plot("Treatment","IV option")
cdi_sens_df %>% util_sens_plot("Treatment","CDI")
tox_sens_df %>% util_sens_plot("Treatment","Toxicity")
highcost_sens_df %>% util_sens_plot("Treatment","Cost")
reserve_sens_df %>% util_sens_plot("Treatment","Reserve category")
R_sens_df %>% R_util_sens_plot("Treatment","Illness severity")

###Sensitivity analysis with variation of CDI and toxicity risk
cdi_prob_df <- ur_util %>% cdi_prob_sens(util_probs_df,Rx_utility,"cdi_prob","prob_CDI",cdi_prob,prob_CDI,"CDI")
tox_prob_df <- ur_util %>% tox_prob_sens(util_probs_df,Rx_utility,"tox_prob","prob_tox",tox_prob,prob_tox,"toxicity")

cdi_prob_df %>% dens_sens_plot_2("CDI probability","Treatment",cdi_prob)
tox_prob_df %>% dens_sens_plot_2("toxicity probability","Treatment",tox_prob)

###Sensitivity analysis with variation of AUC
auc_prob_df <- util_sens_auc(ur_util,util_probs_df,Rx_utility,0,1)
auc_prob_df %>% auc_sens_plot("Treatment")

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

util_itu %>% utility_plot(Rx_utility,"Treatment",", intensive care")
util_itu %>% utility_plot(Rx_utility,"Treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Rx_utility,"Treatment", " (combinations)")
util_itu %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, intensive care)")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", intensive care")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, intensive care)")
util_itu %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, intensive care)")

util_medicine %>% utility_plot(Rx_utility,"Treatment",", medical specialties")
util_medicine %>% utility_plot(Rx_utility,"Treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Rx_utility,"Treatment", " (combinations, medical specialties)")
util_medicine %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, medical specialties)")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", medical specialties")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, medical specialties)")
util_medicine %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, medical specialties)")

util_surgery %>% utility_plot(Rx_utility,"Treatment",", surgical specialties")
util_surgery %>% utility_plot(Rx_utility,"Treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Rx_utility,"Treatment", " (combinations, surgical specialties)")
util_surgery %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, surgical specialties)")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", surgical specialties")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, surgical specialties)")
util_surgery %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, surgical specialties)")

util_gp %>% utility_plot(Rx_utility,"Treatment",", general practice")
util_gp %>% utility_plot(Rx_utility,"Treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Rx_utility,"Treatment", " (combinations, general practice)")
util_gp %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment", " (combinations, general practice)")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment",", general practice")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (single agent, general practice)")
util_gp %>% utility_plot(Outpatient_Rx_utility,"Oral treatment", " (combinations, general practice)")




