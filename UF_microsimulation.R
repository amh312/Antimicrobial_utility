#MAIN ANALYSIS

set.seed(123)

##Functions

###Utility calculation
calculate_utilities_illsens <- function(df,formulary_list=c(),R_weight=1) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox + util_CDI + util_iv*exp(acuity)*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util,
                      imp_S_utility=imp_S*overall_util,
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
                      Rx_utility = S_utility,
                      imp_Rx_utility = imp_S_utility)
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~0,
        TRUE~Rx_utility),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~0,
        TRUE~Rx_utility),
      imp_Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~0,
        TRUE~imp_Rx_utility),
      imp_Outpatient_Rx_utility = case_when(
        util_oral ==0 ~0,
        TRUE~imp_Rx_utility)
    )
  
}
calculate_utilities <- function(df,formulary_list=c(),R_weight=1) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox + util_CDI + util_iv*exp(acuity)*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util,
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
abs_calc <- function(val,prob) {
  
  ifelse(val>=0,val*prob,abs(val)*(1-prob))
  
}

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                overall_tox = factor(overall_tox),
                sepsis_ae=factor(sepsis_ae))
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
  
  sens_plots <- list()
  
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
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("dens_sens_1_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave("resistance_sens_plots_grid.pdf", plot = grid_plot, width = 20, height = 30)
  
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
    score_rownames <- rownames(scores)
    
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
      probs_df_2 <- probs_df_2 %>% calculate_utilities(R_weight = 1)
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
    score_rownames <- rownames(scores)
    
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
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
util_sens_illness <- function(df,probs_df,uf,min_val,max_val) {
  
  uf <- enquo(uf)
  
  util_cum <- data.frame(matrix(nrow = 0,ncol=7))
  colnames(util_cum) <- c("med_util","lower_iqr","upper_iqr","spec_value","overall_med","Antimicrobial")
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% iterabs)
  
  for (j in seq_along(iterabs)) {
    
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
    score_rownames <- rownames(scores)
    
    scores$Coefficient <- factor(scores$Coefficient, levels= scores %>% arrange(Value) %>% select(Coefficient) %>% unlist())
    scores$Value <- scores$OR
    weight_sq <- seq(min_val,max_val)
    
    for(i in seq_along(weight_sq)) {
      
      probs_df_2 <- probs_df %>% mutate(acuity=weight_sq[i])
      
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
      
      ###Calculate overall utility score
      probs_df_2 <- probs_df_2 %>% calculate_utilities(R_weight=R_wt_value)
      print(probs_df_2 %>% filter(Antimicrobial==iterabs[j]) %>% summarise(medac=median(Rx_utility)))
      
      ###Filter out combination predictions not present in training dataset
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
  sens_plots <- list()
  
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
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("util_sens_{value}_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave(glue("{value}_sens_plots_grid.pdf"), plot = grid_plot, width = 20, height = 30)
  
}
R_util_sens_plot <- function(df,measure,value) {
  
  iterabs <- all_singles %>% ab_name() %>% str_replace("/","-")
  sens_plots <- list()
  
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
      scale_x_continuous(breaks = seq(min(df$spec_value), max(df$spec_value), by = 5))
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("dens_sens_{value}_{iterabs[i]}_Rsens.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave("sens_plots_grid.pdf", plot = grid_plot, width = 20, height = 30)
  
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
               prob_CDI=prob_vector,
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
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
               prob_tox=prob_vector,
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
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
      abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
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
  sens_plots <- list()
  
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
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("dens_sens_2_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave(glue("{characteristic}_sens_plots_grid.pdf"), plot = grid_plot, width = 20, height = 30)
  
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
  scores_rownames <- rownames(scores)
  
  # Approximate confidence intervals using bootstrap approach
  boot_coefs <- matrix(NA, nrow = 1000, ncol = length(coef(mlogit(formula_no_int, data = mlogit_data, method = "bfgs"))))
  for (i in 1:1000) {
    set.seed(i)
    boot_index <- sample(unique(mlogit_data$id), replace = TRUE)
    boot_mlogit <- mlogit_data %>% subset(id %in% boot_index)
    boot_model <- mlogit(formula_no_int, data = boot_mlogit, method = "bfgs")
    boot_coefs[i, ] <- coef(boot_model)
  }
  iqs <- data.frame(Coefficient=factor(names(coef(boot_model))),Q5=NA,Q95=NA)
  iqs <- iqs %>% mutate(Coefficient=case_when(
    grepl("CDI",Coefficient)~ "High CDI risk",
    grepl("Toxicity",Coefficient)~ "High toxicity risk",
    grepl("UTI",Coefficient)~ "UTI-specific",
    grepl("Access",Coefficient)~ "Access category",
    grepl("Reserve",Coefficient)~ "Reserve category",
    TRUE~Coefficient
  )) %>% 
    mutate(Coefficient=str_replace_all(Coefficient,"_"," "))
  for (i in 1:ncol(boot_coefs)) {
    
    iqs$Q5[i] <- quantile(boot_coefs[,i],probs=seq(0,1,0.05))['5%']
    iqs$Q95[i] <- quantile(boot_coefs[,i],probs=seq(0,1,0.05))['95%']
    
  }
  scores <- scores %>% left_join(iqs)
  rownames(scores) <- scores_rownames
  
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
    geom_hline(aes(yintercept=0)) +
    ylab("Drug property") +
    xlab("Coefficient value for drug selection probability") +
    ggtitle(glue("The effect of different antimicrobial drug properties on\n{specialty} clinician prescribing preference in the UTI\nscenario discrete choice experiment"))+
    geom_vline(xintercept = 0,colour="grey40")+
    geom_errorbar(aes(xmin = Q5, xmax = Q95), width = 0.1,colour="grey40")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          legend.position = "None")
  
  ggsave(glue("{specialty}_importances.pdf"), plot = sens_features, device = "pdf", width = 8, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  print(sens_features)
  
  
  ###Antimicrobial dummy variables in probability prediction dataframe
  df$ab_name_ <- as.factor(df$Antimicrobial)
  df <- df %>% mutate(
    ab_name_ = str_replace_all(ab_name_,"-",".")
  )
  dummy_vars <- model.matrix(~ ab_name_ - 1, data = df)
  df <- cbind(df, dummy_vars) %>% tibble() %>% 
    select(-ab_name_)
  
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

###Susceptibility plot
susc_plotter_overall <- function(df,df2,subset="",measure,suffix="",agent_col1,agent_name1,agent_col2,agent_name2,variable) {
  
  agent_col1 <- enquo(agent_col1)
  agent_col2 <- enquo(agent_col2)
  
  df$Metric <- factor(df$Metric, levels=c("All agents","Access agents","IV agents","Oral agents"))
  
  if (grepl("severity", variable)) {
    
    if (grepl("Access", suffix)&!grepl("inpatient", suffix)) {
      
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        geom_hline(linetype = "dashed", yintercept = 70, color = "darkgreen") +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)") +
        annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen") +
        scale_x_continuous(breaks = weightseq)
      
    } else if (!grepl("Access", suffix)&!grepl("inpatient", suffix)) {
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)") +
        scale_x_continuous(breaks = weightseq)
      
    } else {
      
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        geom_hline(linetype = "dashed", yintercept = 70, color = "darkgreen") +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)") +
        annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen") +
        scale_x_continuous(breaks = weightseq)
    }
    
  } else {
    
    if (grepl("Access", suffix)&!grepl("inpatient", suffix)) {
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        geom_hline(linetype = "dashed", yintercept = 70, color = "darkgreen") +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)") +
        annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen")
      
    } else if (!grepl("Access", suffix)&!grepl("inpatient", suffix)) {
      
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)")
      
    } else {
      
      susplot <- ggplot(df, aes(x = Weight, y = Percentage, group = Metric, color = Metric, fill = Metric)) +
        geom_area(position = "identity", alpha = 0.6) +
        ylim(0, 100) +
        theme_minimal() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        ggtitle(glue("{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Urine pathogens covered (%)") +
        scale_x_continuous(breaks = weightseq)
      
    }
  }
  
  
  
  ggsave(glue("{suffix}.pdf"), plot = susplot, device = "pdf", width = 10, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  susplot <<-susplot
  
  print(susplot)
  
}

###Prioritisation by treatment utility
util_mk1 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(Rx_utility)) %>% select(Antimicrobial,Rx_utility) %>% 
    mutate(Rx_utility = round(Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "Rx_utility")
  
}

###Prioritisation by treatment utility
util_mk1a = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(imp_Rx_utility)) %>% select(Antimicrobial,imp_Rx_utility) %>% 
    mutate(imp_Rx_utility = round(imp_Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`imp_Rx Utility` = "imp_Rx_utility")
  
}

###Prioritisation by AST utility
util_mk2 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(AST_utility)) %>% select(Antimicrobial,AST_utility) %>% 
    mutate(AST_utility = round(AST_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`AST Utility` = "AST_utility")
  
}

###Prioritisation by Intravenous utility
util_mk3 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>% 
    arrange(desc(Urosepsis_Rx_utility)) %>% select(Antimicrobial,Urosepsis_Rx_utility) %>% 
    mutate(Urosepsis_Rx_utility = round(Urosepsis_Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Intravenous Rx Utility` = "Urosepsis_Rx_utility")
  
}

###Prioritisation by Oral utility
util_mk4 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(Outpatient_Rx_utility)) %>% select(Antimicrobial,Outpatient_Rx_utility) %>% 
    mutate(Outpatient_Rx_utility = round(Outpatient_Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Oral Rx Utility` = "Outpatient_Rx_utility")
  
}

###Prioritisation by Intravenous utility (improved)
util_mk3a = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>% 
    arrange(desc(imp_Urosepsis_Rx_utility)) %>% select(Antimicrobial,imp_Urosepsis_Rx_utility) %>% 
    mutate(imp_Urosepsis_Rx_utility = round(imp_Urosepsis_Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`imp_Intravenous Rx Utility` = "imp_Urosepsis_Rx_utility")
  
}

###Prioritisation by Oral utility (improved)
util_mk4a = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(imp_Outpatient_Rx_utility)) %>% select(Antimicrobial,imp_Outpatient_Rx_utility) %>% 
    mutate(imp_Outpatient_Rx_utility = round(imp_Outpatient_Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`imp_Oral Rx Utility` = "imp_Outpatient_Rx_utility")
  
}

###Assigning treatment recommendations
assign_PDRx <- function(df,probab_df,method_used,ab_list1=ab_singles) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk1(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Intravenous treatment recommendations
assign_Intravenous <- function(df,probab_df,method_used,ab_list1=iv_ab_singles) {
  
  
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk3(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Oral treatment recommendations
assign_Oral <- function(df,probab_df,method_used,ab_list1=oral_ab_singles) {
  
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk4(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning treatment recommendations
assign_PDRxa <- function(df,probab_df,method_used,ab_list1=ab_singles) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk1a(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Intravenous treatment recommendations
assign_Intravenousa <- function(df,probab_df,method_used,ab_list1=iv_ab_singles) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk3a(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(ab_list1)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning Oral treatment recommendations
assign_Orala <- function(df,probab_df,method_used,ab_list1=oral_ab_singles) {
  
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk4a(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
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
  
  data.frame(vec,Metric=label,Weight=weightseq)
  
}

###Antimicrobial name mapping
replace_values <- function(column, map) {
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(map)) map[[x]] else x)
}
reverse_values <- function(column, map) {
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% map) names(map)[map == x] else x)
}

###Illness severity result checker
result_checker <- function(df,probs_df,acuity_score) {
  
  df <- df %>% filter(acuity==acuity_score)  
  probs_df <- probs_df %>% semi_join(
    df,by="micro_specimen_id"
  )
  
  
  urines_abx <- df
  
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='S'|PDIVRx_1_result=='I'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='S'|PDPORx_1_result=='I'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='S'|PDRx_1_result=='I'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter((PDIVRx_1_result=='S'|PDIVRx_1_result=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter((PDPORx_1_result=='S'|PDPORx_1_result=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df))*100
  overall_s_oral <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
  overall_s_iv <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df))*100
  
  urkey <- urines_abx %>% select(micro_specimen_id,PDRx_1) %>% 
    mutate(Antimicrobial=ab_name(PDRx_1) %>% str_replace("/","-")) %>% 
    select(-PDRx_1)
  median_CDI <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost <- probs_df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_percs <- abrx1_df %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  iv_perc <<- iv_perc
  po_perc <<- po_perc
  overall_perc <<- overall_perc
  iv_s_access <<- iv_s_access
  po_s_access <<- po_s_access
  overall_s_access <<- overall_s_access
  overall_s_oral <<- overall_s_oral
  overall_s_iv <<- overall_s_iv
  median_CDI <<- median_CDI
  Q1_CDI <<- Q1_CDI
  Q3_CDI <<- Q3_CDI
  median_tox <<- median_tox
  Q1_tox <<- Q1_tox
  Q3_tox <<- Q3_tox
  median_cost <<- median_cost
  Q1_cost <<- Q1_cost
  Q3_cost <<- Q3_cost
  abrx1_df <<- abrx1_df
  abrx1_percs <<- abrx1_percs
  
  ###2nd-lines
  
  iv_res_2 <- nrow(df %>% filter(PDIVRx_2_result=='S'|PDIVRx_2_result=='I'))
  po_res_2 <- nrow(df %>% filter(PDPORx_2_result=='S'|PDPORx_2_result=='I'))
  overall_res_2 <- nrow(df %>% filter(PDRx_2_result=='S'|PDRx_2_result=='I'))
  
  iv_perc_2 <- iv_res_2/nrow(df)*100
  po_perc_2 <- po_res_2/nrow(df)*100
  overall_perc_2 <- overall_res_2/nrow(df)*100
  
  iv_s_access_2 <- (nrow(df %>% filter((PDIVRx_2_result=='S'|PDIVRx_2_result=='I')&
                                         Intravenous_2 %in% access_singles))/
                      nrow(df)) * 100
  po_s_access_2 <- (nrow(df %>% filter((PDPORx_2_result=='S'|PDPORx_2_result=='I') &
                                         Oral_2 %in% access_singles))/
                      nrow(df))*100
  overall_s_access_2 <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                              PDRx_2 %in% access_singles))/
                           nrow(df))*100
  overall_s_oral_2 <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                            PDRx_2 %in% oral_singles))/
                         nrow(df))*100
  overall_s_iv_2 <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                          PDRx_2 %in% iv_singles))/
                       nrow(df))*100
  
  urkey_2 <- urines_abx %>% select(micro_specimen_id,PDRx_2) %>% 
    mutate(Antimicrobial=ab_name(PDRx_2) %>% str_replace("/","-")) %>% 
    select(-PDRx_2)
  median_CDI_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost_2 <- probs_df %>% semi_join(urkey_2,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df_2 <- df %>% count(PDRx_2) %>% arrange(desc(n))
  abrx1_percs_2 <- abrx1_df_2 %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  iv_perc_2 <<- iv_perc_2
  po_perc_2 <<- po_perc_2
  overall_perc_2 <<- overall_perc_2
  iv_s_access_2 <<- iv_s_access_2
  po_s_access_2 <<- po_s_access_2
  overall_s_access_2 <<- overall_s_access_2
  overall_s_oral_2 <<- overall_s_oral_2
  overall_s_iv_2 <<- overall_s_iv_2
  median_CDI_2 <<- median_CDI_2
  Q1_CDI_2 <<- Q1_CDI_2
  Q3_CDI_2 <<- Q3_CDI_2
  median_tox_2 <<- median_tox_2
  Q1_tox_2 <<- Q1_tox_2
  Q3_tox_2 <<- Q3_tox_2
  median_cost_2 <<- median_cost_2
  Q1_cost_2 <<- Q1_cost_2
  Q3_cost_2 <<- Q3_cost_2
  abrx1_df_2 <<- abrx1_df_2
  abrx1_percs_2 <<- abrx1_percs_2
  
  ###Improved probability predictions
  
  iv_res_imp <- nrow(df %>% filter(PDIVRx_1_result_imp=='S'|PDIVRx_1_result_imp=='I'))
  po_res_imp <- nrow(df %>% filter(PDPORx_1_result_imp=='S'|PDPORx_1_result_imp=='I'))
  overall_res_imp <- nrow(df %>% filter(PDRx_1_result_imp=='S'|PDRx_1_result_imp=='I'))
  
  iv_perc_imp <- iv_res_imp/nrow(df)*100
  po_perc_imp <- po_res_imp/nrow(df)*100
  overall_perc_imp <- overall_res_imp/nrow(df)*100
  
  iv_s_access_imp <- (nrow(df %>% filter((PDIVRx_1_result_imp=='S'|PDIVRx_1_result_imp=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access_imp <- (nrow(df %>% filter((PDPORx_1_result_imp=='S'|PDPORx_1_result_imp=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access_imp <- (nrow(df %>% filter((PDRx_1_result_imp=='S'|PDRx_1_result_imp=='I') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df))*100
  overall_s_oral_imp <- (nrow(df %>% filter((PDRx_1_result_imp=='S'|PDRx_1_result_imp=='I') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
  overall_s_iv_imp <- (nrow(df %>% filter((PDRx_1_result_imp=='S'|PDRx_1_result_imp=='I') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df))*100
  
  urkey_imp <- urines_abx %>% select(micro_specimen_id,imp_PDRx_1) %>% 
    mutate(Antimicrobial=ab_name(imp_PDRx_1) %>% str_replace("/","-")) %>% 
    select(-imp_PDRx_1)
  median_CDI_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost_imp <- probs_df %>% semi_join(urkey_imp,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df_imp <- df %>% count(imp_PDRx_1) %>% arrange(desc(n))
  abrx1_percs_imp <- abrx1_df_imp %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  iv_perc_imp <<- iv_perc_imp
  po_perc_imp <<- po_perc_imp
  overall_perc_imp <<- overall_perc_imp
  iv_s_access_imp <<- iv_s_access_imp
  po_s_access_imp <<- po_s_access_imp
  overall_s_access_imp <<- overall_s_access_imp
  overall_s_oral_imp <<- overall_s_oral_imp
  overall_s_iv_imp <<- overall_s_iv_imp
  median_CDI_imp <<- median_CDI_imp
  Q1_CDI_imp <<- Q1_CDI_imp
  Q3_CDI_imp <<- Q3_CDI_imp
  median_tox_imp <<- median_tox_imp
  Q1_tox_imp <<- Q1_tox_imp
  Q3_tox_imp <<- Q3_tox_imp
  median_cost_imp <<- median_cost_imp
  Q1_cost_imp <<- Q1_cost_imp
  Q3_cost_imp <<- Q3_cost_imp
  abrx1_df_imp <<- abrx1_df_imp
  abrx1_percs_imp <<- abrx1_percs_imp
  
  ###Combinations
  
  iv_res_comb <- nrow(df %>% filter(PDIVRx_comb_result=='S'|PDIVRx_comb_result=='I'))
  po_res_comb <- nrow(df %>% filter(PDPORx_comb_result=='S'|PDPORx_comb_result=='I'))
  overall_res_comb <- nrow(df %>% filter(PDRx_comb_result=='S'|PDRx_comb_result=='I'))
  
  iv_perc_comb <- iv_res_comb/nrow(df)*100
  po_perc_comb <- po_res_comb/nrow(df)*100
  overall_perc_comb <- overall_res_comb/nrow(df)*100
  
  iv_s_access_comb <- (nrow(df %>% filter((PDIVRx_comb_result=='S'|PDIVRx_comb_result=='I')&
                                         Intravenous_comb_1 %in% all_access))/
                      nrow(df)) * 100
  po_s_access_comb <- (nrow(df %>% filter((PDPORx_comb_result=='S'|PDPORx_comb_result=='I') &
                                         Oral_comb_1 %in% all_access))/
                      nrow(df))*100
  overall_s_access_comb <- (nrow(df %>% filter((PDRx_comb_result=='S'|PDRx_comb_result=='I') &
                                              PDRx_comb_1 %in% all_access))/
                           nrow(df))*100
  overall_s_oral_comb <- (nrow(df %>% filter((PDRx_comb_result=='S'|PDRx_comb_result=='I') &
                                            PDRx_comb_1 %in% all_orals))/
                         nrow(df))*100
  overall_s_iv_comb <- (nrow(df %>% filter((PDRx_comb_result=='S'|PDRx_comb_result=='I') &
                                          PDRx_comb_1 %in% all_ivs))/
                       nrow(df))*100
  
  urkey_comb <- urines_abx %>% select(micro_specimen_id,PDRx_comb_1) %>% 
    mutate(Antimicrobial=replace_values(PDRx_comb_1,combined_antimicrobial_map) %>% str_replace("/","-")) %>% 
    select(-PDRx_comb_1)
  median_CDI_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost_comb <- probs_df %>% semi_join(urkey_comb,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df_comb <- df %>% count(PDRx_comb_1) %>% arrange(desc(n))
  abrx1_percs_comb <- abrx1_df_comb %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  iv_perc_comb <<- iv_perc_comb
  po_perc_comb <<- po_perc_comb
  overall_perc_comb <<- overall_perc_comb
  iv_s_access_comb <<- iv_s_access_comb
  po_s_access_comb <<- po_s_access_comb
  overall_s_access_comb <<- overall_s_access_comb
  overall_s_oral_comb <<- overall_s_oral_comb
  overall_s_iv_comb <<- overall_s_iv_comb
  median_CDI_comb <<- median_CDI_comb
  Q1_CDI_comb <<- Q1_CDI_comb
  Q3_CDI_comb <<- Q3_CDI_comb
  median_tox_comb <<- median_tox_comb
  Q1_tox_comb <<- Q1_tox_comb
  Q3_tox_comb <<- Q3_tox_comb
  median_cost_comb <<- median_cost_comb
  Q1_cost_comb <<- Q1_cost_comb
  Q3_cost_comb <<- Q3_cost_comb
  abrx1_df_comb <<- abrx1_df_comb
  abrx1_percs_comb <<- abrx1_percs_comb
  
  ###Antibiotics prescribed
  overall_res_px_abx <- nrow(df %>% filter(Px_Abx_result=='S'|Px_Abx_result=='I'))
  overall_perc_px_abx <- overall_res_px_abx/nrow(df)*100
  
  overall_s_access_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                                   Px_Abx %in% access_modified_abx_map))/
                                nrow(df))*100
  overall_s_oral_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                                 Px_Abx %in% oral_modified_abx_map))/
                              nrow(df))*100
  overall_s_iv_px_abx <- (nrow(df %>% filter((Px_Abx_result=='S'|Px_Abx_result=='I') &
                                               Px_Abx %in% iv_modified_abx_map))/
                            nrow(df))*100
  
  urkey_abx_px <- urines_abx %>% select(micro_specimen_id,Prescribed_abx) %>% 
    mutate(Antimicrobial=Prescribed_abx) %>% 
    select(-Prescribed_abx)
  median_CDI_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost_px_abx <- probs_df %>% semi_join(urkey_abx_px,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df_px_abx <- df %>% count(Px_Abx) %>% arrange(desc(n))
  abrx1_percs_px_abx <- abrx1_df_px_abx %>% mutate(
    Percentage = (n/nrow(df))*100) %>% select(-n)
  
  overall_perc_px_abx <<- overall_perc_px_abx
  overall_s_access_px_abx <<- overall_s_access_px_abx
  overall_s_oral_px_abx <<- overall_s_oral_px_abx
  overall_s_iv_px_abx <<- overall_s_iv_px_abx
  median_CDI_px_abx <<- median_CDI_px_abx
  Q1_CDI_px_abx <<- Q1_CDI_px_abx
  Q3_CDI_px_abx <<- Q3_CDI_px_abx
  median_tox_px_abx <<- median_tox_px_abx
  Q1_tox_px_abx <<- Q1_tox_px_abx
  Q3_tox_px_abx <<- Q3_tox_px_abx
  median_cost_px_abx <<- median_cost_px_abx
  Q1_cost_px_abx <<- Q1_cost_px_abx
  Q3_cost_px_abx <<- Q3_cost_px_abx
  abrx1_df_px_abx <<- abrx1_df_px_abx
  abrx1_percs_px_abx <<- abrx1_percs_px_abx
  
}

###Illness severity sensitivity analysis
illness_sens <- function(df,df2,acuity_score) {
  
  df <- df %>% mutate(acuity=acuity_score) %>% 
    calculate_utilities(R_weight = 1)
    
  urines_abx <- df2
  
  ###Filter out combination predictions not present in training dataset
  abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
    str_replace_all("/","-")
  df <- df %>% filter(Antimicrobial %in% abx_in_train)
  
  ##Individual treatment recommendations
  all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
               "MEM","CIP","GEN","SXT","NIT","VAN")
  long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
  all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  ab_singles <- names(singles_map)
  iv_ab_singles <- names(iv_singles_map)
  oral_ab_singles <- names(oral_singles_map)
  
  ##Overall
  df2 <- df2 %>% assign_PDRx(df,"PDRx_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ##IV
  df2 <- df2 %>% assign_Intravenous(df,"Intravenous_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral
  df2 <- df2 %>% assign_Oral(df,"Oral_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Standard panel treatment & AST recommendations
  df2 <- df2 %>% assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO")
  
  df2 <- df2 %>% mutate(STANDARD_IV_1="TZP",
                                STANDARD_IV_2="MEM",
                                STANDARD_PO_1="NIT",
                                STANDARD_PO_2="SXT")
  
  ##Microsimulation analysis
  df2 <- df2 %>%
    rowwise() %>%
    mutate(PDIVRx_1_result = get(Intravenous_1),
           PDPORx_1_result = get(Oral_1),
           PDRx_1_result = get(PDRx_1)) %>%
    ungroup()
  
  iv_res <- nrow(df2 %>% filter(PDIVRx_1_result=='S'|PDIVRx_1_result=='I'))
  po_res <- nrow(df2 %>% filter(PDPORx_1_result=='S'|PDPORx_1_result=='I'))
  overall_res <- nrow(df2 %>% filter(PDRx_1_result=='S'|PDRx_1_result=='I'))
  
  iv_perc <- iv_res/nrow(df2)*100
  po_perc <- po_res/nrow(df2)*100
  overall_perc <- overall_res/nrow(df2)*100
  
  iv_s_access <- (nrow(df2 %>% filter((PDIVRx_1_result=='S'|PDIVRx_1_result=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df2)) * 100
  po_s_access <- (nrow(df2 %>% filter((PDPORx_1_result=='S'|PDPORx_1_result=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df2))*100
  overall_s_access <- (nrow(df2 %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df2))*100
  overall_s_oral <- (nrow(df2 %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df2))*100
  overall_s_iv <- (nrow(df2 %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df2))*100
  
  urkey <- df2 %>% select(micro_specimen_id,PDRx_1) %>% 
    mutate(Antimicrobial=ab_name(PDRx_1) %>% str_replace("/","-")) %>% 
    select(-PDRx_1)
  median_CDI <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_CDI)) %>% unlist()
  Q1_CDI <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[2]) %>% unlist()
  Q3_CDI <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_CDI)[4]) %>% unlist()
  median_tox <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(prob_tox)) %>% unlist()
  Q1_tox <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[2]) %>% unlist()
  Q3_tox <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(prob_tox)[4]) %>% unlist()
  median_cost <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=median(orig_cost)) %>% unlist()
  Q1_cost <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[2]) %>% unlist()
  Q3_cost <- df %>% semi_join(urkey,by=c("micro_specimen_id","Antimicrobial")) %>%
    summarise(med=quantile(orig_cost)[4]) %>% unlist()
  
  abrx1_df <- df2 %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_percs <- abrx1_df %>% mutate(
    Percentage = (n/nrow(df2))*100) %>% select(-n)
  
  iv_perc_ill <<- iv_perc
  po_perc_ill <<- po_perc
  overall_perc_ill <<- overall_perc
  iv_s_access_ill <<- iv_s_access
  po_s_access_ill <<- po_s_access
  overall_s_access_ill <<- overall_s_access
  overall_s_oral_ill <<- overall_s_oral
  overall_s_iv_ill <<- overall_s_iv
  median_CDI_ill <<- median_CDI
  Q1_CDI_ill <<- Q1_CDI
  Q3_CDI_ill <<- Q3_CDI
  median_tox_ill <<- median_tox
  Q1_tox_ill <<- Q1_tox
  Q3_tox_ill <<- Q3_tox
  median_cost_ill <<- median_cost
  Q1_cost_ill <<- Q1_cost
  Q3_cost_ill <<- Q3_cost
  abrx1_df_ill <<- abrx1_df
  abrx1_percs_ill <<- abrx1_percs
  
}

###Resistance weighting sensitivity analysis (resistance rates increased)
res_sens_analysis_3 <- function(df,probs_df,abx,abcol,a_val,b_val,R_value=1) { 
  
  abcol <- enquo(abcol)
  
  probs_df <- probs_df %>% dist_replace(df,abx,"R","S",a_val,b_val) %>%
    calculate_utilities(R_weight = R_value)
  
  abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
    str_replace_all("/","-")
  probs_df <- probs_df %>% filter(Antimicrobial %in% abx_in_train)
  senskey <- df %>% select(micro_specimen_id,AMP:VAN,AMP_SAM:NIT_VAN)
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
  combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
  iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
              "GEN","SXT","VAN") %>% ab_name() %>% 
    str_replace("/","-") 
  iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  iv_abs <- c(iv_abs, iv_combos)
  oral_abs <- c("AMP","SAM","CIP",
                "SXT","NIT") %>% ab_name() %>% 
    str_replace("/","-")
  oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  oral_abs <- c(oral_abs, oral_combos)
  iv_combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% iv_abs]
  oral_combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% oral_abs]
  ast_combined_antimicrobial_map <- combined_antimicrobial_map[!grepl("_",names(combined_antimicrobial_map))]
  ###Individual treatment recommendations
  all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
               "MEM","CIP","GEN","SXT","NIT","VAN")
  long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
  all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  replace_values <- function(column, map) {
    column %>%
      as.character() %>%
      sapply(function(x) if (x %in% names(map)) map[[x]] else x)
  }
  
  df <- df %>% mutate(!!abcol:=rbern(nrow(df),prob=rbeta(nrow(df),a_val,b_val))) %>% 
    mutate(!!abcol:=case_when(!!abcol==1~"R",TRUE~"S"))
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_",ab_list1 = ab_name(ab_singles) %>% str_replace("/","-")) %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Standard panel treatment & AST recommendations
  df <- df %>% assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO")
  df <- df %>% assign_standard_IV("CRO","TZP","GEN")
  df <- df %>% assign_standard_oral("NIT","SXT","CIP")
  
  urines_abx <- df
  
  ###Reference lists
  all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                   "MEM","CIP","GEN","SXT","NIT","VAN")
  all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
  all_abs <- c(all_singles,all_combos)
  iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                  "GEN","SXT","VAN")
  iv_combos <- combn(iv_singles, 2, FUN = function(x) paste(x, collapse = "_"))
  all_ivs <- c(iv_singles, iv_combos)
  oral_singles <- c("AMP","SAM","CIP",
                    "SXT","NIT")
  oral_combos <- combn(oral_singles, 2, FUN = function(x) paste(x, collapse = "_"))
  all_orals <- c(oral_singles, oral_combos)
  access_singles <- c("AMP","SAM","GEN",
                      "SXT","NIT","CZO")
  access_combos <- combn(access_singles, 2, FUN = function(x) paste(x, collapse = "_"))
  all_access <- c(access_singles, access_combos)
  watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
  watch_combos <- combn(watch_singles, 2, FUN = function(x) paste(x, collapse = "_"))
  all_watch <- c(watch_singles, watch_combos)
  
  df <- df %>%
    rowwise() %>%
    mutate(PDIVRx_1_result = get(Intravenous_1),
           PDIVRx_2_result = get(Intravenous_2),
           PDPORx_1_result = get(Oral_1),
           PDPORx_2_result = get(Oral_2),
           PDRx_1_result = get(PDRx_1),
           PDRx_2_result = get(PDRx_2),
           STIVRx_1_result = get(STANDARD_IV_1),
           STIVRx_2_result = get(STANDARD_IV_2),
           STPORx_1_result = get(STANDARD_PO_1),
           STPORx_2_result = get(STANDARD_PO_2)) %>%
    ungroup()
  
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='S'|PDIVRx_1_result=='I'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='S'|PDPORx_1_result=='I'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='S'|PDRx_1_result=='I'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter((PDIVRx_1_result=='S'|PDIVRx_1_result=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter((PDPORx_1_result=='S'|PDPORx_1_result=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                            PDRx_1 %in% access_singles))/
                         nrow(df))*100
  overall_s_oral <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
  overall_s_iv <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                        PDRx_1 %in% iv_singles))/
                     nrow(df))*100
  
  
  glue("Where an IV agent was required, the personalised approach recommended a switch from a
     prescribed Watch category agent to an effective Access category agent in 
     {round((sum(urines_abx$Wa_Ac_S_IV)/nrow(urines_abx))*100,1)}% (n={sum(urines_abx$Wa_Ac_S_IV)}) of prescriptions.
     
     Where an oral agent was required, the personalised approach recommended a switch from a
     prescribed Watch category agent to an effective Access category agent in 
     {round((sum(urines_abx$Wa_Ac_S_PO)/nrow(urines_abx))*100,1)}% (n={sum(urines_abx$Wa_Ac_S_PO)}) of prescriptions.
     
     Where an IV agent was required, the personalised approach recommended an agent with a subsequent resistant result in
     {round(iv_res/nrow(df)*100,1)}% (n={iv_res}) of cases where urine culture was sent.
     
     Where an oral agent was required, the personalised approach recommended an agent with a subsequent resistant result in
     {round(po_res/nrow(df)*100,1)}% (n={po_res}) of cases where urine culture was sent.
     
     ") %>% print()
  
  iv_perc <<- iv_perc
  po_perc <<- po_perc
  overall_perc <<- overall_perc
  iv_s_access <<- iv_s_access
  po_s_access <<- po_s_access
  overall_s_access <<- overall_s_access
  overall_s_oral <<- overall_s_oral
  overall_s_iv <<- overall_s_iv
  
  print(df %>% count(Intravenous_1) %>% arrange(desc(n)))
  print(df %>% count(Oral_1) %>% arrange(desc(n)))
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  print(df %>% filter(PDRx_1%in%access_singles) %>% nrow())
  print(df %>% filter(PDRx_1%in%oral_singles) %>% nrow())
  print(df %>% filter(PDRx_1%in%iv_singles) %>% nrow())
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
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
  
  ggsave(glue("UF_{left_label}_{right_label}_plot.pdf"), plot = df_plot, device = "pdf", width = 6, height = 6,
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

###Reference lists
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                 "MEM","CIP","GEN","SXT","NIT","VAN")
ab_singles <- all_singles
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                "GEN","SXT","VAN")
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
    sapply(map_combinations, function(x) paste(x, collapse = " & "))
  )
)

singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%all_singles]
iv_singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%iv_singles]
oral_singles_map <- combined_antimicrobial_map[combined_antimicrobial_map%in%oral_singles]
modified_abx_map <- combined_antimicrobial_map
names(modified_abx_map) <- str_replace_all(names(modified_abx_map),
                                           " & ", "_")
modified_abx_map_brief <- modified_abx_map[intersect(names(modified_abx_map),(util_probs_df$Antimicrobial))]
iv_mod_brief <- intersect(modified_abx_map_brief,iv_combos)
oral_mod_brief <- intersect(modified_abx_map_brief,oral_combos)
access_mod_brief <- intersect(modified_abx_map_brief,all_access)

##Cleaning observations information
staykey <- edstays %>% select(stay_id,intime) %>% distinct(stay_id,.keep_all = T) %>% 
  rename(charttime="intime")
triage <- triage %>% left_join(staykey) %>% relocate(charttime,.before = "temperature") %>% 
  mutate(chartdate=as.Date(charttime)) %>% select(subject_id,chartdate,acuity) %>% 
  filter(!is.na(acuity))
ur_util <- ur_util %>% select(-acuity)
ur_util <- ur_util %>% semi_join(triage,by=c("subject_id","chartdate"))
util_probs_df <- util_probs_df %>% semi_join(ur_util,by="micro_specimen_id")
ur_util <- ur_util %>% left_join(triage) %>%
  distinct(micro_specimen_id,.keep_all = T)
acuitykey <- ur_util %>% select(micro_specimen_id,acuity)
util_probs_df <- util_probs_df %>% left_join(acuitykey)

pyxis <- MIMER::clean_antibiotics(pyxis,drug_col=drug)

pyxis <- pyxis %>% filter(is_abx) %>% 
  mutate(chartdate=as.Date(charttime))

pyxis <- pyxis %>% rename(ab_name="abx_name")

pyxis <- pyxis %>% mutate(ab_name=case_when(grepl("Piperacillin",ab_name)~
                                               "Piperacillin/tazobactam",
                                             TRUE~ab_name))
pyxis <- pyxis %>% semi_join(ur_util,by=c("chartdate","subject_id"))

pyxis <- pyxis %>%
  mutate(order = match(ab_name, all_singles %>% ab_name())) %>%
  distinct(subject_id,chartdate,ab_name,.keep_all = T) %>% 
  arrange(subject_id,chartdate, order) %>%
  group_by(subject_id,chartdate) %>%
  summarize(ab_name = str_c(ab_name, collapse = "_"), .groups = "drop") %>% 
  ungroup() %>% 
  mutate(ab_name=str_replace_all(ab_name,"/","-"))

ur_util <- ur_util %>% semi_join(pyxis,by=c("subject_id","chartdate"))
pyxis_key <- pyxis %>% select(ab_name,subject_id) %>% rename(
  Prescribed_abx = "ab_name"
)
ur_util <- ur_util %>% left_join(pyxis_key)
modified_abx_map <- combined_antimicrobial_map
names(modified_abx_map) <- str_replace_all(names(modified_abx_map),
                                           " & ", "_")
iv_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_ivs]
oral_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_orals]
access_modified_abx_map <- modified_abx_map[modified_abx_map%in%all_access]

ur_util <- ur_util %>% filter(Prescribed_abx%in%names(modified_abx_map))
ur_util <- ur_util %>% mutate(Px_Abx = replace_values(Prescribed_abx,modified_abx_map))

util_probs_df <- util_probs_df %>% semi_join(ur_util,by="micro_specimen_id")

##Write to CSV
write_csv(util_probs_df,"util_probs_df_pre_util.csv")
write_csv(ur_util,"ur_util_pre_recs.csv")

##Sensitivity analysis for illness severity score
iv_perclist_ill <- c()
po_perclist_ill <- c()
overall_perclist_ill <- c()
ivac_list_ill <- c()
poac_list_ill <- c()
ovac_list_ill <- c()
ovor_list_ill <- c()
oviv_list_ill <- c()
pdrx1_list_ill <- data.frame(PDRx_1=ab_singles)
weightseq <- c()

overall_perc_ill <- 1
last_perc_ill <- 0
weight_iter <- 0
last_abrx_ill <- tibble(1)
this_abrx_ill <- tibble(0)

while (!identical(last_abrx_ill,this_abrx_ill)) {
  
  last_abrx_ill <- this_abrx_ill
  last_perc_ill <- overall_perc_ill
  
  illness_sens(util_probs_df,ur_util,acuity_score = weight_iter)
  
  iv_perclist_ill <- c(iv_perclist_ill,iv_perc_ill)
  po_perclist_ill <- c(po_perclist_ill,po_perc_ill)
  overall_perclist_ill <- c(overall_perclist_ill,overall_perc_ill)
  ivac_list_ill <- c(ivac_list_ill,iv_s_access_ill)
  poac_list_ill <- c(poac_list_ill,po_s_access_ill)
  ovac_list_ill <- c(ovac_list_ill,overall_s_access_ill)
  ovor_list_ill <- c(ovor_list_ill,overall_s_oral_ill)
  oviv_list_ill <- c(oviv_list_ill,overall_s_iv_ill)
  pdrx1_list_ill <- pdrx1_list_ill %>% left_join(abrx1_df_ill,by="PDRx_1")
  colnames(pdrx1_list_ill)[weight_iter+2] <- weight_iter
  
  this_abrx_ill <- abrx1_df_ill
  
  weightseq <- c(weightseq,weight_iter)
  weight_iter <- weight_iter+1
  print(weight_iter)
  print(overall_perc_ill)
  
}

iv_perclist_ill <- iv_perclist_ill %>% label_binder("All agents")
po_perclist_ill <- po_perclist_ill %>% label_binder("All agents")
overall_perclist_ill <- overall_perclist_ill %>% label_binder("All agents")
ivac_list_ill <- ivac_list_ill %>% label_binder("Access agents")
poac_list_ill <- poac_list_ill %>% label_binder("Access agents")
ovac_list_ill <- ovac_list_ill %>% label_binder("Access agents")
ovor_list_ill <- ovor_list_ill %>% label_binder("Oral agents")
oviv_list_ill <- oviv_list_ill %>% label_binder("IV agents")
iv_perclist_ill <- iv_perclist_ill %>% rename(Percentage = "vec")
po_perclist_ill <- po_perclist_ill %>% rename(Percentage = "vec")
overall_perclist_ill <- overall_perclist_ill %>% rename(Percentage = "vec")
ivac_list_ill <- ivac_list_ill %>% rename(Percentage = "vec")
poac_list_ill <- poac_list_ill %>% rename(Percentage = "vec")
ovac_list_ill <- ovac_list_ill %>% rename(Percentage = "vec")
ovor_list_ill <- ovor_list_ill %>% rename(Percentage = "vec")
oviv_list_ill <- oviv_list_ill %>% rename(Percentage = "vec")

iv_xg_plot_df_ill <- data.frame(rbind(
  iv_perclist_ill,ivac_list_ill
))
po_xg_plot_df_ill <- data.frame(rbind(
  po_perclist_ill,poac_list_ill
))
overall_xg_plot_df_ill <- data.frame(rbind(
  overall_perclist_ill,ovac_list_ill,ovor_list_ill,oviv_list_ill
))

pdrx1_list_ill[is.na(pdrx1_list_ill)] <- 0
pdrx1_df_ill <- pdrx1_list_ill %>% filter(rowSums(select(.,2:ncol(pdrx1_list_ill)))!=0)
pdrx1_df_ill <- melt(pdrx1_df_ill)
colnames(pdrx1_df_ill) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df_ill <- pdrx1_df_ill %>% mutate(Antimicrobial = ab_name(Antimicrobial))
pdrx1_df_ill <- pdrx1_df_ill %>% mutate(`Percentage of first-line recommendations`=
                                  (`Percentage of first-line recommendations`/nrow(ur_util))*100)

write_csv(pdrx1_df_ill,"abplot_df_ill.csv")

abplot_ill <- ggplot(pdrx1_df_ill,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Illness severity sensitivity analysis for first-line recommendations")+
  scale_x_discrete(breaks = weightseq)+
  ylim(0,100)
abplot_ill
ggsave(glue("illness_abplot_ill.pdf"), plot = abplot_ill, device = "pdf", width = 10, height = 6,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

write_csv(iv_xg_plot_df_ill,"iv_xg_plot_df_ill.csv")
write_csv(po_xg_plot_df_ill,"po_xg_plot_df_ill.csv")
write_csv(overall_xg_plot_df_ill,"overall_xg_plot_df_ill.csv")

overall_xg_plot_iv_ill <- overall_xg_plot_df_ill %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv_ill %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                            agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                            suffix="Utility function IV activity against pathogen by ascending illness severity")
overall_xg_plot_oral_ill <- overall_xg_plot_df_ill %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral_ill %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="Utility function oral activity against pathogen by ascending illness severity")
overall_xg_plot_access_ill <- overall_xg_plot_df_ill %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access_ill %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="Utility function Access activity against pathogen by ascending illness severity")

##Nitrofurantoin resistance sensitivity analysis
weightseq <- c()
iv_perclist_nit <- c()
po_perclist_nit <- c()
overall_perclist_nit <- c()
ivac_list_nit <- c()
poac_list_nit <- c()
ovac_list_nit <- c()
oviv_list_nit <- c()
ovor_list_nit <- c()
pdrx1_list_nit <- data.frame(PDRx_1=ab_singles)

a <- 0.5
b <- 18.5

for(i in 1:9) {
  
  densy <- util_probs_df %>% dens_check("Nitrofurantoin",R,a,b)
  print(densy)
  
  res_sens_analysis_3(ur_util,util_probs_df,abx="Nitrofurantoin",NIT,a_val=a,b_val=b)
  
  weightseq[i] <- util_probs_df %>% dist_replace(ur_util,"Nitrofurantoin","R","S",a,b) %>%
    filter(Antimicrobial=="Nitrofurantoin") %>% summarize(medR=median(R)) %>% unlist()
  iv_perclist_nit <- c(iv_perclist_nit,iv_perc)
  po_perclist_nit <- c(po_perclist_nit,po_perc)
  overall_perclist_nit <- c(overall_perclist_nit,overall_perc)
  ivac_list_nit <- c(ivac_list_nit,iv_s_access)
  poac_list_nit <- c(poac_list_nit,po_s_access)
  ovac_list_nit <- c(ovac_list_nit,overall_s_access)
  ovor_list_nit <- c(ovor_list_nit,overall_s_oral)
  oviv_list_nit <- c(oviv_list_nit,overall_s_iv)
  pdrx1_list_nit <- pdrx1_list_nit %>% left_join(abrx1_df,by="PDRx_1")
  colnames(pdrx1_list_nit)[i+1] <- weightseq[i]
  
  a <- a+i
  b <- b-(i/2)
  
}

iv_perclist_nit <- iv_perclist_nit %>% label_binder("All agents")
po_perclist_nit <- po_perclist_nit %>% label_binder("All agents")
overall_perclist_nit <- overall_perclist_nit %>% label_binder("All agents")
ivac_list_nit <- ivac_list_nit %>% label_binder("Access agents")
poac_list_nit <- poac_list_nit %>% label_binder("Access agents")
ovac_list_nit <- ovac_list_nit %>% label_binder("Access agents")
ovor_list_nit <- ovor_list_nit %>% label_binder("Oral agents")
oviv_list_nit <- oviv_list_nit %>% label_binder("IV agents")
iv_perclist_nit <- iv_perclist_nit %>% rename(Percentage = "vec")
po_perclist_nit <- po_perclist_nit %>% rename(Percentage = "vec")
overall_perclist_nit <- overall_perclist_nit %>% rename(Percentage = "vec")
ivac_list_nit <- ivac_list_nit %>% rename(Percentage = "vec")
poac_list_nit <- poac_list_nit %>% rename(Percentage = "vec")
ovac_list_nit <- ovac_list_nit %>% rename(Percentage = "vec")
ovor_list_nit <- ovor_list_nit %>% rename(Percentage = "vec")
oviv_list_nit <- oviv_list_nit %>% rename(Percentage = "vec")

pdrx1_list_nit[is.na(pdrx1_list_nit)] <- 0

pdrx1_df_nit <- pdrx1_list_nit %>% filter(rowSums(select(.,2:ncol(pdrx1_list_nit)))!=0)
pdrx1_df_nit <- melt(pdrx1_df_nit)
colnames(pdrx1_df_nit) <- c("Antimicrobial","Median probability of nitrofurantoin resistance","Percentage of first-line recommendations")
pdrx1_df_nit <- pdrx1_df_nit %>% mutate(Antimicrobial = ab_name(Antimicrobial))
pdrx1_df_nit <- pdrx1_df_nit %>% mutate(`Percentage of first-line recommendations`=
                                          (`Percentage of first-line recommendations`/nrow(ur_util))*100)

pdrx1_df_nit$`Median probability of nitrofurantoin resistance` <- as.numeric(as.character(pdrx1_df_nit$`Median probability of nitrofurantoin resistance`))

write_csv(pdrx1_df_nit,"nitplot_df.csv")

nitplot <- ggplot(pdrx1_df_nit,aes(x=`Median probability of nitrofurantoin resistance`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("First-line automated recommendations according to probability of nitrofurantoin resistance")+
  ylim(0,100)+
  xlim(0,1)

ggsave(glue("nitres_abplot.pdf"), plot = nitplot, device = "pdf", width = 10, height = 6,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

print(nitplot)

overall_plot_df_nit <- data.frame(rbind(
  overall_perclist_nit,ovac_list_nit,ovor_list_nit,oviv_list_nit
))

write_csv(overall_plot_df_nit,"overall_plot_df_nit.csv")

overall_nit_plot_iv <- overall_plot_df_nit %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_nit_plot_iv %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Nitrofurantoin median resistance probability",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                             agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                             suffix="Utility function IV activity against pathogen by nitrofurantoin resistance rate")
overall_nit_plot_oral <- overall_plot_df_nit %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_nit_plot_oral %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Nitrofurantoin median resistance probability",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                               agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                               suffix="Utility function oral activity against pathogen by nitrofurantoin resistance rate")
overall_nit_plot_access <- overall_plot_df_nit %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_nit_plot_access %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Nitrofurantoin median resistance probability",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                 agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                 suffix="Utility function Access activity against pathogen by nitrofurantoin resistance rate")

##Improved probability predictions
imp_cols <- ur_util %>% select(c(micro_specimen_id,AMP:VAN,AMP_SAM:NIT_VAN))
colnames(imp_cols) <- reverse_values(colnames(imp_cols),modified_abx_map)
imp_cols[imp_cols=="S"|imp_cols=="I"] <- "1"
imp_cols[imp_cols=="R"|imp_cols=="NT"] <- "0"
util_probs_df <- util_probs_df %>% left_join(imp_cols,by="micro_specimen_id")
util_probs_df <- util_probs_df %>% rowwise() %>% mutate(
  ab_result_s = as.numeric(get(Antimicrobial)),
  ab_result_r = 1-as.numeric(get(Antimicrobial))
) %>% ungroup() %>% mutate(
  imp_S=(1+ab_result_s)/2,
  imp_R=(1+ab_result_r)/2
)

##Calculate utilities
ur_util <- ur_util %>% mutate(acuity=5-acuity) %>% 
  mutate(acuity=acuity-1)
util_probs_df <- util_probs_df %>% mutate(acuity=5-acuity) %>% 
  mutate(acuity=acuity-1) 

util_probs_df <- util_probs_df %>% 
  calculate_utilities_illsens(R_weight = 1)

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(ab_name) %>% unlist() %>% 
  str_replace_all("/","-")
util_probs_df <- util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

###Dataframe for formulary agent sensitivity analysis
form_util_probs_df <- util_probs_df %>% calculate_utilities(
  formulary_list = c("Ceftriaxone","Ciprofloxacin"))
form_util_probs_df <- form_util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

##Utility analysis

###Utility data visualisation
util_probs_df %>% group_by(Antimicrobial) %>% filter(Antimicrobial%in%names(singles_map)) %>% 
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

##Individual treatment recommendations
all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN")
long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")
all_combos <- combn(all_abs, 2, FUN = function(x) paste(x, collapse = "_"))
ab_singles <- names(singles_map)
iv_ab_singles <- names(iv_singles_map)
oral_ab_singles <- names(oral_singles_map)

##Overall
ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_") %>% 
  mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))

##IV
ur_util <- ur_util %>% assign_Intravenous(util_probs_df,"Intravenous_") %>% 
  mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))

###Oral
ur_util <- ur_util %>% assign_Oral(util_probs_df,"Oral_") %>% 
  mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))

##Combinations overall
ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_comb_",names(modified_abx_map_brief)) %>% 
  mutate(across(starts_with("PDRx_comb_"), ~ replace_values(., modified_abx_map_brief)))

##IV
ur_util <- ur_util %>% assign_Intravenous(util_probs_df,"Intravenous_comb_",names(iv_mod_brief)) %>% 
  mutate(across(starts_with("Intravenous_comb_"), ~ replace_values(., iv_mod_brief)))

###Oral
ur_util <- ur_util %>% assign_Oral(util_probs_df,"Oral_comb_",names(oral_mod_brief)) %>% 
  mutate(across(starts_with("Oral_comb_"), ~ replace_values(., oral_mod_brief)))

##Overall (improved)
ur_util <- ur_util %>% assign_PDRxa(util_probs_df,"imp_PDRx_") %>% 
  mutate(across(starts_with("imp_PDRx_"), ~ replace_values(., combined_antimicrobial_map)))

##IV (improved)
ur_util <- ur_util %>% assign_Intravenousa(util_probs_df,"imp_Intravenous_") %>% 
  mutate(across(starts_with("imp_Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))

###Oral (improved)
ur_util <- ur_util %>% assign_Orala(util_probs_df,"imp_Oral_") %>% 
  mutate(across(starts_with("imp_Oral_"), ~ replace_values(., combined_antimicrobial_map)))

###Standard panel treatment & AST recommendations
ur_util <- ur_util %>% assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO")

ur_util <- ur_util %>% mutate(STANDARD_IV_1="TZP",
                  STANDARD_IV_2="MEM",
                  STANDARD_PO_1="NIT",
                  STANDARD_PO_2="SXT")

##Microsimulation analysis
ur_util <- ur_util %>%
  rowwise() %>%
  mutate(PDIVRx_1_result = get(Intravenous_1),
         PDIVRx_2_result = get(Intravenous_2),
         PDPORx_1_result = get(Oral_1),
         PDPORx_2_result = get(Oral_2),
         PDRx_1_result = get(PDRx_1),
         PDRx_2_result = get(PDRx_2),
         STIVRx_1_result = get(STANDARD_IV_1),
         STIVRx_2_result = get(STANDARD_IV_2),
         STPORx_1_result = get(STANDARD_PO_1),
         STPORx_2_result = get(STANDARD_PO_2),
         Px_Abx_result = get(Px_Abx),
         PDRx_comb_result=get(PDRx_comb_1),
         PDIVRx_comb_result=get(Intravenous_comb_1),
         PDPORx_comb_result=get(Oral_comb_1),
         PDRx_1_result_imp = get(imp_PDRx_1),
         PDIVRx_1_result_imp = get(imp_Intravenous_1),
         PDPORx_1_result_imp = get(imp_Oral_1)) %>%
  ungroup()

##Write to CSVs
write_csv(util_probs_df,"util_probs_df_final.csv")
write_csv(ur_util,"ur_util_final.csv")

###Significance tests
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


###Illness severity primary analysis
iv_perclist <- c()
po_perclist <- c()
overall_perclist <- c()
ivac_list <- c()
poac_list <- c()
ovac_list <- c()
ovor_list <- c()
oviv_list <- c()
pdrx1_list <- data.frame(PDRx_1=as.ab(ab_singles))
percs_list <- data.frame(PDRx_1=as.ab(ab_singles))

iv_perclist_2 <- c()
po_perclist_2 <- c()
overall_perclist_2 <- c()
ivac_list_2 <- c()
poac_list_2 <- c()
ovac_list_2 <- c()
ovor_list_2 <- c()
oviv_list_2 <- c()
pdrx1_list_2 <- data.frame(PDRx_2=as.ab(ab_singles))
percs_list_2 <- data.frame(PDRx_2=as.ab(ab_singles))

iv_perclist_imp <- c()
po_perclist_imp <- c()
overall_perclist_imp <- c()
ivac_list_imp <- c()
poac_list_imp <- c()
ovac_list_imp <- c()
ovor_list_imp <- c()
oviv_list_imp <- c()
pdrx1_list_imp <- data.frame(imp_PDRx_1=as.ab(ab_singles))
percs_list_imp <- data.frame(imp_PDRx_1=as.ab(ab_singles))

iv_perclist_comb <- c()
po_perclist_comb <- c()
overall_perclist_comb <- c()
ivac_list_comb <- c()
poac_list_comb <- c()
ovac_list_comb <- c()
ovor_list_comb <- c()
oviv_list_comb <- c()
pdrx1_list_comb <- data.frame(PDRx_comb_1=modified_abx_map_brief %>% unlist())
percs_list_comb <- data.frame(PDRx_comb_1=modified_abx_map_brief %>% unlist())

cam <- combined_antimicrobial_map %>% unlist()
names(cam) <- NULL
overall_perclist_px_abx <- c()
ovac_list_px_abx <- c()
ovor_list_px_abx <- c()
oviv_list_px_abx <- c()
pdrx1_list_px_abx <- data.frame(Px_Abx=cam)
percs_list_px_abx <- data.frame(Px_Abx=cam)
weightseq <- seq(0,3)

cdi_medlist <- c()
cdi_q1list <- c()
cdi_q3list <- c()
tox_medlist <- c()
tox_q1list <- c()
tox_q3list <- c()
cost_medlist <- c()
cost_q1list <- c()
cost_q3list <- c()
cdi_medlist_2 <- c()
cdi_q1list_2 <- c()
cdi_q3list_2 <- c()
tox_medlist_2 <- c()
tox_q1list_2 <- c()
tox_q3list_2 <- c()
cost_medlist_2 <- c()
cost_q1list_2 <- c()
cost_q3list_2 <- c()
cdi_medlist_px_abx <- c()
cdi_q1list_px_abx <- c()
cdi_q3list_px_abx <- c()
tox_medlist_px_abx <- c()
tox_q1list_px_abx <- c()
tox_q3list_px_abx <- c()
cost_medlist_px_abx <- c()
cost_q1list_px_abx <- c()
cost_q3list_px_abx <- c()

for(weight in seq_along(weightseq)) {
  
  result_checker(ur_util,util_probs_df,weightseq[weight])
  
  iv_perclist <- c(iv_perclist,iv_perc)
  po_perclist <- c(po_perclist,po_perc)
  overall_perclist <- c(overall_perclist,overall_perc)
  ivac_list <- c(ivac_list,iv_s_access)
  poac_list <- c(poac_list,po_s_access)
  ovac_list <- c(ovac_list,overall_s_access)
  ovor_list <- c(ovor_list,overall_s_oral)
  oviv_list <- c(oviv_list,overall_s_iv)
  pdrx1_list <- pdrx1_list %>% left_join(abrx1_df,by="PDRx_1")
  colnames(pdrx1_list)[weightseq[weight]+2] <- weightseq[weight]
  percs_list <- percs_list %>% left_join(abrx1_percs,by="PDRx_1")
  colnames(percs_list)[weightseq[weight]+2] <- weightseq[weight]
  
  iv_perclist_2 <- c(iv_perclist_2,iv_perc_2)
  po_perclist_2 <- c(po_perclist_2,po_perc_2)
  overall_perclist_2 <- c(overall_perclist_2,overall_perc_2)
  ivac_list_2 <- c(ivac_list_2,iv_s_access_2)
  poac_list_2 <- c(poac_list_2,po_s_access_2)
  ovac_list_2 <- c(ovac_list_2,overall_s_access_2)
  ovor_list_2 <- c(ovor_list_2,overall_s_oral_2)
  oviv_list_2 <- c(oviv_list_2,overall_s_iv_2)
  pdrx1_list_2 <- pdrx1_list_2 %>% left_join(abrx1_df_2,by="PDRx_2")
  colnames(pdrx1_list_2)[weightseq[weight]+2] <- weightseq[weight]
  percs_list_2 <- percs_list_2 %>% left_join(abrx1_percs_2,by="PDRx_2")
  colnames(percs_list_2)[weightseq[weight]+2] <- weightseq[weight]
  
  iv_perclist_imp <- c(iv_perclist_imp,iv_perc_imp)
  po_perclist_imp <- c(po_perclist_imp,po_perc_imp)
  overall_perclist_imp <- c(overall_perclist_imp,overall_perc_imp)
  ivac_list_imp <- c(ivac_list_imp,iv_s_access_imp)
  poac_list_imp <- c(poac_list_imp,po_s_access_imp)
  ovac_list_imp <- c(ovac_list_imp,overall_s_access_imp)
  ovor_list_imp <- c(ovor_list_imp,overall_s_oral_imp)
  oviv_list_imp <- c(oviv_list_imp,overall_s_iv_imp)
  pdrx1_list_imp <- pdrx1_list_imp %>% left_join(abrx1_df_imp,by="imp_PDRx_1")
  colnames(pdrx1_list_imp)[weightseq[weight]+2] <- weightseq[weight]
  percs_list_imp <- percs_list_imp %>% left_join(abrx1_percs_imp,by="imp_PDRx_1")
  colnames(percs_list_imp)[weightseq[weight]+2] <- weightseq[weight]
  
  iv_perclist_comb <- c(iv_perclist_comb,iv_perc_comb)
  po_perclist_comb <- c(po_perclist_comb,po_perc_comb)
  overall_perclist_comb <- c(overall_perclist_comb,overall_perc_comb)
  ivac_list_comb <- c(ivac_list_comb,iv_s_access_comb)
  poac_list_comb <- c(poac_list_comb,po_s_access_comb)
  ovac_list_comb <- c(ovac_list_comb,overall_s_access_comb)
  ovor_list_comb <- c(ovor_list_comb,overall_s_oral_comb)
  oviv_list_comb <- c(oviv_list_comb,overall_s_iv_comb)
  pdrx1_list_comb <- pdrx1_list_comb %>% left_join(abrx1_df_comb,by="PDRx_comb_1")
  colnames(pdrx1_list_comb)[weightseq[weight]+2] <- weightseq[weight]
  percs_list_comb <- percs_list_comb %>% left_join(abrx1_percs_comb,by="PDRx_comb_1")
  colnames(percs_list_comb)[weightseq[weight]+2] <- weightseq[weight]

  overall_perclist_px_abx <- c(overall_perclist_px_abx,overall_perc_px_abx)
  ovac_list_px_abx <- c(ovac_list_px_abx,overall_s_access_px_abx)
  ovor_list_px_abx <- c(ovor_list_px_abx,overall_s_oral_px_abx)
  oviv_list_px_abx <- c(oviv_list_px_abx,overall_s_iv_px_abx)
  pdrx1_list_px_abx <- pdrx1_list_px_abx %>% left_join(abrx1_df_px_abx,by="Px_Abx")
  colnames(pdrx1_list_px_abx)[weightseq[weight]+2] <- weightseq[weight]
  percs_list_px_abx <- percs_list_px_abx %>% left_join(abrx1_percs_px_abx,by="Px_Abx")
  colnames(percs_list_px_abx)[weightseq[weight]+2] <- weightseq[weight]
  
  cdi_medlist <- c(cdi_medlist,median_CDI)
  cdi_q1list <- c(cdi_q1list,Q1_CDI)
  cdi_q3list <- c(cdi_q3list,Q3_CDI)
  tox_medlist <- c(tox_medlist,median_tox)
  tox_q1list <- c(tox_q1list,Q1_tox)
  tox_q3list <- c(tox_q3list,Q3_tox)
  cost_medlist <- c(cost_medlist,median_cost)
  cost_q1list <- c(cost_q1list,Q1_cost)
  cost_q3list <- c(cost_q3list,Q3_cost)
  
  cdi_medlist_2 <- c(cdi_medlist_2,median_CDI_2)
  cdi_q1list_2 <- c(cdi_q1list_2,Q1_CDI_2)
  cdi_q3list_2 <- c(cdi_q3list_2,Q3_CDI_2)
  tox_medlist_2 <- c(tox_medlist_2,median_tox_2)
  tox_q1list_2 <- c(tox_q1list_2,Q1_tox_2)
  tox_q3list_2 <- c(tox_q3list_2,Q3_tox_2)
  cost_medlist_2 <- c(cost_medlist_2,median_cost_2)
  cost_q1list_2 <- c(cost_q1list_2,Q1_cost_2)
  cost_q3list_2 <- c(cost_q3list_2,Q3_cost_2)
  
  cdi_medlist_px_abx <- c(cdi_medlist_px_abx,median_CDI_px_abx)
  cdi_q1list_px_abx <- c(cdi_q1list_px_abx,Q1_CDI_px_abx)
  cdi_q3list_px_abx <- c(cdi_q3list_px_abx,Q3_CDI_px_abx)
  tox_medlist_px_abx <- c(tox_medlist_px_abx,median_tox_px_abx)
  tox_q1list_px_abx <- c(tox_q1list_px_abx,Q1_tox_px_abx)
  tox_q3list_px_abx <- c(tox_q3list_px_abx,Q3_tox_px_abx)
  cost_medlist_px_abx <- c(cost_medlist_px_abx,median_cost_px_abx)
  cost_q1list_px_abx <- c(cost_q1list_px_abx,Q1_cost_px_abx)
  cost_q3list_px_abx <- c(cost_q3list_px_abx,Q3_cost_px_abx)
  
}

iv_perclist <- iv_perclist %>% label_binder("All agents")
po_perclist <- po_perclist %>% label_binder("All agents")
overall_perclist <- overall_perclist %>% label_binder("All agents")
ivac_list <- ivac_list %>% label_binder("Access agents")
poac_list <- poac_list %>% label_binder("Access agents")
ovac_list <- ovac_list %>% label_binder("Access agents")
ovor_list <- ovor_list %>% label_binder("Oral agents")
oviv_list <- oviv_list %>% label_binder("IV agents")
iv_perclist <- iv_perclist %>% rename(Percentage = "vec")
po_perclist <- po_perclist %>% rename(Percentage = "vec")
overall_perclist <- overall_perclist %>% rename(Percentage = "vec")
ivac_list <- ivac_list %>% rename(Percentage = "vec")
poac_list <- poac_list %>% rename(Percentage = "vec")
ovac_list <- ovac_list %>% rename(Percentage = "vec")
ovor_list <- ovor_list %>% rename(Percentage = "vec")
oviv_list <- oviv_list %>% rename(Percentage = "vec")

iv_perclist_2 <- iv_perclist_2 %>% label_binder("All agents")
po_perclist_2 <- po_perclist_2 %>% label_binder("All agents")
overall_perclist_2 <- overall_perclist_2 %>% label_binder("All agents")
ivac_list_2 <- ivac_list_2 %>% label_binder("Access agents")
poac_list_2 <- poac_list_2 %>% label_binder("Access agents")
ovac_list_2 <- ovac_list_2 %>% label_binder("Access agents")
ovor_list_2 <- ovor_list_2 %>% label_binder("Oral agents")
oviv_list_2 <- oviv_list_2 %>% label_binder("IV agents")
iv_perclist_2 <- iv_perclist_2 %>% rename(Percentage = "vec")
po_perclist_2 <- po_perclist_2 %>% rename(Percentage = "vec")
overall_perclist_2 <- overall_perclist_2 %>% rename(Percentage = "vec")
ivac_list_2 <- ivac_list_2 %>% rename(Percentage = "vec")
poac_list_2 <- poac_list_2 %>% rename(Percentage = "vec")
ovac_list_2 <- ovac_list_2 %>% rename(Percentage = "vec")
ovor_list_2 <- ovor_list_2 %>% rename(Percentage = "vec")
oviv_list_2 <- oviv_list_2 %>% rename(Percentage = "vec")

iv_perclist_imp <- iv_perclist_imp %>% label_binder("All agents")
po_perclist_imp <- po_perclist_imp %>% label_binder("All agents")
overall_perclist_imp <- overall_perclist_imp %>% label_binder("All agents")
ivac_list_imp <- ivac_list_imp %>% label_binder("Access agents")
poac_list_imp <- poac_list_imp %>% label_binder("Access agents")
ovac_list_imp <- ovac_list_imp %>% label_binder("Access agents")
ovor_list_imp <- ovor_list_imp %>% label_binder("Oral agents")
oviv_list_imp <- oviv_list_imp %>% label_binder("IV agents")
iv_perclist_imp <- iv_perclist_imp %>% rename(Percentage = "vec")
po_perclist_imp <- po_perclist_imp %>% rename(Percentage = "vec")
overall_perclist_imp <- overall_perclist_imp %>% rename(Percentage = "vec")
ivac_list_imp <- ivac_list_imp %>% rename(Percentage = "vec")
poac_list_imp <- poac_list_imp %>% rename(Percentage = "vec")
ovac_list_imp <- ovac_list_imp %>% rename(Percentage = "vec")
ovor_list_imp <- ovor_list_imp %>% rename(Percentage = "vec")
oviv_list_imp <- oviv_list_imp %>% rename(Percentage = "vec")

iv_perclist_comb <- iv_perclist_comb %>% label_binder("All agents")
po_perclist_comb <- po_perclist_comb %>% label_binder("All agents")
overall_perclist_comb <- overall_perclist_comb %>% label_binder("All agents")
ivac_list_comb <- ivac_list_comb %>% label_binder("Access agents")
poac_list_comb <- poac_list_comb %>% label_binder("Access agents")
ovac_list_comb <- ovac_list_comb %>% label_binder("Access agents")
ovor_list_comb <- ovor_list_comb %>% label_binder("Oral agents")
oviv_list_comb <- oviv_list_comb %>% label_binder("IV agents")
iv_perclist_comb <- iv_perclist_comb %>% rename(Percentage = "vec")
po_perclist_comb <- po_perclist_comb %>% rename(Percentage = "vec")
overall_perclist_comb <- overall_perclist_comb %>% rename(Percentage = "vec")
ivac_list_comb <- ivac_list_comb %>% rename(Percentage = "vec")
poac_list_comb <- poac_list_comb %>% rename(Percentage = "vec")
ovac_list_comb <- ovac_list_comb %>% rename(Percentage = "vec")
ovor_list_comb <- ovor_list_comb %>% rename(Percentage = "vec")
oviv_list_comb <- oviv_list_comb %>% rename(Percentage = "vec")

overall_perclist_px_abx <- overall_perclist_px_abx %>% label_binder("All agents")
ovac_list_px_abx <- ovac_list_px_abx %>% label_binder("Access agents")
ovor_list_px_abx <- ovor_list_px_abx %>% label_binder("Oral agents")
oviv_list_px_abx <- oviv_list_px_abx %>% label_binder("IV agents")
overall_perclist_px_abx <- overall_perclist_px_abx %>% rename(Percentage = "vec")
ovac_list_px_abx <- ovac_list_px_abx %>% rename(Percentage = "vec")
ovor_list_px_abx <- ovor_list_px_abx %>% rename(Percentage = "vec")
oviv_list_px_abx <- oviv_list_px_abx %>% rename(Percentage = "vec")

iv_xg_plot_df <- data.frame(rbind(
  iv_perclist,ivac_list
))
po_xg_plot_df <- data.frame(rbind(
  po_perclist,poac_list
))
overall_xg_plot_df <- data.frame(rbind(
  overall_perclist,ovac_list,ovor_list,oviv_list
))

overall_xg_plot_df_2 <- data.frame(rbind(
  overall_perclist_2,ovac_list_2,ovor_list_2,oviv_list_2
))

overall_xg_plot_df_imp <- data.frame(rbind(
  overall_perclist_imp,ovac_list_imp,ovor_list_imp,oviv_list_imp
))

overall_xg_plot_df_comb <- data.frame(rbind(
  overall_perclist_comb,ovac_list_comb,ovor_list_comb,oviv_list_comb
))

overall_xg_plot_df_px_abx <- data.frame(rbind(
  overall_perclist_px_abx,ovac_list_px_abx,ovor_list_px_abx,oviv_list_px_abx
))



###Recommendations (1st line)
pdrx1_list[is.na(pdrx1_list)] <- 0
percs_list[is.na(percs_list)] <- 0
pdrx1_df <- percs_list %>% filter(rowSums(select(.,2:ncol(percs_list)))!=0)
pdrx1_df <- melt(pdrx1_df)
colnames(pdrx1_df) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df <- pdrx1_df %>% mutate(Antimicrobial = ab_name(Antimicrobial))

write_csv(pdrx1_df,"abplot_df.csv")

abplot <- ggplot(pdrx1_df,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("First-line utility function recommendations according to illness severity")+
  scale_x_discrete(breaks = seq(0,3))+
  ylim(0,100)
abplot
ggsave(glue("illness_abplot.pdf"), plot = abplot, device = "pdf", width = 10, height = 5,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Recommendations (2nd line)
pdrx1_list_2[is.na(pdrx1_list_2)] <- 0
percs_list_2[is.na(percs_list_2)] <- 0
pdrx1_df_2 <- percs_list_2 %>% filter(rowSums(select(.,2:ncol(percs_list_2)))!=0)
pdrx1_df_2 <- melt(pdrx1_df_2)
colnames(pdrx1_df_2) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df_2 <- pdrx1_df_2 %>% mutate(Antimicrobial = ab_name(Antimicrobial))

write_csv(pdrx1_df_2,"abplot_df_2.csv")

abplot_2 <- ggplot(pdrx1_df_2,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Second-line utility function recommendations according to illness severity")+
  scale_x_discrete(breaks = seq(0,3))+
  ylim(0,100)
abplot_2
ggsave(glue("illness_abplot_2.pdf"), plot = abplot_2, device = "pdf", width = 10, height = 5,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Recommendations (improved predictions)
pdrx1_list_imp[is.na(pdrx1_list_imp)] <- 0
percs_list_imp[is.na(percs_list_imp)] <- 0
pdrx1_df_imp <- percs_list_imp %>% filter(rowSums(select(.,2:ncol(percs_list_imp)))!=0)
pdrx1_df_imp <- melt(pdrx1_df_imp)
colnames(pdrx1_df_imp) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df_imp <- pdrx1_df_imp %>% mutate(Antimicrobial = reverse_values(Antimicrobial,combined_antimicrobial_map))

write_csv(pdrx1_df_imp,"abplot_df_imp.csv")

abplot_imp <- ggplot(pdrx1_df_imp,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Utility function recommendations according to illness severity\n(improved probability predictions)")+
  scale_x_discrete(breaks = seq(0,3))+
  ylim(0,100)
abplot_imp
ggsave(glue("illness_abplot_imp.pdf"), plot = abplot_imp, device = "pdf", width = 10, height = 5,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Recommendations (combinations)
pdrx1_list_comb[is.na(pdrx1_list_comb)] <- 0
percs_list_comb[is.na(percs_list_comb)] <- 0
pdrx1_df_comb <- percs_list_comb %>% filter(rowSums(select(.,2:ncol(percs_list_comb)))!=0)
pdrx1_df_comb <- melt(pdrx1_df_comb)
colnames(pdrx1_df_comb) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df_comb <- pdrx1_df_comb %>% mutate(Antimicrobial = reverse_values(Antimicrobial,combined_antimicrobial_map))

write_csv(pdrx1_df_comb,"abplot_df_comb.csv")

abplot_comb <- ggplot(pdrx1_df_comb,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Utility function recommendations according to illness severity\n(combinations included)")+
  scale_x_discrete(breaks = seq(0,3))+
  ylim(0,100)
abplot_comb
ggsave(glue("illness_abplot_comb.pdf"), plot = abplot_comb, device = "pdf", width = 10, height = 5,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Actually prescribed
pdrx1_list_px_abx[is.na(pdrx1_list_px_abx)] <- 0
percs_list_px_abx[is.na(percs_list_px_abx)] <- 0
pdrx1_df_px_abx <- percs_list_px_abx %>% mutate(Px_Abx=case_when((`0`+`1`+`2`+`3`)<15~"Other",
                                              TRUE~Px_Abx))
pdrx1_df_px_abx <- pdrx1_df_px_abx %>% filter(rowSums(select(.,2:ncol(pdrx1_df_px_abx)))!=0)
pdrx1_df_px_abx <- melt(pdrx1_df_px_abx)
colnames(pdrx1_df_px_abx) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df_px_abx <- pdrx1_df_px_abx %>% mutate(Antimicrobial=
                                                reverse_values(Antimicrobial,
                                                               combined_antimicrobial_map)) %>% 
  mutate(Antimicrobial=str_replace_all(Antimicrobial,"-","/"))

pdrx1_df_px_abx <- pdrx1_df_px_abx %>% group_by(Antimicrobial,`Illness severity`) %>% 
  summarise(`Percentage of first-line recommendations`=sum(`Percentage of first-line recommendations`)) %>% 
  ungroup()

write_csv(pdrx1_df_px_abx,"abplot_df_px_abx.csv")

abplot_px_abx <- ggplot(pdrx1_df_px_abx,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Antibiotics prescribed in ED according to illness severity")+
  scale_x_discrete(breaks = seq(0,3))+
  ylim(0,100)
abplot_px_abx
ggsave(glue("illness_abplot_px_abx.pdf"), plot = abplot_px_abx, device = "pdf", width = 10, height = 5,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")
abplot_px_abx


###Proportion plot dataframes

write_csv(iv_xg_plot_df,"iv_xg_plot_df.csv")
write_csv(po_xg_plot_df,"po_xg_plot_df.csv")
write_csv(overall_xg_plot_df,"overall_xg_plot_df.csv")

overall_xg_plot_iv <- overall_xg_plot_df %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                            agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                            suffix="Utility function IV agent activity against pathogen by severity score")
illplot_1 <- susplot
overall_xg_plot_oral <- overall_xg_plot_df %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="Utility function oral agent activity against pathogen by severity score")
illplot_2 <- susplot
overall_xg_plot_access <- overall_xg_plot_df %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="Utility function Access agent activity against pathogen by severity score")
illplot_3 <- susplot

###2nd line coverage
write_csv(overall_xg_plot_df_2,"overall_xg_plot_df_2.csv")

overall_xg_plot_iv_2 <- overall_xg_plot_df_2 %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv_2 %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                            agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                            suffix="Utility function second-line IV activity against pathogen by severity score")
overall_xg_plot_oral_2 <- overall_xg_plot_df_2 %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral_2 %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="Utility function second-line oral activity against pathogen by severity score")
overall_xg_plot_access_2 <- overall_xg_plot_df_2 %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access_2 %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="Utility function second-line Access activity against pathogen by severity score")

##Improved predictions
write_csv(overall_xg_plot_df_imp,"overall_xg_plot_df_imp.csv")

overall_xg_plot_iv_imp <- overall_xg_plot_df_imp %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv_imp %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="Utility function IV activity against pathogen by severity score (improved predictions)")
overall_xg_plot_oral_imp <- overall_xg_plot_df_imp %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral_imp %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="Utility function oral activity against pathogen by severity score (improved predictions)")
overall_xg_plot_access_imp <- overall_xg_plot_df_imp %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access_imp %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                  agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                  suffix="Utility function Access activity against pathogen by severity score (improved predictions)")

###Combination coverage
write_csv(overall_xg_plot_df_comb,"overall_xg_plot_df_comb.csv")

overall_xg_plot_iv_comb <- overall_xg_plot_df_comb %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv_comb %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="Utility function IV combination activity against pathogen by severity score")
overall_xg_plot_oral_comb <- overall_xg_plot_df_comb %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral_comb %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="Utility function oral combination activity against pathogen by severity score")
overall_xg_plot_access_comb <- overall_xg_plot_df_comb %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access_comb %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                  agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                  suffix="Utility function Access combination activity against pathogen by severity score")

###Prescribed abx coverage
write_csv(overall_xg_plot_df_px_abx,"overall_xg_plot_df_px_abx.csv")

overall_xg_plot_iv_px_abx <- overall_xg_plot_df_px_abx %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv_px_abx %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                            agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                            suffix="ED prescription IV agent activity against pathogen by severity score")
illplot_4 <- susplot
overall_xg_plot_oral_px_abx <- overall_xg_plot_df_px_abx %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral_px_abx %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="ED prescription oral agent activity against pathogen by severity score")
illplot_5 <- susplot
overall_xg_plot_access_px_abx <- overall_xg_plot_df_px_abx %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access_px_abx %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity score",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="ED prescription Access agent activity against pathogen by severity score")

##CDI, toxicity, and cost

###First-lines
cdi_df <- data.frame(
  illness_severity = weightseq,
  med=cdi_medlist,
  q1=cdi_q1list,
  q3=cdi_q3list
)
tox_df <- data.frame(
  illness_severity = weightseq,
  med=tox_medlist,
  q1=tox_q1list,
  q3=tox_q3list
)
cost_df <- data.frame(
  illness_severity = weightseq,
  med=cost_medlist,
  q1=cost_q1list,
  q3=cost_q3list
)

cdiplot <- ggplot(cdi_df, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted CDI probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of CDI")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("cdiplot.pdf"), plot = cdiplot, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

toxplot <- ggplot(tox_df, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted toxicity probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of toxicity")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("toxplot.pdf"), plot = toxplot, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

costplot <- ggplot(cost_df, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median cost of recommended agents")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median minimum agent cost (USD)")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("costplot.pdf"), plot = costplot, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

write_csv(cdi_df,"cdi_df.csv")
write_csv(tox_df,"tox_df.csv")
write_csv(cost_df,"cost_df.csv")

###Second-lines
cdi_df_2 <- data.frame(
  illness_severity = weightseq,
  med=cdi_medlist_2,
  q1=cdi_q1list_2,
  q3=cdi_q3list_2
)
tox_df_2 <- data.frame(
  illness_severity = weightseq,
  med=tox_medlist_2,
  q1=tox_q1list_2,
  q3=tox_q3list_2
)
cost_df_2 <- data.frame(
  illness_severity = weightseq,
  med=cost_medlist_2,
  q1=cost_q1list_2,
  q3=cost_q3list_2
)

cdiplot_2 <- ggplot(cdi_df_2, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted CDI probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of CDI")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("cdiplot_2.pdf"), plot = cdiplot_2, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

toxplot_2 <- ggplot(tox_df_2, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted toxicity probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of toxicity")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("toxplot.pdf_2"), plot = toxplot_2, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

costplot_2 <- ggplot(cost_df_2, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median cost of recommended agents")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median minimum agent cost (USD)")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("costplot_2.pdf"), plot = costplot_2, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

write_csv(cdi_df_2,"cdi_df_2.csv")
write_csv(tox_df_2,"tox_df_2.csv")
write_csv(cost_df_2,"cost_df_2.csv")

###Actual
cdi_df_px_abx <- data.frame(
  illness_severity = weightseq,
  med=cdi_medlist_px_abx,
  q1=cdi_q1list_px_abx,
  q3=cdi_q3list_px_abx
)
tox_df_px_abx <- data.frame(
  illness_severity = weightseq,
  med=tox_medlist_px_abx,
  q1=tox_q1list_px_abx,
  q3=tox_q3list_px_abx
)
cost_df_px_abx <- data.frame(
  illness_severity = weightseq,
  med=cost_medlist_px_abx,
  q1=cost_q1list_px_abx,
  q3=cost_q3list_px_abx
)

cdiplot_px_abx <- ggplot(cdi_df_px_abx, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted CDI probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of CDI")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("cdiplot_px_abx.pdf"), plot = cdiplot_px_abx, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

toxplot_px_abx <- ggplot(tox_df_px_abx, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median predicted toxicity probability of recommended agent")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median predicted probability of toxicity")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("toxplot.pdf_px_abx"), plot = toxplot_px_abx, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

costplot_px_abx <- ggplot(cost_df_px_abx, aes(x = illness_severity)) +
  geom_line(aes(y = as.numeric(med))) +
  geom_ribbon(aes(y = as.numeric(med),
                  ymin = as.numeric(q1),
                  ymax = as.numeric(q3)), alpha = 0.3) +
  ggtitle(glue("Effect of varying illness severity on median cost of recommended agents")) +
  xlab(glue("Illness severity")) +
  ylab(glue("Median minimum agent cost (USD)")) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank())

ggsave(glue("costplot_px_abx.pdf"), plot = costplot_px_abx, device = "pdf", width = 8, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

write_csv(cdi_df_px_abx,"cdi_df_px_abx.csv")
write_csv(tox_df_px_abx,"tox_df_px_abx.csv")
write_csv(cost_df_px_abx,"cost_df_px_abx.csv")

##AST performance analysis - number of results per panel

###Number of AST results per panel
ur_util <- ur_util %>% rpp_ast()

###Assemble data frame for dot plot data visualisation
acs_df <- ur_util %>% assemble_dotplot_df()
acs_df %>% count(AWaRe_results)
write_csv(acs_df,"uf_ast_sourcedata_aware_dotplot.csv")

###Dot plot of number of all S results and Access S results per panel
main_aware_plot <- acs_df %>% main_dotplotter("PDRx\nsingle S","Standard\nsingle S","PDRx\nsingle Access S","Standard\nsingle Access S",
                                              "All agents","WHO access agents")

###Dot plot of number of oral S results and IV S results per panel
main_ivo_plot <- acs_df %>% main_dotplotter("PDRx\nsingle IV S","Standard\nsingle IV S","PDRx\nsingle oral S","Standard\nsingle oral S",
                                              "IV agents","Oral agents")

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

###Significance testing
p_acs <- round(wilcox.test(ur_util %>% pull(n_acS_PDRx6), ur_util %>% pull(n_acS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_all <- round(wilcox.test(ur_util %>% pull(n_allS_PDRx6), ur_util %>% pull(n_allS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_pos <- round(wilcox.test(ur_util %>% pull(n_poS_PDRx6), ur_util %>% pull(n_poS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)
p_ivs <- round(wilcox.test(ur_util %>% pull(n_ivS_PDRx6), ur_util %>% pull(n_ivS_standard6), paired = TRUE, conf.int = TRUE)$p.value, 3)

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

write_csv(util_probs_df,"util_probs_df_manuscript.csv")
write_csv(ur_util,"ur_util_manuscript.csv")

##Utility value sensitivity analysis

util_probs_df <- read_csv("util_probs_df_pre_util.csv")
ur_util <- read_csv("ur_util_pre_recs.csv")

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

###Sensitivity analysis varying illness severity
illness_sens_df <- ur_util %>% util_sens_illness(util_probs_df,Rx_utility,0,3)
illness_sens_df %>% util_sens_plot("Treatment","Illness severity score")

###Sensitivity analysis with variation of CDI and toxicity risk
cdi_prob_df <- ur_util %>% cdi_prob_sens(util_probs_df,Rx_utility,"cdi_prob","prob_CDI",cdi_prob,prob_CDI,"CDI")
tox_prob_df <- ur_util %>% tox_prob_sens(util_probs_df,Rx_utility,"tox_prob","prob_tox",tox_prob,prob_tox,"toxicity")

cdi_prob_df %>% dens_sens_plot_2("CDI probability","Treatment",cdi_prob)
tox_prob_df %>% dens_sens_plot_2("toxicity probability","Treatment",tox_prob)

