#EDGE CASE ANALYSIS

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
util_sens_auc_cdi <- function(df,probs_df,uf,min_val,max_val) {
  
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
      
      probs_df_2 <- probs_df_2 %>% mutate(CDI_AUC=weight_sq[i])
      
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
util_sens_auc_tox <- function(df,probs_df,uf,min_val,max_val) {
  
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
      
      probs_df_2 <- probs_df_2 %>% mutate(tox_AUC=weight_sq[i])
      
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
auc_sens_plot <- function(df,measure) {
  
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
      xlim(0,1) +
      ggtitle(glue("Effect of varying {iterabs[i]} prediction\nAUROC on {iterabs[i]} {measure} utility")) +
      xlab(glue("AUROC for {iterabs[i]} susceptibility prediction")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("auc_sensplot_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave(glue("auc_sens_plots_grid.pdf"), plot = grid_plot, width = 20, height = 30)
  
}
auc_sens_plot_2 <- function(df,measure,auc_type) {
  
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
      xlim(0,1) +
      ggtitle(glue("Effect of varying {auc_type} prediction\nAUROC on {iterabs[i]} {measure} utility")) +
      xlab(glue("AUROC for {iterabs[i]} susceptibility prediction")) +
      ylab(glue("{measure} utility")) +
      theme_minimal()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    sens_plots[[i]] <- df_plot
    
    ggsave(glue("auc_{auc_type}_sensplot_{iterabs[i]}_{measure}.pdf"), plot = df_plot, device = "pdf", width = 6, height = 4,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(df_plot)
    
  }
  
  grid_plot <- plot_grid(plotlist = sens_plots, ncol = 3)
  
  ggsave(glue("{auc_type}_auc_sens_plots_grid.pdf"), plot = grid_plot, width = 20, height = 30)
  
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

###Assembling antimicrobial plot dataframe
abs_df_assemble <- function(pd_df,standard_df) {
  
  minuser <- function(df,abx) {
    
    df %>% filter(ind==ab_name(abx)) %>% arrange(Approach) %>% select(values)
    
    if(nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
       abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
      
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
        abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2)
      
    } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
               abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="Standard"){
      
      -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
          abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2))
      
    } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==1 &
               abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
      
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1)
      
    } else {
      
      -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1))
      
    }
    
  }
  
  abs_df <- bind_rows(
    standard_df %>% data.frame() %>% mutate(Approach = "Standard"),
    pd_df %>% data.frame() %>% mutate(Approach = "PDAST")
  ) %>% mutate(ind = ab_name(ind))
  abs_diffs <- map_df(all_singles, function(abs) {
    abs_df %>% 
      minuser(abs) %>% 
      tibble() %>% 
      mutate(ind = ab_name(abs))
  }) %>% 
    mutate(better = if_else(values > 0, "PDAST", "Standard"),
           values2 = abs(values)) %>% 
    left_join(
      abs_df %>% filter(Approach == "PDAST") %>% select(values, ind) %>% rename(PDAST = values), 
      by = "ind"
    ) %>% 
    left_join(
      abs_df %>% filter(Approach == "Standard") %>% select(values, ind) %>% rename(Standard = values), 
      by = "ind"
    ) %>% 
    mutate(
      Standard = if_else(is.na(Standard), 0, Standard),
      values = if_else(better == "PDAST", PDAST + 200, Standard + 200)
    )
  abs_diffs <<- abs_diffs
  standard_levels <- abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) 
  pdast_levels <- abs_df %>% filter(Approach == "PDAST") %>% filter(!ind %in% standard_levels) %>% arrange(values) %>% pull(ind)
  abs_df <- abs_df %>% 
    anti_join(
      abs_df %>% filter(Approach == "Standard") %>% select(ind), 
      by = "ind"
    ) %>% 
    mutate(Approach = "Standard", values = 0) %>% 
    bind_rows(abs_df) %>% 
    mutate(
      Approach = factor(Approach, levels = c("PDAST", "Standard")),
      ind = factor(ind, levels = c(pdast_levels,standard_levels)),
      aware = if_else(ind %in% ab_name(access_singles), "Access", "Watch")
    )
  axiscols <- if_else(
    abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) %in% ab_name(access_singles),
    "seagreen", "darkorange"
  )
  
  axiscols <<- axiscols
  
  abs_df
  
}

###Cleveland dot plot of AST results
cleveland_ab_plot <- function(ab_df,resulttype,addendum="") {
  
  results_by_ab <- ggplot(ab_df,aes(x=ind,y=values))+
    geom_line(aes(group=ind),alpha=0.5)+
    geom_point(aes(color=Approach),size=4) +
    coord_flip() +
    scale_color_manual(values=c("#00BFC4","#F8766D"))+
    geom_text(data = abs_diffs, aes(color = better, 
                                    label = as.character(glue("+{values2}"))),
              size = 3,hjust=0.5) +
    ggtitle(glue("Total number of {resulttype} AST results by antimicrobial agent\n{addendum}"))+
    xlab("") +
    ylab(glue("Total number of {resulttype} results")) +
    theme_minimal() +
    theme(axis.text.y = element_text(
      colour = axiscols))
  
  ggsave(glue("{resulttype}_ab_cleveplot.pdf"), plot = results_by_ab, device = "pdf", width = 10, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  results_by_ab
  
}

###Illness severity weighting sensitivity analysis
res_sens_analysis <- function(df,probs_df,MEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  median_CDI <<- median_CDI
  Q1_CDI <<- Q1_CDI
  Q3_CDI <<- Q3_CDI
  median_tox <<- median_tox
  Q1_tox <<- Q1_tox
  Q3_tox <<- Q3_tox
  median_cost <<- median_cost
  Q1_cost <<- Q1_cost
  Q3_cost <<- Q3_cost
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}

###Illness severity weighting sensitivity analysis (better predictions)
res_sens_analysis_2 <- function(df,probs_df,MEWS_variable=1) { 
  
  probs_df <- probs_df %>% mutate(S=better_S)
  probs_df <- probs_df %>% calculate_utilities(MEWS = MEWS_variable)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}

###Resistance weighting sensitivity analysis (resistance rates increased)
res_sens_analysis_3 <- function(df,probs_df,abx,abcol,a_val,b_val) { 
  
  abcol <- enquo(abcol)
  
  probs_df <- probs_df %>% dist_replace(df,abx,"R","S",a_val,b_val) %>%
    calculate_utilities()
  
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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

###Resistance weighting sensitivity analysis (2nd-line)
res_sens_analysis_4 <- function(df,probs_df,MEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
    str_replace_all("/","-")
  combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
  iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
              "GEN","SXT","NIT","VAN") %>% ab_name() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  
  iv_res <- nrow(df %>% filter(PDIVRx_2_result=='S'|PDIVRx_2_result=='I'))
  po_res <- nrow(df %>% filter(PDPORx_2_result=='S'|PDPORx_2_result=='I'))
  overall_res <- nrow(df %>% filter(PDRx_2_result=='S'|PDRx_2_result=='I'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter((PDIVRx_2_result=='S'|PDIVRx_2_result=='I')&
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter((PDPORx_2_result=='S'|PDPORx_2_result=='I') &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                            PDRx_2 %in% access_singles))/
                         nrow(df))*100
  overall_s_oral <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                          PDRx_2 %in% oral_singles))/
                       nrow(df))*100
  overall_s_iv <- (nrow(df %>% filter((PDRx_2_result=='S'|PDRx_2_result=='I') &
                                        PDRx_2 %in% iv_singles))/
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
  
  print(df %>% count(PDRx_2) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_2) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}

###Susceptibility plots
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
        annotate("text", x = Inf, y = (nrow(df2 %>% filter(!!agent_col1 == "S" | !!agent_col1 == "I")) / nrow(df2)) * 100 + 2, label = glue(agent_name1), hjust = 1.1, color = "gray16") +
        annotate("text", x = Inf, y = (nrow(df2 %>% filter(!!agent_col2 == "S" | !!agent_col2 == "I")) / nrow(df2)) * 100 + 2, label = glue(agent_name2), hjust = 1.1, color = "gray16")
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
        ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}")) +
        xlab(glue("{variable}")) +
        ylab("Percentage of isolates susceptible to recommendation") +
        scale_x_continuous(breaks = weightseq)
      
    }
  }
  
  
  
  ggsave(glue("{subset}_{measure}_{suffix}_susplot.pdf"), plot = susplot, device = "pdf", width = 10, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  susplot <<-susplot
  
  print(susplot)
  
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

###Combination Resistance weighting sensitivity analysis
combo_res_sens_analysis <- function(df,probs_df,MEWS_variable=0,R_value=1) { 
  
  probs_df <- probs_df %>% calculate_utilities(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_",ab_list1 = names(combined_antimicrobial_map)) %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_",ab_list1 = names(iv_combined_antimicrobial_map)) %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., iv_combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_",ab_list1 = names(oral_combined_antimicrobial_map)) %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., oral_combined_antimicrobial_map)))
  
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
                                       Intravenous_1 %in% all_access))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter((PDPORx_1_result=='S'|PDPORx_1_result=='I') &
                                       Oral_1 %in% all_access))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                            PDRx_1 %in% all_access))/
                         nrow(df))*100
  overall_s_oral <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                          PDRx_1 %in% all_orals))/
                       nrow(df))*100
  overall_s_iv <- (nrow(df %>% filter((PDRx_1_result=='S'|PDRx_1_result=='I') &
                                        PDRx_1 %in% all_ivs))/
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
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}

###Prioritisation by treatment utility
util_mk1 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(Rx_utility)) %>% select(Antimicrobial,Rx_utility) %>% 
    mutate(Rx_utility = round(Rx_utility,1)) %>% dplyr::slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "Rx_utility")
  
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

###Assigning treatment recommendations
assign_PDRx <- function(df,probab_df,method_used,ab_list1=ab_singles %>% ab_name() %>% str_replace("/","-")) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk1(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list = ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
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

###Assigning AST recommendations
assign_PDAST <- function(df,probab_df,method_used,ab_list1=ab_singles %>% ab_name() %>% str_replace("/","-")) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk2(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list = ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
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
assign_Intravenous <- function(df,probab_df,method_used,ab_list1=iv_ab_singles %>% ab_name() %>% str_replace("/","-")) {
  
  
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk3(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list=ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
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
assign_Oral <- function(df,probab_df,method_used,ab_list1=oral_ab_singles %>% ab_name() %>% str_replace("/","-")) {
  
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk4(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
                                  ab_list = ab_list1) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
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

###Utility calculation
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

###Utility calculation for sensitivity analyses
calculate_utilities_recall <- function(df,formulary_list=c(),R_weight=1,MEWS=0) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox*tox_recall + util_CDI*CDI_recall + util_iv*MEWS*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util*recall,
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
res_sens_analysis_recall <- function(df,probs_df,MEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities_recall(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  median_CDI <<- median_CDI
  Q1_CDI <<- Q1_CDI
  Q3_CDI <<- Q3_CDI
  median_tox <<- median_tox
  Q1_tox <<- Q1_tox
  Q3_tox <<- Q3_tox
  median_cost <<- median_cost
  Q1_cost <<- Q1_cost
  Q3_cost <<- Q3_cost
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}
calculate_utilities_precision <- function(df,formulary_list=c(),R_weight=1,MEWS=0) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox*tox_precision + util_CDI*CDI_precision + util_iv*MEWS*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util*precision,
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
res_sens_analysis_precision <- function(df,probs_df,MEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities_precision(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  median_CDI <<- median_CDI
  Q1_CDI <<- Q1_CDI
  Q3_CDI <<- Q3_CDI
  median_tox <<- median_tox
  Q1_tox <<- Q1_tox
  Q3_tox <<- Q3_tox
  median_cost <<- median_cost
  Q1_cost <<- Q1_cost
  Q3_cost <<- Q3_cost
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
}
calculate_utilities_accuracy <- function(df,formulary_list=c(),R_weight=1,MEWS=0) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox*tox_accuracy + util_CDI*CDI_accuracy + util_iv*MEWS*R_weight,
                      short_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      S_utility = S*overall_util*accuracy,
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
res_sens_analysis_accuracy <- function(df,probs_df,MEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities_accuracy(MEWS = MEWS_variable,R_weight = R_value)
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
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
  
  df <- df %>% assign_PDRx(probs_df,"PDRx_") %>% 
    mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Intravenous treatment recommendations
  df <- df %>% assign_Intravenous(probs_df,"Intravenous_") %>% 
    mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Oral treatment recommendations
  df <- df %>% assign_Oral(probs_df,"Oral_") %>% 
    mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))
  
  ###Individual AST recommendations
  df <- df %>% assign_PDAST(probs_df,"PDAST_") %>%
    mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
  
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
  median_CDI <<- median_CDI
  Q1_CDI <<- Q1_CDI
  Q3_CDI <<- Q3_CDI
  median_tox <<- median_tox
  Q1_tox <<- Q1_tox
  Q3_tox <<- Q3_tox
  median_cost <<- median_cost
  Q1_cost <<- Q1_cost
  Q3_cost <<- Q3_cost
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
  abrx1_df <- df %>% count(PDRx_1) %>% arrange(desc(n))
  abrx1_df <<- abrx1_df
  
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
    theme(legend.position = "None",axis.text.y = element_text(
      colour = axiscols))+
    xlab(glue("{application} utility{modification}"))+
    theme()
  
  print(thisplot)
  
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
  
  df$PDAST_rpp_ass <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_singles,"S","I")
  print(1/60)
  df$PDAST_rpp_acs <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_combos,"S","I")
  print(2/60)
  df$PDAST_rpp_aas <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_abs,"S","I")
  print(3/60)
  df$PDAST_rpp_iss <- df %>% number_results_per_panel(PDAST_1,PDAST_6,iv_singles,"S","I")
  print(4/60)
  df$PDAST_rpp_ics <- df %>% number_results_per_panel(PDAST_1,PDAST_6,iv_combos,"S","I")
  print(5/60)
  df$PDAST_rpp_ais <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_ivs,"S","I")
  print(6/60)
  df$PDAST_rpp_oss <- df %>% number_results_per_panel(PDAST_1,PDAST_6,oral_singles,"S","I")
  print(7/60)
  df$PDAST_rpp_ocs <- df %>% number_results_per_panel(PDAST_1,PDAST_6,oral_combos,"S","I")
  print(8/60)
  df$PDAST_rpp_aos <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_orals,"S","I")
  print(9/60)
  df$PDAST_rpp_acss <- df %>% number_results_per_panel(PDAST_1,PDAST_6,access_singles,"S","I")
  print(10/60)
  df$PDAST_rpp_accs <- df %>% number_results_per_panel(PDAST_1,PDAST_6,access_combos,"S","I")
  print(11/60)
  df$PDAST_rpp_aacs <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_access,"S","I")
  print(12/60)
  df$PDAST_rpp_wass <- df %>% number_results_per_panel(PDAST_1,PDAST_6,watch_singles,"S","I")
  print(13/60)
  df$PDAST_rpp_wacs <- df %>% number_results_per_panel(PDAST_1,PDAST_6,watch_combos,"S","I")
  print(14/60)
  df$PDAST_rpp_awas <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_watch,"S","I")
  print(15/60)
  df$PDAST_rpp_asr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_singles,"R","NT")
  print(16/60)
  df$PDAST_rpp_acr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_combos,"R","NT")
  print(17/60)
  df$PDAST_rpp_aar <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_abs,"R","NT")
  print(18/60)
  df$PDAST_rpp_isr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,iv_singles,"R","NT")
  print(19/60)
  df$PDAST_rpp_icr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,iv_combos,"R","NT")
  print(20/60)
  df$PDAST_rpp_air <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_ivs,"R","NT")
  print(21/60)
  df$PDAST_rpp_osr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,oral_singles,"R","NT")
  print(22/60)
  df$PDAST_rpp_ocr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,oral_combos,"R","NT")
  print(23/60)
  df$PDAST_rpp_aor <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_orals,"R","NT")
  print(24/60)
  df$PDAST_rpp_acsr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,access_singles,"R","NT")
  print(25/60)
  df$PDAST_rpp_accr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,access_combos,"R","NT")
  print(26/60)
  df$PDAST_rpp_aacr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_access,"R","NT")
  print(27/60)
  df$PDAST_rpp_wasr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,watch_singles,"R","NT")
  print(28/60)
  df$PDAST_rpp_wacr <- df %>% number_results_per_panel(PDAST_1,PDAST_6,watch_combos,"R","NT")
  print(29/60)
  df$PDAST_rpp_awar <- df %>% number_results_per_panel(PDAST_1,PDAST_6,all_watch,"R","NT")
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
  
  PDAST_rpp_ass <- df %>% create_df(PDAST_rpp_ass,"PDAST\nsingle S","PDAST","Single","S")
  PDAST_rpp_acs <- df %>% create_df(PDAST_rpp_acs,"PDAST\ncombo S","PDAST","Combo","S")
  PDAST_rpp_aas <- df %>% create_df(PDAST_rpp_aas,"PDAST\nall S","PDAST","All","S")
  PDAST_rpp_iss <- df %>% create_df(PDAST_rpp_iss,"PDAST\nsingle IV S","PDAST","Single IV","S")
  PDAST_rpp_ics <- df %>% create_df(PDAST_rpp_ics,"PDAST\ncombo IV S","PDAST","Combo IV","S")
  PDAST_rpp_ais <- df %>% create_df(PDAST_rpp_ais,"PDAST\nall IV S","PDAST","All IV","S")
  PDAST_rpp_oss <- df %>% create_df(PDAST_rpp_oss,"PDAST\nsingle oral S","PDAST","Single Oral","S")
  PDAST_rpp_ocs <- df %>% create_df(PDAST_rpp_ocs,"PDAST\ncombo oral S","PDAST","Combo Oral","S")
  PDAST_rpp_aos <- df %>% create_df(PDAST_rpp_aos,"PDAST\nall oral S","PDAST","All Oral","S")
  PDAST_rpp_acss <- df %>% create_df(PDAST_rpp_acss,"PDAST\nsingle Access S","PDAST","Single Access","S")
  PDAST_rpp_accs <- df %>% create_df(PDAST_rpp_accs,"PDAST\ncombo Access S","PDAST","Combo Access","S")
  PDAST_rpp_aacs <- df %>% create_df(PDAST_rpp_aacs,"PDAST\nall Access S","PDAST","All Access","S")
  PDAST_rpp_wass <- df %>% create_df(PDAST_rpp_wass,"PDAST\nsingle Watch S","PDAST","Single Watch","S")
  PDAST_rpp_wacs <- df %>% create_df(PDAST_rpp_wacs,"PDAST\ncombo watch S","PDAST","Combo Watch","S")
  PDAST_rpp_awas <- df %>% create_df(PDAST_rpp_awas,"PDAST\nall watch S","PDAST","All Watch","S")
  PDAST_rpp_asr <- df %>% create_df(PDAST_rpp_asr,"PDAST\nsingle R","PDAST","Single","R")
  PDAST_rpp_acr <- df %>% create_df(PDAST_rpp_acr,"PDAST\ncombo R","PDAST","Combo","R")
  PDAST_rpp_aar <- df %>% create_df(PDAST_rpp_aar,"PDAST\nall R","PDAST","All","R")
  PDAST_rpp_isr <- df %>% create_df(PDAST_rpp_isr,"PDAST\nsingle IV R","PDAST","Single IV","R")
  PDAST_rpp_icr <- df %>% create_df(PDAST_rpp_icr,"PDAST\ncombo IV R","PDAST","Combo IV","R")
  PDAST_rpp_air <- df %>% create_df(PDAST_rpp_air,"PDAST\nall IV R","PDAST","All IV","R")
  PDAST_rpp_osr <- df %>% create_df(PDAST_rpp_osr,"PDAST\nsingle oral R","PDAST","Single Oral","R")
  PDAST_rpp_ocr <- df %>% create_df(PDAST_rpp_ocr,"PDAST\ncombo oral R","PDAST","Combo Oral","R")
  PDAST_rpp_aor <- df %>% create_df(PDAST_rpp_aor,"PDAST\nall oral R","PDAST","All Oral","R")
  PDAST_rpp_acsr <- df %>% create_df(PDAST_rpp_acsr,"PDAST\nsingle Access R","PDAST","Single Access","R")
  PDAST_rpp_accr <- df %>% create_df(PDAST_rpp_accr,"PDAST\ncombo Access R","PDAST","Combo Access","R")
  PDAST_rpp_aacr <- df %>% create_df(PDAST_rpp_aacr,"PDAST\nall Access R","PDAST","All Access","R")
  PDAST_rpp_wasr <- df %>% create_df(PDAST_rpp_wasr,"PDAST\nsingle Watch R","PDAST","Single Watch","R")
  PDAST_rpp_wacr <- df %>% create_df(PDAST_rpp_wacr,"PDAST\ncombo Watch R","PDAST","Combo Watch","R")
  PDAST_rpp_awar <- df %>% create_df(PDAST_rpp_awar,"PDAST\nall Watch R","PDAST","All Watch","R")
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
  
  acs_df <- data.frame(rbind(PDAST_rpp_ass,PDAST_rpp_acs,PDAST_rpp_aas,PDAST_rpp_iss,PDAST_rpp_ics,
                             PDAST_rpp_ais,PDAST_rpp_oss,PDAST_rpp_ocs,PDAST_rpp_aos,PDAST_rpp_acss,
                             PDAST_rpp_accs,PDAST_rpp_aacs,PDAST_rpp_wass,PDAST_rpp_wacs,PDAST_rpp_awas,
                             PDAST_rpp_asr,PDAST_rpp_acr,PDAST_rpp_aar,PDAST_rpp_isr,PDAST_rpp_icr,
                             PDAST_rpp_air,PDAST_rpp_osr,PDAST_rpp_ocr,PDAST_rpp_aor,PDAST_rpp_acsr,
                             PDAST_rpp_accr,PDAST_rpp_aacr,PDAST_rpp_wasr,PDAST_rpp_wacr,PDAST_rpp_awar,
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
main_dotplotter <- function(df,pdast_1,standard_1,pdast_2,standard_2,
                            left_label,right_label,sens_addendum="") {
  
  plot_df <- df %>% filter(grepl(pdast_1,AWaRe_results) |
                             grepl(standard_1,AWaRe_results) |
                             grepl(pdast_2,AWaRe_results) |
                             grepl(standard_2,AWaRe_results) )
  
  plot_df$AWaRe_results <- factor(plot_df$AWaRe_results,levels = c(pdast_1,
                                                                   standard_1,
                                                                   pdast_2,
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

###Assembling antimicrobial plot dataframe
abs_df_assemble <- function(pd_df,standard_df) {
  
  minuser <- function(df,abx) {
    
    df %>% filter(ind==ab_name(abx)) %>% arrange(Approach) %>% select(values)
    
    if(nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
       abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
      
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
        abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2)
      
    } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
               abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="Standard"){
      
      -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1) -
          abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(2))
      
    } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==1 &
               abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% dplyr::slice(1) =="PDAST") {
      
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1)
      
    } else {
      
      -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% dplyr::slice(1))
      
    }
    
  }
  
  abs_df <- bind_rows(
    standard_df %>% data.frame() %>% mutate(Approach = "Standard"),
    pd_df %>% data.frame() %>% mutate(Approach = "PDAST")
  ) %>% mutate(ind = ab_name(ind))
  abs_diffs <- map_df(all_singles, function(abs) {
    abs_df %>% 
      minuser(abs) %>% 
      tibble() %>% 
      mutate(ind = ab_name(abs))
  }) %>% 
    mutate(better = if_else(values > 0, "PDAST", "Standard"),
           values2 = abs(values)) %>% 
    left_join(
      abs_df %>% filter(Approach == "PDAST") %>% select(values, ind) %>% rename(PDAST = values), 
      by = "ind"
    ) %>% 
    left_join(
      abs_df %>% filter(Approach == "Standard") %>% select(values, ind) %>% rename(Standard = values), 
      by = "ind"
    ) %>% 
    mutate(
      Standard = if_else(is.na(Standard), 0, Standard),
      values = if_else(better == "PDAST", PDAST + 200, Standard + 200)
    )
  abs_diffs <<- abs_diffs
  standard_levels <- abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) 
  pdast_levels <- abs_df %>% filter(Approach == "PDAST") %>% filter(!ind %in% standard_levels) %>% arrange(values) %>% pull(ind)
  abs_df <- abs_df %>% 
    anti_join(
      abs_df %>% filter(Approach == "Standard") %>% select(ind), 
      by = "ind"
    ) %>% 
    mutate(Approach = "Standard", values = 0) %>% 
    bind_rows(abs_df) %>% 
    mutate(
      Approach = factor(Approach, levels = c("PDAST", "Standard")),
      ind = factor(ind, levels = c(pdast_levels,standard_levels)),
      aware = if_else(ind %in% ab_name(access_singles), "Access", "Watch")
    )
  axiscols <- if_else(
    abs_df %>% filter(Approach == "Standard") %>% arrange(values) %>% pull(ind) %in% ab_name(access_singles),
    "seagreen", "darkorange"
  )
  
  axiscols <<- axiscols
  
  abs_df
  
}

###Cleveland dot plot of AST results
cleveland_ab_plot <- function(ab_df,resulttype,addendum="") {
  
  results_by_ab <- ggplot(ab_df,aes(x=ind,y=values))+
    geom_line(aes(group=ind),alpha=0.5)+
    geom_point(aes(color=Approach),size=4) +
    coord_flip() +
    scale_color_manual(values=c("#00BFC4","#F8766D"))+
    geom_text(data = abs_diffs, aes(color = better, 
                                    label = as.character(glue("+{values2}"))),
              size = 3,hjust=0.5) +
    ggtitle(glue("Total number of {resulttype} AST results by antimicrobial agent\n{addendum}"))+
    xlab("") +
    ylab(glue("Total number of {resulttype} results")) +
    theme_minimal() +
    theme(axis.text.y = element_text(
      colour = axiscols))
  
  ggsave(glue("{resulttype}_ab_cleveplot.pdf"), plot = results_by_ab, device = "pdf", width = 10, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  results_by_ab
  
}

###Counting antimicrobial-result matches per paneL
totaller <- function(panel) {
  abx_graph %>% filter(Panel==panel) %>% summarise(Total=sum(n)) %>% unlist()
}

###Counting Access antimicrobial-result matches per paneL
awaretotaller <- function(panel) {
  awaresum %>% filter(Panel==panel & AWaRe=="Access" & Result=="S") %>%
    select(n_results) %>% unlist()
}

###Assembling dataframe for treatment recommendations bar plot
ab_rec_df <- function(df,pd_rec,quote_pdrec,pd_res,st_res,st_rec) {
  
  pd_rec <- enquo(pd_rec)
  pd_res <- enquo(pd_res)
  st_res <- enquo(st_res)
  
  urines_abx <- urines_abx %>% filter(!!pd_rec %in% all_singles)
  
  pdast_graph <- urines_abx %>% count(!!pd_rec,!!pd_res,!!st_res) %>% 
    mutate(AWaRe = case_when(!!pd_rec %in% access_singles ~ "Access",
                             TRUE~"Watch")) %>% 
    mutate(!!pd_rec:=ab_name(!!pd_rec),
           Panel="PDAST",
           Result=!!pd_res) %>% 
    rename(abx_name=quote_pdrec)
  standard_graph <- pdast_graph %>% mutate(Panel="Standard",
                                           Result=!!st_res)
  ab_gr_df <- data.frame(rbind(pdast_graph,standard_graph))
  ab_gr_df <- ab_gr_df %>% mutate(abx_name = case_when(Panel=="Standard" ~ st_rec,
                                                       TRUE ~ abx_name))
  
  axiscols <- ifelse(ab_gr_df %>% group_by(abx_name) %>% 
                       mutate(n=sum(n)) %>% 
                       arrange(n) %>% ungroup() %>% 
                       distinct(abx_name) %>% unlist() %in% ab_name(access_abs),"seagreen",
                     "darkorange")
  
  axiscols <<- axiscols
  
  ab_gr_df$abx_name <- factor(ab_gr_df$abx_name,
                              levels=ab_gr_df %>% group_by(abx_name) %>% 
                                mutate(n=sum(n)) %>% 
                                arrange(n) %>% ungroup() %>% 
                                distinct(abx_name) %>% unlist())
  max_count <- ceiling((ab_gr_df %>% group_by(abx_name,Panel) %>% summarise(sumn=sum(n)) %>% 
                          ungroup() %>% summarise(max=max(sumn)) %>% unlist() %>% 
                          as.numeric() /25) * 25)
  
  max_count <<- max_count
  
  ab_gr_df
  
}

###Bar plot of treatment recommendations - results
ab_result_graph <- function(df,line,route,title_height,title_width) {
  
  ab_plot <- ggplot(df, aes(x = abx_name, y = if_else(Panel == "Standard", -n, n),
                            fill = Result,color=Result)) +
    geom_bar(stat = "identity", width=0.85,position = "stack") +
    scale_y_continuous(
      limits = c(-ceiling(max_count / 100) * 100, ceiling(max_count / 100) * 100), 
      breaks = seq(-ceiling(max_count / 500) * 500, ceiling(max_count / 500) * 500,by=500), 
      labels = abs(seq(-ceiling(max_count / 500) * 500, ceiling(max_count / 500) * 500,by=500)) 
    ) +
    coord_flip() +
    labs(y = "Number of results", x = "Antimicrobial agent recommended",
         fill = "Result") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
    geom_hline(yintercept = 0,linetype="solid",color="black") +
    ggtitle(glue("Comparison of subsequent urine susceptibility results for antibiotics\nrecommended as {line}-line {route} treatment")) +
    geom_text(x=title_height,y=-title_width,label=glue("Standard {line}-line {route} agent"),color="#3C3C3C",size=4) +
    geom_text(x=title_height,y=title_width,label=glue("Personalised {line}-line {route} agent"),color="#3C3C3C",size=4) +
    theme(plot.title = element_text(size = 16, margin = margin(b = 20)),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y = element_text(
            colour = axiscols))
  
  
  
  ggsave(glue("uf_{line}-line_{route}_recommended.pdf"), plot = abx_prescribed, device = "pdf", width = 10, height = 4,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  ab_plot
  
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
util_ref <- read_csv("ur_util_pre_recs.csv")
ur_util <- ur_util %>% semi_join(util_ref,by="micro_specimen_id")
util_probs_df <- util_probs_df %>% semi_join(util_ref,by="micro_specimen_id")

##Discrete choice experiment edge case

###Simulate data with universal first and last options
resvec <- c(
  6,1,2,9,10,7,11,12,3,13,8,4,5
) %>% as.character()
results[2,10:22] <- as.list(resvec)
results <- results %>%
  dplyr::slice(rep(2, 1001))
results$`Respondent ID` <- seq(1,nrow(results))
results[1,] <- NA

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
score_rownames <- rownames(scores)

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
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle("The effect of different antimicrobial drug properties on clinician prescribing\npreference in the UTI scenario discrete choice experiment (edge case)")+
  geom_vline(xintercept = 0,colour="grey40")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "None")

ggsave(glue("ORplot_edgecase.pdf"), plot = ORplot, device = "pdf", width = 10, height = 8,
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

###Calculate overall utility score
uprobdf <- util_probs_df %>% calculate_utilities(R_weight = R_wt_value)

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
uprobdf <- uprobdf %>% filter(Antimicrobial %in% abx_in_train)

uplot <- uprobdf %>% utility_plot(Rx_utility,"Treatment", " (single agent in edge case)")

ggsave(glue("utility_Treatment_ (single agent in edge case).pdf"), plot = uplot, device = "pdf", width = 10, height = 8,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

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
reversed_combomap <- list()
for (i in 1:length(combined_antimicrobial_map)) {
  reversed_combomap[i] <- names(combined_antimicrobial_map)[i]
  names(reversed_combomap)[i] <- combined_antimicrobial_map[i]
}
reversed_combomap
replace_values <- function(column, map) {
  flipped_map <- setNames(names(map), map)
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(flipped_map)) flipped_map[[x]] else x)
}

###Illness severity sensitivity analysis
weightseq <- seq(0,3)
iv_perclist <- c()
po_perclist <- c()
overall_perclist <- c()
ivac_list <- c()
poac_list <- c()
ovac_list <- c()
ovor_list <- c()
oviv_list <- c()
pdrx1_list <- data.frame(PDRx_1=ab_singles)

for(weight in seq_along(weightseq)) {
  
  res_sens_analysis(ur_util,util_probs_df,MEWS_variable=weightseq[weight],R_value=1)
  
  iv_perclist <- c(iv_perclist,iv_perc)
  po_perclist <- c(po_perclist,po_perc)
  overall_perclist <- c(overall_perclist,overall_perc)
  ivac_list <- c(ivac_list,iv_s_access)
  poac_list <- c(poac_list,po_s_access)
  ovac_list <- c(ovac_list,overall_s_access)
  ovor_list <- c(ovor_list,overall_s_oral)
  oviv_list <- c(oviv_list,overall_s_iv)
  pdrx1_list <- pdrx1_list %>% left_join(abrx1_df,by="PDRx_1")
  colnames(pdrx1_list)[weight+1] <- weightseq[weight]
  
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

iv_xg_plot_df <- data.frame(rbind(
  iv_perclist,ivac_list
))
po_xg_plot_df <- data.frame(rbind(
  po_perclist,poac_list
))
overall_xg_plot_df <- data.frame(rbind(
  overall_perclist,ovac_list,ovor_list,oviv_list
))

pdrx1_list[is.na(pdrx1_list)] <- 0
pdrx1_df <- pdrx1_list %>% filter(rowSums(select(.,2:ncol(pdrx1_list)))!=0)
pdrx1_df <- melt(pdrx1_df)
colnames(pdrx1_df) <- c("Antimicrobial","Illness severity","Percentage of first-line recommendations")
pdrx1_df <- pdrx1_df %>% mutate(Antimicrobial = ab_name(Antimicrobial))
pdrx1_df <- pdrx1_df %>% mutate(`Percentage of first-line recommendations`=
                                  (`Percentage of first-line recommendations`/nrow(ur_util))*100)

write_csv(pdrx1_df,"abplot_df_edgecase.csv")

abplot <- ggplot(pdrx1_df,aes(x=`Illness severity`,y=`Percentage of first-line recommendations`,group=Antimicrobial,colour=Antimicrobial))+
  geom_line()+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("First-line automated recommendations according to illness severity")+
  scale_x_discrete(breaks = NULL)+
  ylim(0,100)

ggsave(glue("illness_abplot_edgecase.pdf"), plot = abplot, device = "pdf", width = 10, height = 6,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

print(abplot)

write_csv(iv_xg_plot_df,"iv_xg_plot_df_edgecase.csv")
write_csv(po_xg_plot_df,"po_xg_plot_df_edgecase.csv")
write_csv(overall_xg_plot_df,"overall_xg_plot_df_edgecase.csv")

overall_xg_plot_iv <- overall_xg_plot_df %>% filter(!grepl("(access|oral)",Metric,ignore.case=T))
overall_xg_plot_iv %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                            agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                            suffix="according to illness severity (proportion of IV-administrable agents, edge case)")
overall_xg_plot_oral <- overall_xg_plot_df %>% filter(!grepl("(access|iv)",Metric,ignore.case=T))
overall_xg_plot_oral %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                              agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                              suffix="according to illness severity (proportion of orally-administrable agents, edge case)")
overall_xg_plot_access <- overall_xg_plot_df %>% filter(!grepl("(oral|iv)",Metric,ignore.case=T))
overall_xg_plot_access %>% susc_plotter_overall(ur_util,"overall ", measure="MEWS",variable="Illness severity",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam",
                                                suffix="according to illness severity (proportion of WHO Access agents, edge case)")

##Change antibiotics to make UTI-specificity less correlated with other variables

results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")

###Simulate data with universal first and last options
resvec <- c(
  6,1,2,9,10,7,11,12,3,13,8,4,5
) %>% as.character()
results[2,10:22] <- as.list(resvec)
results <- results %>%
  dplyr::slice(rep(2, 1001))
results$`Respondent ID` <- seq(1,nrow(results))
results[1,] <- NA

###Re-factorising outcome variables on abx dataframes after read-in
train_abx <- train_abx %>% factorise()
test_abx <- test_abx %>% factorise()

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
Endoxolol <- c("Access","Low","Low","No","Yes","No","Low")
Olanzasys <- c("Watch","High","Low","Yes","No","Yes","Low")
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
score_rownames <- rownames(scores)

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
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle("The effect of different antimicrobial drug properties on clinician prescribing\npreference in the UTI scenario discrete choice experiment (edge case 2)")+
  geom_vline(xintercept = 0,colour="grey40")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "None")

ggsave(glue("ORplot_edgecase2.pdf"), plot = ORplot, device = "pdf", width = 10, height = 8,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")
print(ORplot)
