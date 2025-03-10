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
  formulary_agents <- c()
  
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

###Separated-specialty sensitivity analysis
spec_sens_analysis <- function(df,results_df,specialty) {
  
  ###Engineer scores dataframe
  rankings <- results_df
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
  
  ###Calculate overall utility score
  df <- df %>% calculate_utilities(R_weight = 1)
  
  ##Utility analysis
  
  ###Utility data visualisation
  df %>% group_by(Antimicrobial) %>% 
    summarise(Median_util=median(Rx_utility)) %>% 
    arrange(desc(Median_util))
  formulary_agents <- c()
  
  df
  
}

##Read-in and filter
results <- read_csv("labelled_ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
util_probs_df_2 <- read_csv("probs_df_overall.csv")
ref_util<- read_csv("util_probs_df_final.csv")

acuities <- ref_util %>% select(micro_specimen_id,acuity) %>% distinct(
  micro_specimen_id,.keep_all = T
)

util_probs_df_2 <- util_probs_df_2 %>% semi_join(ref_util,by="micro_specimen_id") %>% 
  left_join(acuities,by="micro_specimen_id")

###Read-in and filter questionnaire dataframes

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


