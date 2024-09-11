#RECOMMENDATION FUNCTION

##Functions

###Factorise training and testing datasets
factorise <- function(df) {
  df %>% mutate(CDI = factor(CDI),
                 overall_tox = factor(overall_tox),
                 sepsis_ae=factor(sepsis_ae))
}

###Prioritisation by treatment utility
util_mk1 = function(df,spec_id,panel_size) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(desc(Rx_utility)) %>% select(Antimicrobial,Rx_utility) %>% 
    mutate(Rx_utility = round(Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "Rx_utility")
  
}

###Prioritisation by AST utility
util_mk2 = function(df,spec_id,panel_size) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(desc(AST_utility)) %>% select(Antimicrobial,AST_utility) %>% 
    mutate(AST_utility = round(AST_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`AST Utility` = "AST_utility")
  
}

###Prioritisation by urosepsis utility
util_mk3 = function(df,spec_id,panel_size) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(desc(Urosepsis_Rx_utility)) %>% select(Antimicrobial,Urosepsis_Rx_utility) %>% 
    mutate(Urosepsis_Rx_utility = round(Urosepsis_Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Urosepsis Rx Utility` = "Urosepsis_Rx_utility")
  
}

###Prioritisation by outpatient utility
util_mk4 = function(df,spec_id,panel_size) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(desc(Outpatient_Rx_utility)) %>% select(Antimicrobial,Outpatient_Rx_utility) %>% 
    mutate(Outpatient_Rx_utility = round(Outpatient_Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Oupatient Rx Utility` = "Outpatient_Rx_utility")
  
}

###Assigning treatment recommendations
assign_PDRx <- function(df,probab_df,method_used) {
  
  
  
  test_recs <-  data.frame(matrix(nrow=length(all_abs),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk1(spec_id = df$micro_specimen_id[i], panel_size = length(all_abs)) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(all_abs)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning AST recommendations
assign_PDAST <- function(df,probab_df,method_used) {
  
  test_recs <-  data.frame(matrix(nrow=length(all_abs),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk2(spec_id = df$micro_specimen_id[i], panel_size = length(all_abs)) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(all_abs)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning urosepsis treatment recommendations
assign_urosepsis <- function(df,probab_df,method_used) {
  
  
  
  test_recs <-  data.frame(matrix(nrow=length(iv_abs),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk3(spec_id = df$micro_specimen_id[i], panel_size = length(iv_abs)) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(iv_abs)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning outpatient treatment recommendations
assign_outpatient <- function(df,probab_df,method_used) {
  
  
  test_recs <-  data.frame(matrix(nrow=length(oral_abs),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk4(spec_id = df$micro_specimen_id[i], panel_size = length(oral_abs)) %>% 
      select(1)
    
    test_recs <- cbind(test_recs,rec)
    
    print(glue("{round((i/nrow(df)) * 100,0)}%"))
    
  }
  
  test_recs <- data.frame(t(test_recs))
  test_recs <- data.frame(cbind(df$micro_specimen_id,test_recs))
  testrec_cols <- c("micro_specimen_id")
  
  for(i in 1:length(oral_abs)) {
    
    testrec_cols[i+1] <- paste0(method_used,i)
    
  }
  
  colnames(test_recs) <- testrec_cols
  
  df %>% 
    left_join(test_recs,by="micro_specimen_id") 
  
}

###Assigning standard recommendations
assign_standard <- function(df,probab_df,micro_df,method_used) {
  
  ab_vector <- probab_df %>% mutate(Antimicrobial=as.ab(Antimicrobial)) %>% 
    distinct(Antimicrobial) %>% pull(Antimicrobial)
  standard_panel <- micro_df %>% filter(!is.na(org_name) & test_name == "URINE CULTURE") %>%
    count(ab_name) %>% arrange(desc(n)) %>% mutate(ab_name=as.ab(ab_name)) %>% pull(ab_name) %>% 
    intersect(ab_vector)
  print(standard_panel)
  standard_columns <- paste0(method_used, seq_along(standard_panel))
  df <- df %>%
    bind_cols(setNames(as.data.frame(matrix(NA, nrow = nrow(df), ncol = length(standard_columns))), standard_columns))
  for (i in seq_along(standard_panel)) {
    df[[standard_columns[i]]] <- standard_panel[i]
  }
  
  df
  
}

###Utility score calculator
calculate_utilities <- function(df,formulary_list=NULL) {
  
  df <- df %>% mutate(overall_util = exp(util_CDI + util_tox +
                                     util_uti + util_access +
                                     util_oral + util_iv +
                                     util_reserve + util_highcost),
                S_utility = S*overall_util,
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
                Rx_utility = (overall_util * S) -
                  (R*prob_sepsisae))
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv == 0 ~ min(df$Rx_utility)-0.01, TRUE ~ Rx_utility
      ),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01, TRUE ~ Rx_utility
      )
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
    
    df <- df %>% filter(AST_utility != min(AST_utility))
    
  }
  
  thisplot <- ggplot(df %>% mutate(Antimicrobial = 
                                     factor(Antimicrobial,
                                            levels = df %>% group_by(Antimicrobial) %>% 
                                              summarise(Median_util=median(!!variable)) %>% 
                                              arrange(Median_util) %>% select(Antimicrobial) %>% unlist())), 
                     aes(x=!!variable,y=Antimicrobial,fill=Antimicrobial)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position = "None",axis.text.y = element_text(
      colour = axiscols))+
    xlab(glue("{application} utility{modification}"))+
    theme()
  
  print(thisplot)
  
}

###Sensitivity analysis probability distribution check
dens_check <- function(df,abx_agent,result,alpha,beta) {
  
  result <- enquo(result)
  
  test_dist <- tibble(prob=rbeta(nrow(ur_util),alpha,beta),dist="test")
  actual_dist <- tibble(
    prob=df %>% filter(Antimicrobial==abx_agent) %>% pull(!!result),
    dist="actual")
  dist_df <- tibble(rbind(test_dist,actual_dist))
  
  thisplot <- ggplot(dist_df,
                     aes(x=prob,color=dist)) + 
    geom_density()+
    ggtitle(abx_agent)
  
  print(thisplot)
  
}

###Replacing probability distribution with simulated distribution
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

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_pdast <- function(df) {
  
  all_si <- c()
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(PDAST_1:PDAST_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_pdrx <- function(df) {
  
  all_si <- c()
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(PDRx_1:PDRx_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(PDRx_1:PDRx_6) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Counting the number of S results per antimicrobial provided by the personalised approach
number_abs_standard <- function(df) {
  
  all_si <- c()
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="S") %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(STANDARD_1,STANDARD_2,STANDARD_3,
                                                      STANDARD_7,STANDARD_8,STANDARD_11) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. =="I") %>% rownames()
    
    
    ac_si <- all_s %>% append(all_i)
    
    all_si <- all_si %>% append(ac_si)
    
  }
  
  all_si %>% table() %>% stack()
  
}

###Comparison of differences between number of results per antimicrobial with each approach
minuser <- function(df,abx) {
  
  df %>% filter(ind==ab_name(abx)) %>% arrange(Approach) %>% select(values)
  
  if(nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
     abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1) -
      abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(2)
    
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==2 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="Standard"){
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1) -
        abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(2))
    
  } else if (nrow(abs_df %>% filter(ind==ab_name(abx)) %>% select(1)) ==1 &
             abs_df %>% filter(ind==ab_name(abx)) %>% select(Approach) %>% slice(1) =="PDAST") {
    
    abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1)
    
  } else {
    
    -(abs_df %>% filter(ind==ab_name(abx)) %>% select(1) %>% slice(1))
    
  }
  
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

###Factorising outcome variables on abx dataframes
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
ggplot(scores,aes(x=OR_dif,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Odds ratio for drug selection") +
  ggtitle("The effect of different antimicrobial drug properties\non clinician prescribing preference in UTI scenario")+
  scale_x_continuous(labels = function(x) x+1)+
  geom_vline(xintercept = 0)

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
curr_service_columns <- grep("^curr_service_", names(train_abx), value = TRUE)
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
  ) %>% left_join(Rval_key)

###Calculate overall utility score
formulary_agents <- c() ####Populate with formulary agent full names
util_probs_df <- util_probs_df %>% calculate_utilities(
  formulary_list = formulary_agents)

###Filter out combination predictions not present in training dataset
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
util_probs_df <- util_probs_df %>% filter(Antimicrobial %in% abx_in_train)

##Utility analysis

###Utility data visualisation
util_probs_df %>% group_by(Antimicrobial) %>% 
  summarise(Median_util=median(Rx_utility)) %>% 
  arrange(desc(Median_util))
util_probs_df %>% utility_plot(Rx_utility,"Treatment")
util_probs_df %>% utility_plot(Urosepsis_Rx_utility,"Intravenous treatment")
util_probs_df %>% utility_plot(Outpatient_Rx_utility,"Oral treatment")
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

###Specific probability density check
util_probs_df %>% dens_check("Nitrofurantoin",R,12,3) ####Add full name of antimicrobial of interest

###Replace probability distribution with simulated distribution
sens_df <- util_probs_df %>% dist_replace(ur_util,"Ampicillin_Vancomycin","R","S",3,12) %>%
  calculate_utilities() ####Add atimicrobial agent of interest
sens_df %>% utility_plot(Rx_utility,"Intravenous treatment"," (Vanc R down)")
sens_df %>% utility_plot(AST_utility,"Test utility")

###Antimicrobial formulary recommendations
senskey <- ur_util %>% select(micro_specimen_id,AMP:VAN,AMP_SAM:NIT_VAN)
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
form_recs <- util_probs_df %>%
  left_join(senskey) %>% 
  rowwise() %>%
  mutate(
    Confirmed_S = ifelse(combined_antimicrobial_map[[Antimicrobial]] %in% names(.) && 
                           get(combined_antimicrobial_map[[Antimicrobial]]) == "S", TRUE, FALSE),
    Confirmed_R = ifelse(combined_antimicrobial_map[[Antimicrobial]] %in% names(.) && 
                           get(combined_antimicrobial_map[[Antimicrobial]]) == "R", TRUE, FALSE)
  ) %>%
  ungroup() %>% mutate(
    Urosepsis_Formutil = case_when(util_iv==0 ~ 0, TRUE ~ (Confirmed_S*overall_util) - (Confirmed_R*sepsis_ae)),
    Outpatient_Formutil = case_when(util_oral==0~0, TRUE ~ (Confirmed_S*overall_util) - (Confirmed_R*sepsis_ae)),
  ) %>% group_by(Antimicrobial) %>% 
  summarise(Formulary_util_urosepsis=mean(Urosepsis_Formutil),
            Formulary_util_outpatient=mean(Outpatient_Formutil)
  ) %>% ungroup()

form_recs_urosepsis <- form_recs %>% 
  arrange(desc(Formulary_util_urosepsis)) %>% 
  slice(1:length(iv_abs)) %>% select(Antimicrobial) %>% unlist()

form_recs_outpatient <- form_recs %>% 
  arrange(desc(Formulary_util_outpatient)) %>% 
  slice(1:length(oral_abs)) %>% select(Antimicrobial) %>% unlist()

for (i in 1:length(form_recs_urosepsis)) {
  ur_util <- ur_util %>%
    mutate(!!paste0("PDFo_Urosepsis_", i) := form_recs_urosepsis[[i]])
}

for (i in 1:length(form_recs_outpatient)) {
  ur_util <- ur_util %>%
    mutate(!!paste0("PDFo_Outpatient_", i) := form_recs_outpatient[[i]])
}


###Individual treatment recommendations
all_abs <- all_combos
replace_values <- function(column, map) {
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(map)) map[[x]] else x)
}

ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_") %>% 
  mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))

###Urosepsis treatment recommendations
ur_util <- ur_util %>% assign_urosepsis(util_probs_df,"Urosepsis_") %>% 
  mutate(across(starts_with("Urosepsis_"), ~ replace_values(., combined_antimicrobial_map)))

###Outpatient treatment recommendations
ur_util <- ur_util %>% assign_outpatient(util_probs_df,"Outpatient_") %>% 
  mutate(across(starts_with("Outpatient_"), ~ replace_values(., combined_antimicrobial_map)))

###Individual AST recommendations
ur_util <- ur_util %>% assign_PDAST(util_probs_df,"PDAST_") %>%
  mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))

###Standard panel treatment & AST recommendations
micro_raw <- read_csv("microbiologyevents.csv")
ur_util <- ur_util %>% assign_standard(util_probs_df,micro_raw,"STANDARD_")

##Treatment recommendation analysis





