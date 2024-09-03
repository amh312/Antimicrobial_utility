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
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "AST_utility")
  
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
  
  df %>% mutate(overall_util = exp(util_CDI + util_tox +
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
                AST_utility = S_utility+
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
                  Formulary_utility,
                Rx_utility = (overall_util * S) -
                  (R*prob_sepsisae))
  
}

###Utility data visualisation
utility_plot <- function(df, variable,utility) {
  
  variable <- enquo(variable)
  
  axiscols <- if_else(
    df %>% pull(Antimicrobial) %in% formulary_agents,
    "seagreen", "black")
  
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
    xlab(utility)+
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
util_probs_df <- util_probs_df %>% 
  mutate(abx_name_Ampicillin.sulbactam=case_when(
    Antimicrobial=="Ampicillin-sulbactam" ~ 1,
    TRUE ~ 0),
    abx_name_Cefazolin=case_when(
      Antimicrobial=="Cefazolin" ~ 1,
      TRUE ~ 0),
    abx_name_Cefepime=case_when(
      Antimicrobial=="Cefepime" ~ 1,
      TRUE ~ 0),
    abx_name_Ceftazidime=case_when(
      Antimicrobial=="Ceftazidime" ~ 1,
      TRUE ~ 0),
    abx_name_Ceftriaxone=case_when(
      Antimicrobial=="Ceftriaxone" ~ 1,
      TRUE ~ 0),
    abx_name_Ciprofloxacin=case_when(
      Antimicrobial=="Ciprofloxacin" ~ 1,
      TRUE ~ 0),
    abx_name_Gentamicin=case_when(
      Antimicrobial=="Gentamicin" ~ 1,
      TRUE ~ 0),
    abx_name_Meropenem=case_when(
      Antimicrobial=="Meropenem" ~ 1,
      TRUE ~ 0),
    abx_name_Nitrofurantoin=case_when(
      Antimicrobial=="Nitrofurantoin" ~ 1,
      TRUE ~ 0),
    abx_name_Piperacillin.tazobactam=case_when(
      Antimicrobial=="Piperacillin-tazobactam" ~ 1,
      TRUE ~ 0),
    abx_name_Trimethoprim.sulfamethoxazole=case_when(
      Antimicrobial=="Trimethoprim-sulfamethoxazole" ~ 1,
      TRUE ~ 0),
    abx_name_Ampicillin=case_when(
      Antimicrobial=="Ampicillin" ~ 1,
      TRUE ~ 0))

##CDI prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
cdi_fit <- log_reg_spec %>%
  fit(CDI ~ pCDI+pHADM+age65+pCKD+pDIAB+pLIVER+pCARD+pCVA+pCA+MALE+
        abx_name_Ampicillin.sulbactam+abx_name_Cefazolin+
        abx_name_Cefepime+abx_name_Ceftazidime+
        abx_name_Ceftriaxone+abx_name_Ciprofloxacin+
        abx_name_Gentamicin+abx_name_Meropenem+
        abx_name_Nitrofurantoin+abx_name_Piperacillin.tazobactam+
        abx_name_Trimethoprim.sulfamethoxazole+
        curr_service_CSURG+curr_service_ENT+curr_service_GU+
        curr_service_GYN+curr_service_MED+curr_service_NMED+
        curr_service_NSURG+curr_service_OBS+curr_service_OMED+
        curr_service_ORTHO+curr_service_PSURG+curr_service_PSYCH+
        curr_service_SURG+curr_service_TRAUM+curr_service_TSURG+
        curr_service_VSURG+pICU+pSEPSIS,
      data = train_abx)

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
tox_fit <- log_reg_spec %>%
  fit(overall_tox ~ prAKI+pHADM+age65+pCKD+pDIAB+pLIVER+pCARD+pCVA+pCA+MALE+
        abx_name_Ampicillin.sulbactam+abx_name_Cefazolin+
        abx_name_Cefepime+abx_name_Ceftazidime+
        abx_name_Ceftriaxone+abx_name_Ciprofloxacin+
        abx_name_Gentamicin+abx_name_Meropenem+
        abx_name_Nitrofurantoin+abx_name_Piperacillin.tazobactam+
        abx_name_Trimethoprim.sulfamethoxazole +
        curr_service_CSURG+curr_service_ENT+curr_service_GU+
        curr_service_GYN+curr_service_MED+curr_service_NMED+
        curr_service_NSURG+curr_service_OBS+curr_service_OMED+
        curr_service_ORTHO+curr_service_PSURG+curr_service_PSYCH+
        curr_service_SURG+curr_service_TRAUM+curr_service_TSURG+
        curr_service_VSURG+pICU+pSEPSIS,
      data = train_abx)

underlying_tox <- extract_fit_engine(tox_fit)
coef(underlying_tox) %>% round(2)

###Evaluate model on test data
tox_test_probs <- predict(underlying_tox, test_abx,type="response")

###Attach to probability prediction dataframe
tox_util_key <- ur_util %>% select(micro_specimen_id,overall_tox)
util_probs_df <- util_probs_df %>% 
  left_join(tox_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

##Sepsis adverse outcomes prediction model

###Initialise model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

###Fit model
sepsis_fit <- log_reg_spec %>%
  fit(sepsis_ae ~ pCDI+pHADM+age65+pCKD+pDIAB+pLIVER+pCARD+pCVA+pCA+MALE+
        abx_name_Ampicillin.sulbactam+abx_name_Cefazolin+
        abx_name_Cefepime+abx_name_Ceftazidime+
        abx_name_Ceftriaxone+abx_name_Ciprofloxacin+
        abx_name_Gentamicin+abx_name_Meropenem+
        abx_name_Nitrofurantoin+abx_name_Piperacillin.tazobactam+
        abx_name_Trimethoprim.sulfamethoxazole+
        curr_service_CSURG+curr_service_ENT+curr_service_GU+
        curr_service_GYN+curr_service_MED+curr_service_NMED+
        curr_service_NSURG+curr_service_OBS+curr_service_OMED+
        curr_service_ORTHO+curr_service_PSURG+curr_service_PSYCH+
        curr_service_SURG+curr_service_TRAUM+curr_service_TSURG+
        curr_service_VSURG+pICU+pSEPSIS,
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
access_value <- scores[rownames(scores)=="Access",] %>% 
  select(Value) %>% unlist()

###Oral option utilituy
oral_abs <- c("AMP","SAM","CIP",
              "GEN","SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
  select(Value) %>% unlist()

###IV option utility
iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-")
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
  ) %>% left_join(Rval_key)

###Calculate overall utility score
formulary_agents <- c() ####Populate with formulary agent full names
util_probs_df <- util_probs_df %>% calculate_utilities(
  formulary_list = formulary_agents)

##Utility analysis

###Utility data visualisation
util_probs_df %>% group_by(Antimicrobial) %>% 
  summarise(Median_util=median(Rx_utility)) %>% 
  arrange(desc(Median_util))
util_probs_df %>% utility_plot(Rx_utility,"Treatment utility")
util_probs_df %>% utility_plot(AST_utility,"Test utility")

##Sensitivity analysis

###Probability density check across all antimicrobials
all_abs <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN")
long_allabs <- all_abs %>% ab_name() %>% str_replace("/","-")

for (i in seq_along(long_allabs)) {
  
  util_probs_df %>% dens_check(long_allabs[i],S,1,4) 
  
}

###Specific probability density check
util_probs_df %>% dens_check("Nitrofurantoin",R,12,3) ####Add full name of antimicrobial of interest

###Replace probability distribution with simulated distribution
nitro_sens_df <- util_probs_df %>% dist_replace(ur_util,"Nitrofurantoin","R","S",12,3) %>%
  calculate_utilities() ####Add atimicrobial agent of interest
nitro_sens_df %>% utility_plot(Rx_utility,"Treatment utility (Nitro R increased)")
nitro_sens_df %>% utility_plot(AST_utility,"Test utility (Nitro R increased)")

##Recommendations

###Antimicrobial formulary recommendations

form_recs <- util_probs_df %>% group_by(Antimicrobial) %>% 
  summarise(Mean_util=mean(S_utility)) %>% 
  arrange(desc(Mean_util)) %>% select(Antimicrobial) %>% unlist()

for (i in 1:length(form_recs)) {
  ur_util <- ur_util %>%
    mutate(!!paste0("PDFo_", i) := form_recs[[i]])
}

###Individual treatment recommendations
all_abs <- c("AMP","SAM","CZO",
             "GEN","SXT","NIT",
             "TZP","CRO","CAZ",
             "FEP","MEM","CIP","VAN")

ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_") %>% 
  mutate(across(starts_with("PDRx_"), as.ab))

###Individual AST recommendations
ur_util <- ur_util %>% assign_PDAST(util_probs_df,"PDAST_") %>%
  mutate(across(starts_with("PDAST_"), as.ab))

###Standard panel treatment & AST recommendations
micro_raw <- read_csv("microbiologyevents.csv")
ur_util <- ur_util %>% assign_standard(util_probs_df,micro_raw,"STANDARD_")