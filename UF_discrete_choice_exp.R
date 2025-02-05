#DISCRETE CHOICE EXPERIMENT

set.seed(123)

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
calculate_utilities <- function(df,formulary_list=c(),R_weight=1,MEWS=1) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + 
                        util_reserve + util_highcost 
                      + util_tox + util_CDI + util_iv*MEWS^2*R_weight,
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
rownames(scores) <- score_rownames

###Export unstandardised weights for reference
unstan_vals <- scores %>% select(Value) %>% 
  mutate(Name = c("Wc", "Wt", "Wo", "Wu", "Wi", "Wh", "Wa", "Wr"))
write_csv(unstan_vals, "unstand_weights.csv")

###Visualise coefficients for antimicrobial characteristics
scores$Coefficient <- factor(scores$Coefficient, levels=
                         scores %>% arrange(Value) %>% 
                         select(Coefficient) %>% unlist())

ORplot <- ggplot(scores,aes(x=Value,y=Coefficient,fill=colour)) +
  geom_col() +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle("The effect of different antimicrobial drug properties on clinician prescribing\npreference in the UTI scenario discrete choice experiment")+
  geom_vline(xintercept = 0,colour="grey40")+
  geom_errorbar(aes(xmin = Q5, xmax = Q95), width = 0.1,colour="grey40")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "None")

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

##Utility score preprocessing

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
cost_list %>% dplyr::slice(1:13) %>% summarise(
  MED = Median(orig_cost),
  iq1 = quantile(orig_cost)[2],
  iq3 = quantile(orig_cost)[4]
)

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

write_csv(util_probs_df,"dce_util_probs_df.csv")


