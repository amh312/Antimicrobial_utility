#DISCRETE CHOICE EXPERIMENT

set.seed(123)

##Functions

###Preprocessing DCE results
dce_cleanup <- function(rankdf,chardf) {
  
  scoredf <- rankdf %>% pivot_longer(Abelfenide:Protestryl)
  colnames(scoredf) <- c("Antibiotic","Rank")
  scoredf <- scoredf %>% left_join(chardf,by="Antibiotic")
  scoredf <- scoredf %>% mutate(Access = case_when(AWaRe=="Access"~1,TRUE~0),
                                Watch = case_when(AWaRe=="Watch"~1,TRUE~0),
                                Reserve = case_when(AWaRe=="Reserve"~1,TRUE~0))
  scoredf <- scoredf %>% rename(CDI_highrisk = "CDI",
                                Toxicity_highrisk = "Toxicity",
                                UTI_specific = "UTI",
                                Oral_option = "Oral",
                                IV_option = "IV",
                                High_cost = "Cost")
  scoredf[scoredf=="High"] <- "1"
  scoredf[scoredf=="Low"] <- "0"
  scoredf[scoredf=="Yes"] <- "1"
  scoredf[scoredf=="No"] <- "0"
  scoredf <- scoredf %>% select(-AWaRe) %>%
    mutate(across(c(Rank, CDI_highrisk, Toxicity_highrisk, UTI_specific,Access,
                    Watch,Reserve), as.numeric)) 
  repeat_val <- nrow(scoredf)/13
  seq_rep <- seq(1,nrow(scoredf)/13)
  seq_rep1 <- c()
  
  for (i in seq_rep) {
    
    seq_rep2 <- rep(seq_rep[i],13)
    seq_rep1 <- append(seq_rep1,seq_rep2)
    
  }
  
  scoredf %>% mutate(choice = case_when(Rank==1~1,TRUE~0),
                     id = seq_rep1) %>% 
    mutate(across(Oral_option:High_cost,as.numeric))
  
}

###Preprocessing probability/utility dataframe
utility_preprocess <- function(probsdf,scoredf,costdf,urref_df){
  
  ###CDI risk utility
  cdi_value <- scoredf[rownames(scoredf)=="CDI_highrisk",] %>% 
    select(Value) %>% unlist()
  
  ###Toxicity risk utility
  tox_value <- scoredf[rownames(scoredf)=="Toxicity_highrisk",] %>% 
    select(Value) %>% unlist()
  
  ###UTI-specific utility
  uti_specifics <- c("Nitrofurantoin")
  uti_value <- scoredf[rownames(scoredf)=="UTI_specific",] %>% 
    select(Value) %>% unlist()
  
  ###Access category utility
  access_abs <- c("AMP","SAM","CZO",
                  "GEN","SXT","NIT") %>% ab_name() %>% 
    str_replace("/","-")
  access_combos <- combn(access_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  access_abs <- c(access_abs, access_combos)
  
  access_value <- scoredf[rownames(scoredf)=="Access",] %>% 
    select(Value) %>% unlist()
  
  ###Oral option utilituy
  oral_abs <- c("AMP","SAM","CIP",
                "SXT","NIT") %>% ab_name() %>% 
    str_replace("/","-")
  oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  oral_abs <- c(oral_abs, oral_combos)
  
  oral_value <- scoredf[rownames(scoredf)=="Oral_option",] %>% 
    select(Value) %>% unlist()
  
  ###IV option utility
  iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
              "GEN","SXT","VAN") %>% ab_name() %>% 
    str_replace("/","-") 
  iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
  iv_abs <- c(iv_abs, iv_combos)
  
  iv_value <- scoredf[rownames(scoredf)=="IV_option",] %>% 
    select(Value) %>% unlist()
  
  ###Reserve category utility
  reserve_abs <- c()
  reserve_value <- scoredf[rownames(scoredf)=="Reserve",] %>% 
    select(Value) %>% unlist()
  
  ###High-cost agent utility
  highcost_abs <- c()
  cost_value <- scoredf[rownames(scoredf)=="High_cost",] %>% 
    select(Value) %>% unlist()
  
  drug_order <- costdf %>% distinct(`Generic Name`) %>% unlist()
  costdf <- costdf %>% group_by(`Generic Name`) %>% summarise(orig_cost=min(Cost)) %>% ungroup() %>% 
    mutate(`Generic Name` = factor(`Generic Name`,levels=drug_order)) %>% arrange(`Generic Name`)
  comb <- combn(costdf$`Generic Name`, 2, simplify = FALSE)
  cost_comb <- combn(costdf$orig_cost, 2, function(x) sum(x))
  costdf2 <- data.frame(
    Antimicrobial = sapply(comb, function(x) paste(x, collapse = "_")),
    orig_cost = cost_comb
  )  
  costdf <- costdf %>% rename(Antimicrobial="Generic Name")
  costdf <- data.frame(rbind(costdf,costdf2))
  costdf <- costdf %>% mutate(min_cost = orig_cost/max(orig_cost),
                              Antimicrobial = str_replace_all(Antimicrobial,"/","-"))
  probsdf <- probsdf %>% left_join(costdf)
  costdf %>% dplyr::slice(1:13) %>% summarise(
    MED = Median(orig_cost),
    iq1 = quantile(orig_cost)[2],
    iq3 = quantile(orig_cost)[4]
  )
  
  ###AST R result utility
  Rval_key <- uuref_df %>% select(micro_specimen_id,AMP_R_value:VAN_R_value)
  
  ###Attach individual utilities to dataframe
  probsdf %>% 
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
  
}

##Read_in and final preprocessing

###Read in
abx <- read_csv("interim_abx.csv")
util_probs_df <- read_csv("probs_df_overall.csv")
ur_util <- read_csv("interim_ur_util.csv")
micro <- read_csv("micro_clean2.csv")
mic_ref <- micro %>% anti_join(ur_util,by="subject_id")
rankings <- read_csv("DCE_results.csv")
characteristics <- read_csv("DCE_characteristics.csv")
cost_list <- read_csv("us_drug_cost_list.csv")

##Survey results

###Engineer scores dataframe
scores <- rankings %>% (characteristics)

##Extract characteristic importance weights

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

###Factorise results for plot
scores$Coefficient <- factor(scores$Coefficient, levels=
                         scores %>% arrange(Value) %>% 
                         select(Coefficient) %>% unlist())

###Plot DCE results
ORplot <- ggplot(scores,aes(x=Value,y=Coefficient,fill=colour)) +
  
  ###Bar chart
  geom_col() +
  
  ###Cross line
  geom_hline(aes(yintercept=0)) +
  
  ###Titles and labels
  ylab("Drug property") +
  xlab("Coefficient value for drug selection probability") +
  ggtitle("The effect of different antimicrobial drug properties on clinician prescribing\npreference in the UTI scenario discrete choice experiment")+
  
  ###Central line
  geom_vline(xintercept = 0,colour="grey40")+
  
  ###Confidence intervals
  geom_errorbar(aes(xmin = Q5, xmax = Q95), width = 0.1,colour="grey40")+
  
  ###Theme
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "None")

###Save DCE plot to pdf and print
ggsave(glue("ORplot.pdf"), plot = ORplot, device = "pdf", width = 10, height = 8,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")
print(ORplot)

##Addition to probability dataframe

###Add antimicrobials to 
util_probs_df$ab_name <- as.factor(util_probs_df$Antimicrobial)
util_probs_df <- util_probs_df %>% mutate(
  ab_name = str_replace_all(ab_name,"-",".")
)

###Dummy variables for antimicrobials
abdummy_vars <- model.matrix(~ ab_name - 1, data = util_probs_df)
util_probs_df <- cbind(util_probs_df, abdummy_vars) %>% tibble() %>% 
  select(-ab_name)

###Utility score preprocessing and write to csv
util_probs_df <- util_probs_df %>%
  utility_preprocess(scores,cost_list,ur_util)
write_csv(util_probs_df,"dce_util_probs_df.csv")


