#RECOMMENDATIONS

##Functions

###Prioritisation by treatment utility
util_mk1 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(Rx_utility)) %>% select(Antimicrobial,Rx_utility) %>% 
    mutate(Rx_utility = round(Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "Rx_utility")
  
}

###Prioritisation by AST utility
util_mk2 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(AST_utility)) %>% select(Antimicrobial,AST_utility) %>% 
    mutate(AST_utility = round(AST_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`AST Utility` = "AST_utility")
  
}

###Prioritisation by Intravenous utility
util_mk3 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>% 
    arrange(desc(Urosepsis_Rx_utility)) %>% select(Antimicrobial,Urosepsis_Rx_utility) %>% 
    mutate(Urosepsis_Rx_utility = round(Urosepsis_Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Intravenous Rx Utility` = "Urosepsis_Rx_utility")
  
}

###Prioritisation by Oral utility
util_mk4 = function(df,spec_id,panel_size,ab_list=names(combined_antimicrobial_map)) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    filter(Antimicrobial%in%ab_list) %>%
    arrange(desc(Outpatient_Rx_utility)) %>% select(Antimicrobial,Outpatient_Rx_utility) %>% 
    mutate(Outpatient_Rx_utility = round(Outpatient_Rx_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Oral Rx Utility` = "Outpatient_Rx_utility")
  
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

###Assigning AST recommendations
assign_PDAST <- function(df,probab_df,method_used,ab_list1=ab_singles) {
  
  test_recs <-  data.frame(matrix(nrow=length(ab_list1),ncol=0))
  
  for (i in 1:nrow(df)) {
    
    rec <- probab_df %>% util_mk2(spec_id = df$micro_specimen_id[i], panel_size = length(ab_list1),
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

##Read-in
ur_util <- read_csv("interim_ur_util.csv")
form_ur_util <- read_csv("interim_ur_util.csv")
util_probs_df <- read_csv("utility_dataframe.csv")
form_util_probs_df <- read_csv("form_utility_dataframe.csv")
train_abx <- read_csv("train_abx.csv")

###Map combinations to abbreviations and filter
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
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-") 
iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
iv_abs <- c(iv_abs, iv_combos)
iv_ab_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-") 
oral_abs <- c("AMP","SAM","CIP",
              "SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
oral_abs <- c(oral_abs, oral_combos)
oral_ab_singles <- c("AMP","SAM","CIP",
              "SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
ab_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
"GEN","SXT","NIT","VAN") %>% ab_name() %>% 
  str_replace("/","-")
iv_combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% iv_abs]
oral_combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% oral_abs]
ast_combined_antimicrobial_map <- combined_antimicrobial_map[!grepl("_",names(combined_antimicrobial_map))]

###Antimicrobial formulary recommendations
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
    Intravenous_Formutil = case_when(util_iv==0 ~ 0, TRUE ~ (Confirmed_S*overall_util) - (Confirmed_R*sepsis_ae)),
    Oral_Formutil = case_when(util_oral==0~0, TRUE ~ (Confirmed_S*overall_util) - (Confirmed_R*sepsis_ae)),
  ) %>% group_by(Antimicrobial) %>% 
  summarise(Formulary_util_Intravenous=mean(Intravenous_Formutil),
            Formulary_util_Oral=mean(Oral_Formutil)
  ) %>% ungroup()

form_recs_Intravenous <- form_recs %>% 
  arrange(desc(Formulary_util_Intravenous)) %>% 
  slice(1:length(iv_combined_antimicrobial_map)) %>% select(Antimicrobial) %>% unlist()

form_recs_Oral <- form_recs %>% 
  arrange(desc(Formulary_util_Oral)) %>% 
  slice(1:length(oral_combined_antimicrobial_map)) %>% select(Antimicrobial) %>% unlist()

for (i in 1:length(form_recs_Intravenous)) {
  ur_util <- ur_util %>%
    mutate(!!paste0("PDFo_Intravenous_", i) := form_recs_Intravenous[[i]])
}

for (i in 1:length(form_recs_Oral)) {
  ur_util <- ur_util %>%
    mutate(!!paste0("PDFo_Oral_", i) := form_recs_Oral[[i]])
}


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

ur_util <- ur_util %>% assign_PDRx(util_probs_df,"PDRx_") %>% 
  mutate(across(starts_with("PDRx_"), ~ replace_values(., combined_antimicrobial_map)))

###Intravenous treatment recommendations
ur_util <- ur_util %>% assign_Intravenous(util_probs_df,"Intravenous_") %>% 
  mutate(across(starts_with("Intravenous_"), ~ replace_values(., combined_antimicrobial_map)))

###Oral treatment recommendations
ur_util <- ur_util %>% assign_Oral(util_probs_df,"Oral_") %>% 
  mutate(across(starts_with("Oral_"), ~ replace_values(., combined_antimicrobial_map)))

###Individual AST recommendations
ur_util <- ur_util %>% assign_PDAST(util_probs_df,"PDAST_") %>%
  mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))

###Standard panel treatment & AST recommendations
ur_util <- ur_util %>% assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO")
ur_util <- ur_util %>% assign_standard_IV("CRO","TZP","GEN")
ur_util <- ur_util %>% assign_standard_oral("NIT","SXT","CIP")

##Formulary agent AST recommendation sensitivity analysis
replace_values <- function(column, map) {
  column %>%
    as.character() %>%
    sapply(function(x) if (x %in% names(map)) map[[x]] else x)
}
form_ur_util <- form_ur_util %>% assign_PDAST(form_util_probs_df,"PDAST_") %>%
  mutate(across(starts_with("PDAST_"), ~ replace_values(., combined_antimicrobial_map)))
form_ur_util <- form_ur_util %>% assign_standard_AST("NIT","SXT","CIP","TZP","GEN","CRO")

write_csv(ur_util,"ur_util_final.csv")
write_csv(form_ur_util,"form_ur_util_final.csv")
