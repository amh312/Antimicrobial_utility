#MICROSIMULATION SENSITIVITY ANALYSIS

##Functions

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

###Resistance weighting sensitivity analysis
res_sens_analysis <- function(df,probs_df,NEWS_variable=0,R_value=1,utility_of_interest) { 
  
  utility_of_interest <- enquo(utility_of_interest)
  
  probs_df <- probs_df %>% calculate_utilities(NEWS = NEWS_variable,R_weight = R_value)
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

iv_res <- nrow(df %>% filter(PDIVRx_1_result=='R'))
po_res <- nrow(df %>% filter(PDPORx_1_result=='R'))
overall_res <- nrow(df %>% filter(PDRx_1_result=='R'))

iv_perc <- iv_res/nrow(df)*100
po_perc <- po_res/nrow(df)*100
overall_perc <- overall_res/nrow(df)*100

iv_s_access <- (nrow(df %>% filter(PDIVRx_1_result=='S' &
                                         Intravenous_1 %in% access_singles))/
  nrow(df)) * 100
po_s_access <- (nrow(df %>% filter(PDPORx_1_result=='S' &
                                         Oral_1 %in% access_singles))/
  nrow(df))*100
overall_s_access <- (nrow(df %>% filter(PDRx_1_result=='S' &
                                          PDRx_1 %in% access_singles))/
                  nrow(df))*100
overall_s_oral <- (nrow(df %>% filter(PDRx_1_result=='S' &
                                          PDRx_1 %in% oral_singles))/
                       nrow(df))*100
overall_s_iv <- (nrow(df %>% filter(PDRx_1_result=='S' &
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

}
      
###Resistance weighting sensitivity analysis (better predictions)
res_sens_analysis_2 <- function(df,probs_df,R_variable=1,access_weight=1) { 
  
  probs_df <- probs_df %>% calculate_utilities_2(R_weight = R_variable)
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
  
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='R'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='R'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='R'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter(PDIVRx_1_result=='S' &
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter(PDPORx_1_result=='S' &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter(PDRx_1_result=='S' &
                                            PDRx_1 %in% access_singles))/
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
  
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  
}

###Resistance weighting sensitivity analysis (resistance rates increased)
res_sens_analysis_3 <- function(df,probs_df,abx,abcol,a_val,b_val) { 
  
  abcol <- enquo(abcol)
  
  probs_df <- probs_df %>% dist_replace(df,abx,"R","S",a_val,b_val) %>%
    calculate_utilities_orig()
  
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
  
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='R'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='R'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='R'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter(PDIVRx_1_result=='S' &
                                       Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter(PDPORx_1_result=='S' &
                                       Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter(PDRx_1_result=='S' &
                                            PDRx_1 %in% access_singles))/
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
  
  print(df %>% count(Intravenous_1) %>% arrange(desc(n)))
  print(df %>% count(Oral_1) %>% arrange(desc(n)))
  print(df %>% count(PDRx_1) %>% arrange(desc(n)))
  print(df %>% filter(PDRx_1%in%access_singles) %>% nrow())
  print(df %>% filter(PDRx_1%in%oral_singles) %>% nrow())
  print(df %>% filter(PDRx_1%in%iv_singles) %>% nrow())
  
}

###Susceptibility plots
susc_plotter <- function(df,subset="",measure,suffix="",agent_col1,agent_name1,agent_col2,agent_name2,variable="NEWS") {
  
  agent_col1 <- enquo(agent_col1)
  agent_col2 <- enquo(agent_col2)
  
  df$Metric <- factor(df$Metric, levels=c("All agents","Access agents"))
  
  susplot <- ggplot(df,aes(x=Weight,y=Percentage,group=Metric,color=Metric,fill=Metric)) +
    geom_area(position = "identity", alpha = 0.6) +
    geom_hline(linetype="dashed",yintercept=(nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100,color="gray16") +
    geom_hline(linetype="dashed",yintercept=70,color="darkgreen") +
    ylim(0,100)+
    theme_minimal()+
    ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}"))+
    xlab(glue("{variable}"))+
    ylab("Percentage of isolates susceptible to recommendation")+
    annotate("text", x = Inf, y = (nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100+2, label = glue(agent_name1), hjust = 1.1, color = "gray16")+
    annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen")
  
  ggsave(glue("{subset}_{measure}_{suffix}_susplot.pdf"), plot = susplot, device = "pdf", width = 10, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(susplot)
  
}
susc_plotter_iv <- function(df,subset="",measure,suffix="",agent_col1,agent_name1,variable="NEWS") {
  
  agent_col1 <- enquo(agent_col1)
  agent_col2 <- enquo(agent_col2)
  
  df$Metric <- factor(df$Metric, levels=c("All agents","Access agents"))
  
  susplot <- ggplot(df,aes(x=Weight,y=Percentage,group=Metric,color=Metric,fill=Metric)) +
    geom_area(position = "identity", alpha = 0.6) +
    geom_hline(linetype="dashed",yintercept=(nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100,color="gray16") +
    geom_hline(linetype="dashed",yintercept=70,color="darkgreen") +
    ylim(0,100)+
    theme_minimal()+
    ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}"))+
    xlab(glue("{variable}"))+
    ylab("Percentage of isolates susceptible to recommendation")+
    annotate("text", x = Inf, y = (nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100+2, label = glue(agent_name1), hjust = 1.1, color = "gray16")+
    annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen")
  
  ggsave(glue("{subset}_{measure}_{suffix}_susplot.pdf"), plot = susplot, device = "pdf", width = 10, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
  print(susplot)
  
}
susc_plotter_overall <- function(df,subset="",measure,suffix="",agent_col1,agent_name1,agent_col2,agent_name2,variable="NEWS") {
  
  agent_col1 <- enquo(agent_col1)
  agent_col2 <- enquo(agent_col2)
  
  df$Metric <- factor(df$Metric, levels=c("All agents","Access agents","IV agents","Oral agents"))
  
  susplot <- ggplot(df,aes(x=Weight,y=Percentage,group=Metric,color=Metric,fill=Metric)) +
    geom_area(position = "identity", alpha = 0.6) +
    geom_hline(linetype="dashed",yintercept=(nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100,color="gray16") +
    geom_hline(linetype="dashed",yintercept=(nrow(ur_util %>% filter(!!agent_col2=="S"|!!agent_col2=="I"))/nrow(ur_util))*100,color="gray16") +
    geom_hline(linetype="dashed",yintercept=70,color="darkgreen") +
    ylim(0,100)+
    theme_minimal()+
    ggtitle(glue("Urine isolate susceptibility to {subset}antimicrobial recommended\n{suffix}"))+
    xlab(glue("{variable}"))+
    ylab("Percentage of isolates susceptible to recommendation")+
    annotate("text", x = Inf, y = (nrow(ur_util %>% filter(!!agent_col1=="S"|!!agent_col1=="I"))/nrow(ur_util))*100+2, label = glue(agent_name1), hjust = 1.1, color = "gray16")+
    annotate("text", x = Inf, y = (nrow(ur_util %>% filter(!!agent_col2=="S"|!!agent_col2=="I"))/nrow(ur_util))*100+2, label = glue(agent_name2), hjust = 1.1, color = "gray16")+
    annotate("text", x = Inf, y = 72, label = "UN Access target", hjust = 1.1, color = "darkgreen")
  
  ggsave(glue("{subset}_{measure}_{suffix}_susplot.pdf"), plot = susplot, device = "pdf", width = 10, height = 6,
         path="/Users/alexhoward/Documents/Projects/UDAST_code")
  
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
combo_res_sens_analysis <- function(df,probs_df,NEWS_variable=0,R_value=1) { 
  
  probs_df <- probs_df %>% calculate_utilities(NEWS = NEWS_variable,R_weight = R_value)
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
  
  iv_res <- nrow(df %>% filter(PDIVRx_1_result=='R'))
  po_res <- nrow(df %>% filter(PDPORx_1_result=='R'))
  overall_res <- nrow(df %>% filter(PDRx_1_result=='R'))
  
  iv_perc <- iv_res/nrow(df)*100
  po_perc <- po_res/nrow(df)*100
  overall_perc <- overall_res/nrow(df)*100
  
  iv_s_access <- (nrow(df %>% filter(PDIVRx_1_result=='S' &
                                            Intravenous_1 %in% access_singles))/
                    nrow(df)) * 100
  po_s_access <- (nrow(df %>% filter(PDPORx_1_result=='S' &
                                            Oral_1 %in% access_singles))/
                    nrow(df))*100
  overall_s_access <- (nrow(df %>% filter(PDRx_1_result=='S' &
                                                 PDRx_1 %in% access_singles))/
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
calculate_utilities <- function(df,formulary_list=c(),NEWS=0,R_weight=1) {
  
  df <- df %>% mutate(overall_util=util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI + NEWS*R_weight,
                      S_utility = S*overall_util,
                      ent_S_utility = ent_S*overall_util,
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
                      AST_utility = ((S_utility+
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
                                        Formulary_utility)*single_agent),
                      Rx_utility = S_utility,
                      ent_AMPR_utility = case_when(
                        Antimicrobial=="Ampicillin" ~
                          abs(R_weight)*AMP_R_value*ent_R, TRUE~0
                      ),
                      ent_SAMR_utility = case_when(
                        Antimicrobial=="Ampicillin-sulbactam" ~
                          abs(R_weight)*SAM_R_value*ent_R, TRUE~0
                      ),
                      ent_TZPR_utility = case_when(
                        Antimicrobial=="Piperacillin-tazobactam" ~
                          abs(R_weight)*TZP_R_value*ent_R, TRUE~0
                      ),
                      ent_CZOR_utility = case_when(
                        Antimicrobial=="Cefazolin" ~
                          abs(R_weight)*CZO_R_value*ent_R, TRUE~0
                      ),
                      ent_CROR_utility = case_when(
                        Antimicrobial=="Ceftriaxone" ~
                          abs(R_weight)*CRO_R_value*ent_R, TRUE~0
                      ),
                      ent_CAZR_utility = case_when(
                        Antimicrobial=="Ceftazidime" ~
                          abs(R_weight)*CAZ_R_value*ent_R, TRUE~0
                      ),
                      ent_FEPR_utility = case_when(
                        Antimicrobial=="Cefepime" ~
                          abs(R_weight)*FEP_R_value*ent_R, TRUE~0
                      ),
                      ent_MEMR_utility = case_when(
                        Antimicrobial=="Meropenem" ~
                          abs(R_weight)*MEM_R_value*ent_R, TRUE~0
                      ),
                      ent_CIPR_utility = case_when(
                        Antimicrobial=="Ciprofloxacin" ~
                          abs(R_weight)*CIP_R_value*ent_R, TRUE~0
                      ),
                      ent_GENR_utility = case_when(
                        Antimicrobial=="Gentamicin" ~
                          abs(R_weight)*GEN_R_value*ent_R, TRUE~0
                      ),
                      ent_SXTR_utility = case_when(
                        Antimicrobial=="Trimethoprim-sulfamethoxazole" ~
                          abs(R_weight)*SXT_R_value*ent_R, TRUE~0
                      ),
                      ent_NITR_utility = case_when(
                        Antimicrobial=="Nitrofurantoin" ~
                          abs(R_weight)*NIT_R_value*ent_R, TRUE~0
                      ),
                      ent_VANR_utility = case_when(
                        Antimicrobial=="Vancomycin" ~
                          abs(R_weight)*VAN_R_value*ent_R, TRUE~0
                      ),
                      ent_Formulary_agent = case_when(
                        Antimicrobial%in%formulary_list ~
                          TRUE, TRUE~FALSE
                      ),
                      ent_Formulary_utility = 
                        abs(R_weight)*Formulary_agent*ent_R,
                      ent_AST_utility = ((ent_S_utility+
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
                                            ent_Formulary_utility)*single_agent),
                      ent_Rx_utility = (overall_util * ent_S))
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~0,
        TRUE~Rx_utility),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~0,
        TRUE~Rx_utility),
      ent_Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~0,
        TRUE~(overall_util*ent_S)),
      ent_Outpatient_Rx_utility = case_when(
        util_oral ==0 ~0,
        TRUE~(overall_util*ent_S))
    )
  
}
calculate_utilities_2 <- function(df,formulary_list=c(),R_weight=1) {
  
  df <- df %>% mutate(overall_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      overall_oral_util = util_uti + util_access +
                        util_oral +
                        util_reserve + util_highcost
                      + util_tox + util_CDI,
                      overall_iv_util = util_uti + util_access +
                        util_iv +
                        util_reserve + util_highcost
                      + util_tox + util_CDI,
                      R_penalty=R_weight*better_R*prob_sepsisae,
                      S_utility = S*overall_util,
                      S_PO_utility = S*overall_oral_util,
                      S_IV_utility=S*overall_iv_util,
                      ent_S_utility = ent_S*overall_util,
                      AMPR_utility = case_when(
                        Antimicrobial=="Ampicillin" ~
                          AMP_R_value*better_R, TRUE~0
                      ),
                      SAMR_utility = case_when(
                        Antimicrobial=="Ampicillin-sulbactam" ~
                          SAM_R_value*better_R, TRUE~0
                      ),
                      TZPR_utility = case_when(
                        Antimicrobial=="Piperacillin-tazobactam" ~
                          TZP_R_value*better_R, TRUE~0
                      ),
                      CZOR_utility = case_when(
                        Antimicrobial=="Cefazolin" ~
                          CZO_R_value*better_R, TRUE~0
                      ),
                      CROR_utility = case_when(
                        Antimicrobial=="Ceftriaxone" ~
                          CRO_R_value*better_R, TRUE~0
                      ),
                      CAZR_utility = case_when(
                        Antimicrobial=="Ceftazidime" ~
                          CAZ_R_value*better_R, TRUE~0
                      ),
                      FEPR_utility = case_when(
                        Antimicrobial=="Cefepime" ~
                          FEP_R_value*better_R, TRUE~0
                      ),
                      MEMR_utility = case_when(
                        Antimicrobial=="Meropenem" ~
                          MEM_R_value*better_R, TRUE~0
                      ),
                      CIPR_utility = case_when(
                        Antimicrobial=="Ciprofloxacin" ~
                          CIP_R_value*better_R, TRUE~0
                      ),
                      GENR_utility = case_when(
                        Antimicrobial=="Gentamicin" ~
                          GEN_R_value*better_R, TRUE~0
                      ),
                      SXTR_utility = case_when(
                        Antimicrobial=="Trimethoprim-sulfamethoxazole" ~
                          SXT_R_value*better_R, TRUE~0
                      ),
                      NITR_utility = case_when(
                        Antimicrobial=="Nitrofurantoin" ~
                          NIT_R_value*better_R, TRUE~0
                      ),
                      VANR_utility = case_when(
                        Antimicrobial=="Vancomycin" ~
                          VAN_R_value*better_R, TRUE~0
                      ),
                      Formulary_agent = case_when(
                        Antimicrobial%in%formulary_list ~
                          TRUE, TRUE~FALSE
                      ),
                      Formulary_utility = 
                        R_weight*Formulary_agent*better_R,
                      AST_utility = normalise(((S_utility+
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
                                                  Formulary_utility)*single_agent)),
                      Rx_utility = normalise((S_utility -
                                                R_penalty)),
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
                        R_weight*Formulary_agent*ent_R,
                      ent_AST_utility = normalise(((ent_S_utility+
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
                                                      ent_Formulary_utility)*single_agent)),
                      ent_Rx_utility = normalise((overall_util * ent_S) -
                                                   (R_weight*ent_R*(prob_sepsisae))))
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise(S_IV_utility - R_penalty)),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise(S_PO_utility -  R_penalty)),
      ent_Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise((overall_iv_util*ent_S)-(R_weight*ent_R*(prob_sepsisae - util_CDI -
                                                                  util_tox)))),
      ent_Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise((overall_oral_util*ent_S)-(R_weight*ent_R*(prob_sepsisae - util_CDI -
                                                                    util_tox))))
    )
  
}
calculate_utilities_orig <- function(df,formulary_list=c(),R_weight=1) {
  
  df <- df %>% mutate(overall_util = util_uti + util_access +
                        util_oral + util_iv +
                        util_reserve + util_highcost 
                      + util_tox + util_CDI,
                      overall_oral_util = util_uti + util_access +
                        util_oral +
                        util_reserve + util_highcost
                      + util_tox + util_CDI,
                      overall_iv_util = util_uti + util_access +
                        util_iv +
                        util_reserve + util_highcost
                      + util_tox + util_CDI,
                      R_penalty=R_weight*R*prob_sepsisae,
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
                        R_weight*Formulary_agent*R,
                      AST_utility = normalise(((S_utility+
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
                                                  Formulary_utility)*single_agent)),
                      Rx_utility = normalise((S_utility -
                                                R_penalty)),
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
                        R_weight*Formulary_agent*ent_R,
                      ent_AST_utility = normalise(((ent_S_utility+
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
                                                      ent_Formulary_utility)*single_agent)),
                      ent_Rx_utility = normalise((overall_util * ent_S) -
                                                   (R_weight*ent_R*(prob_sepsisae))))
  
  df %>% 
    mutate(
      Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise(S_IV_utility - R_penalty)),
      Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise(S_PO_utility -  R_penalty)),
      ent_Urosepsis_Rx_utility = case_when(
        util_iv ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise((overall_iv_util*ent_S)-(R_weight*ent_R*(prob_sepsisae - util_CDI -
                                                                  util_tox)))),
      ent_Outpatient_Rx_utility = case_when(
        util_oral ==0 ~min(df$Rx_utility)-0.01,
        TRUE~normalise((overall_oral_util*ent_S)-(R_weight*ent_R*(prob_sepsisae - util_CDI -
                                                                    util_tox))))
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

###Search for previous event across multiple variables
apply_prev_event <- function(df, param,organism) {
  df %>%
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, 365, 1)
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  params <- paste0("p", antibiotics, suffix)
  
  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_event(df, params[i], antibiotics[i])
  }, .init = df1) %>%
    ungroup()
  
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

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  params <- paste0("p", antibiotics, suffix)
  
  apply_prev_event <- function(df, param, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(param), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  df1 <- reduce(seq_along(antibiotics), function(df, i) {
    apply_prev_event(df, params[i], antibiotics[i])
  }, .init = df1) %>%
    ungroup()
  
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(event_type, event)) %>%
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Assigning previous treatment
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx_groupvar <- enquo(abx_groupvar)
  
  drug_df %>%
    select('subject_id', ab_name,charttime='starttime') %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    bind_rows(ur_df) %>% 
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df,icd_df,prefix,codes,time_frame) {
  
  apply_prev_event_assignments <- function(df, code) {
    param_name <- paste0(prefix, code)
    df %>%
      prev_event_type_assign(!!sym(param_name), icd_df, icd_group, code, time_frame, 1)
  }
  
  df <- reduce(codes, function(df, code) {
    apply_prev_event_assignments(df, code)
  }, .init = df) %>%
    mutate(pDIAG_U = FALSE) %>%
    ungroup()
  
}

###Assigning previous event feature variable
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    mutate(event = {{event_var}}) %>% 
    select('subject_id', "event", charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(!is.na(event)) %>% 
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

###Search for previous event across multiple variables
apply_prev_event <- function(df, param,organism,time_frame=365,no_events=1) {
  df %>%
    prev_event_type_assign(!!sym(param), urine_df, org_fullname,organism, time_frame, no_events)
}

###Checking for previous care events
care_event_assigner <- function(df,search_df,search_term,search_column,feature_name,event_date_col,timeframe,n_events=1) {
  
  feature_name <- enquo(feature_name)
  search_column <- enquo(search_column)
  
  care_event <- search_df %>% filter(grepl(search_term,!!search_column,ignore.case=T)) %>% mutate(
    !!search_column:=search_term) %>% rename(admittime=event_date_col)
  df %>% 
    prev_event_type_assign(!!feature_name,care_event,!!search_column,search_term,timeframe,n_events) %>%
    ungroup()
  
}

###Reference lists
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                 "MEM","CIP","GEN","SXT","NIT")
ab_singles <- all_singles
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                "GEN","SXT")
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

###AST

##Sensitivity analysis for inclusion of ciprofloxacin and ceftriaxone as formulary agents

###Read_in
form_ur_util <- read_csv("form_ur_util_final.csv")
util_probs_df <- read_csv("pre_util_probs_df.csv")
ur_util <- read_csv("interim_ur_util.csv")
train_abx <- read_csv("train_abx.csv")
micro <- read_csv("micro_clean2.csv")
drugs <- read_csv("drugs_clean.csv")
hadm <- read_csv("admissions.csv") #Admission data
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv") #icd codes
diagnoses_raw <- read_csv("diagnoses_icd.csv") #icd epi
diagnoses <- read_csv("diagnoses_clean.csv")
procedures <- read_csv("procedures_clean.csv")
urine_df <- read_csv("pos_urines.csv")
poe <- read_csv("poe_clean.csv")

###Count numbers of S and R results per antimicrobial agent
form_pdast_all_singles <- form_ur_util %>% number_ab_results(PDAST_1,PDAST_6,all_singles,"S","I")
form_standard_all_singles <- form_ur_util %>% number_ab_results(STANDARD_AST_1,STANDARD_AST_6,all_singles,"S","I")
form_pdast_all_singles_r <- form_ur_util %>% number_ab_results(PDAST_1,PDAST_6,all_singles,"R","NT")
form_standard_all_singles_r <- form_ur_util %>% number_ab_results(STANDARD_AST_1,STANDARD_AST_6,all_singles,"R","NT")

###Plot S and R results per antimicrobial agent
form_abs_df <- abs_df_assemble(form_pdast_all_singles,form_standard_all_singles)
write_csv(form_abs_df,"form_s_uf_sourcedata_abs_cleveplot.csv")
form_s_results_by_ab <- form_abs_df %>% cleveland_ab_plot("susceptible"," (Ciprofloxacin & Ceftriaxone added to formulary)")
form_abs_df_r <- abs_df_assemble(form_pdast_all_singles_r,form_standard_all_singles_r)
write_csv(form_abs_df_r,"form_r_uf_sourcedata_abs_cleveplot.csv")
form_r_results_by_ab <- form_abs_df_r %>% cleveland_ab_plot("resistant"," (Ciprofloxacin & Ceftriaxone added to formulary)")

###R weight sensitivity analysis
weightseq <- seq(0,20,1)
iv_perclist <- c()
po_perclist <- c()
overall_perclist <- c()
ivac_list <- c()
poac_list <- c()
ovac_list <- c()

for(weight in seq_along(weightseq)) {
  
  res_sens_analysis(ur_util,util_probs_df,weightseq[weight])
  
  iv_perclist <- c(iv_perclist,iv_perc)
  po_perclist <- c(po_perclist,po_perc)
  overall_perclist <- c(overall_perclist,overall_perc)
  ivac_list <- c(ivac_list,iv_s_access)
  poac_list <- c(poac_list,po_s_access)
  ovac_list <- c(ovac_list,overall_s_access)
  
}

label_binder <- function(vec,label) {
  
  data.frame(vec,Metric=label,Weight=weightseq)
  
}

iv_perclist <- iv_perclist %>% label_binder("All agents")
po_perclist <- po_perclist %>% label_binder("All agents")
overall_perclist <- overall_perclist %>% label_binder("All agents")
ivac_list <- ivac_list %>% label_binder("Access agents")
poac_list <- poac_list %>% label_binder("Access agents")
ovac_list <- ovac_list %>% label_binder("Access agents")
iv_perclist <- iv_perclist %>% rename(Percentage = "vec")
po_perclist <- po_perclist %>% rename(Percentage = "vec")
overall_perclist <- overall_perclist %>% rename(Percentage = "vec")
ivac_list <- ivac_list %>% rename(Percentage = "vec")
poac_list <- poac_list %>% rename(Percentage = "vec")
ovac_list <- ovac_list %>% rename(Percentage = "vec")
iv_perclist <- iv_perclist %>% mutate(Percentage=100-Percentage)
po_perclist <- po_perclist %>% mutate(Percentage=100-Percentage)
overall_perclist <- overall_perclist %>% mutate(Percentage=100-Percentage)

iv_plot_df <- data.frame(rbind(
  iv_perclist,ivac_list
))
po_plot_df <- data.frame(rbind(
  po_perclist,poac_list
))
overall_plot_df <- data.frame(rbind(
  overall_perclist,ovac_list
))

write_csv(iv_plot_df,"iv_plot_df.csv")
write_csv(po_plot_df,"po_plot_df.csv")
write_csv(overall_plot_df,"overall_plot_df.csv")

iv_plot_df %>% susc_plotter_iv("IV ","NEWS",agent_col1=GEN,agent_name1="Gentamicin",agent_col2=TZP,agent_name2="Piperacillin-tazobactam")
po_plot_df %>% susc_plotter("oral ","NEWS",agent_col1=NIT,agent_name1="Nitrofurantoin")
overall_plot_df %>% susc_plotter_overall("overall ", "NEWS",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                    agent_col2=TZP,agent_name2="Piperacillin-tazobactam")


###Nitro Resistance sensitivity analysis
weightseq <- c()
iv_perclist_nit <- c()
po_perclist_nit <- c()
overall_perclist_nit <- c()
ivac_list_nit <- c()
poac_list_nit <- c()
ovac_list_nit <- c()
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
  
  a <- a+i
  b <- b-(i/2)
  
}
randomLHS(6,4)

res_sens_analysis_3(ur_util,util_probs_df,abx="Nitrofurantoin",NIT,a_val=9.5,b_val=14)
iv_perclist_nit <- iv_perclist_nit %>% label_binder("All agents")
po_perclist_nit <- po_perclist_nit %>% label_binder("All agents")
overall_perclist_nit <- overall_perclist_nit %>% label_binder("All agents")
ivac_list_nit <- ivac_list_nit %>% label_binder("Access agents")
poac_list_nit <- poac_list_nit %>% label_binder("Access agents")
ovac_list_nit <- ovac_list_nit %>% label_binder("Access agents")
iv_perclist_nit <- iv_perclist_nit %>% rename(Percentage = "vec")
po_perclist_nit <- po_perclist_nit %>% rename(Percentage = "vec")
overall_perclist_nit <- overall_perclist_nit %>% rename(Percentage = "vec")
ivac_list_nit <- ivac_list_nit %>% rename(Percentage = "vec")
poac_list_nit <- poac_list_nit %>% rename(Percentage = "vec")
ovac_list_nit <- ovac_list_nit %>% rename(Percentage = "vec")
iv_perclist_nit <- iv_perclist_nit %>% mutate(Percentage=100-Percentage)
po_perclist_nit <- po_perclist_nit %>% mutate(Percentage=100-Percentage)
overall_perclist_nit <- overall_perclist_nit %>% mutate(Percentage=100-Percentage)

iv_plot_df_nit <- data.frame(rbind(
  iv_perclist_nit,ivac_list_nit
))
po_plot_df_nit <- data.frame(rbind(
  po_perclist_nit,poac_list_nit
))
overall_plot_df_nit <- data.frame(rbind(
  overall_perclist_nit,ovac_list_nit
))

write_csv(iv_plot_df_nit,"iv_plot_df_nit.csv")
write_csv(po_plot_df_nit,"po_plot_df_nit.csv")
write_csv(overall_plot_df_nit,"overall_plot_df_nit.csv")

iv_plot_df_nit %>% susc_plotter_iv("IV ","Resistance",agent_col1=GEN,agent_name1="Gentamicin",agent_col2=TZP,agent_name2="Piperacillin-tazobactam")
po_plot_df_nit %>% susc_plotter("oral ","Resistance",agent_col1=TZP,agent_name1="Piperacillin-tazobactam",variable="Nitrofurantoin resistance rate")
overall_plot_df_nit %>% susc_plotter("overall ", "Resistance",agent_col1=TZP,agent_name1="Piperacillin-tazobactam",variable="Nitrofurantoin resistance rate")



###XGBoost sensitivity analysis
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                 "MEM","CIP","GEN","SXT","NIT")
ab_singles <- all_singles
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                "GEN","SXT")
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
    sapply(map_combinations, function(x) paste(x, collapse = "_"))
  )
)
abx_in_train <- train_abx %>% distinct(abx_name) %>% unlist() %>% 
  str_replace_all("/","-")
fullmap <- combined_antimicrobial_map %>% unlist()
names(fullmap) <- NULL
combined_antimicrobial_map <- combined_antimicrobial_map[names(combined_antimicrobial_map) %in% abx_in_train]
shortmap <- combined_antimicrobial_map %>% unlist()
names(shortmap) <- NULL

urines5 <- read_csv("urines5.csv")
ur_xg <- read_csv("ur_util.csv")
urines5 <- urines5 %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                       TRUE~marital_status))
ur_xg <- ur_xg %>% mutate(marital_status=case_when(is.na(marital_status)~"UNKNOWN",
                                                   TRUE~marital_status))

urines5_outcomes <- urines5 %>%
  select(all_of(shortmap))
ur_xg_outcomes <- ur_xg %>%
  select(all_of(shortmap))

urines5_outcomes <- urines5_outcomes %>%
  mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                 ifelse(. == "S" | . == "I", 1, NA))))
ur_xg_outcomes <- ur_xg_outcomes %>%
  mutate_all(~ as.numeric(ifelse(. == "R" | . == "NT", 0, 
                                 ifelse(. == "S" | . == "I", 1, NA))))

urines5_predictors <- urines5 %>% select(!all_of(fullmap))
ur_xg_predictors <- ur_xg %>%
  select(any_of(names(urines5_predictors)))

set.seed(123)
dummies <- dummyVars(" ~ .", data = urines5_predictors)
urines5_predictors <- predict(dummies, newdata = urines5_predictors)
dummies2 <- dummyVars(" ~ .", data = ur_xg_predictors)
ur_xg_predictors <- predict(dummies2, newdata = ur_xg_predictors)

urines5_combined <- as.data.frame(cbind(urines5_outcomes, urines5_predictors))
ur_xg_combined <- as.data.frame(cbind(ur_xg_outcomes, ur_xg_predictors))

###First round of hyperparameter tuning (max tree depth and min child weight)

num_samples <- 10
max_depth_range <- c(2,9)
min_child_weight_range <- c(1, 10)
lhs_sample <- randomLHS(num_samples, 2)
max_depth <- round(lhs_sample[, 1] * (max_depth_range[2] - max_depth_range[1]) + max_depth_range[1])
min_child_weight <- round(lhs_sample[, 2] * (min_child_weight_range[2] - min_child_weight_range[1]) + min_child_weight_range[1])
parameter_grid <- data.frame(max_depth = max_depth, min_child_weight = min_child_weight)
print(parameter_grid)
max_child_bestparams <- c()
best_auc <- 0

for (outcome in colnames(urines5_outcomes)) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    missing_cols <- setdiff(selected_columns, colnames(ur_xg_combined))
    ur_xg_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(all_of(selected_columns))), 
                                label = ur_xg_combined[[outcome]])
    
    for (i in 1:nrow(parameter_grid)) {
    
    print(glue("Running CV {i} for {outcome}..."))
    
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = parameter_grid %>% select(max_depth) %>% dplyr::slice(i) %>% unlist(),
        min_child_weight = parameter_grid %>% select(min_child_weight) %>% dplyr::slice(i) %>% unlist(),
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
    cv_model <- xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = 50,
      nfold = 5,
      early_stopping_rounds = 50,
      verbose = 1,
    )
    
    best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
    best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
    cv_model$evaluation_log$test_logloss_mean
    if (best_iteration_auc > best_auc) {
      best_auc <- best_iteration_auc
      best_params <- params
      best_nrounds <- best_iteration_index
    }
    
    }
    
    max_child_bestparams[[outcome]] <- best_params
    
  }
}


for (i in 1:length(max_child_bestparams)) {
  
  maxy <- data.frame(max_child_bestparams[i])
  
  write_csv(maxy,glue("max_child_{combined_antimicrobial_map[i]}.csv"))
  
}

for (outcome in colnames(urines5_outcomes)){
print(max_child_bestparams[[outcome]]$min_child_weight)
}

###Second round of hyperparameter tuning

num_samples <- 10
subsample_range <- c(0.5,1)
colsample_bytree_range <- c(0.5, 1)
lhs_sample <- randomLHS(num_samples, 2)
subsample <- round(lhs_sample[, 1] * (subsample_range[2] - subsample_range[1]) + subsample_range[1],2)
colsample_bytree <- round(lhs_sample[, 2] * (colsample_bytree_range[2] - colsample_bytree_range[1]) + colsample_bytree_range[1],2)
parameter_grid <- data.frame(subsample = subsample, colsample_bytree = colsample_bytree)
print(parameter_grid)
col_sub_bestparams <- c()
best_auc <- 0

for (outcome in colnames(urines5_outcomes)) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    missing_cols <- setdiff(selected_columns, colnames(ur_xg_combined))
    ur_xg_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(all_of(selected_columns))), 
                                label = ur_xg_combined[[outcome]])
    
    for (i in 1:nrow(parameter_grid)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.05,
        max_depth = max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams[[outcome]]$min_child_weight,
        subsample = parameter_grid %>% select(subsample) %>% dplyr::slice(i) %>% unlist(),
        colsample_bytree = parameter_grid %>% select(colsample_bytree) %>% dplyr::slice(i) %>% unlist()
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = train_matrix,
        nrounds = 50,
        nfold = 5,
        early_stopping_rounds = 50,
        verbose = 1,
      )
      
      best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
      best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
      cv_model$evaluation_log$test_logloss_mean
      if (best_iteration_auc > best_auc) {
        best_auc <- best_iteration_auc
        best_params <- params
        best_nrounds <- best_iteration_index
      }
      
    }
    
    col_sub_bestparams[[outcome]] <- best_params
    
  }
}

for (i in 1:length(col_sub_bestparams)) {
  
  coly <- data.frame(col_sub_bestparams[i])
  
  write_csv(coly,glue("col_sub_{combined_antimicrobial_map[i]}.csv"))
  
}

###Third round of hyperparameter tuning

parameter_list <- c(0.1,0.05,0.01,0.001)
best_auc <- 0
final_bestparams <- c()

for (outcome in colnames(urines5_outcomes)) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    missing_cols <- setdiff(selected_columns, colnames(ur_xg_combined))
    ur_xg_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(all_of(selected_columns))), 
                                label = ur_xg_combined[[outcome]])
    
    for (i in 1:length(parameter_list)) {
      
      print(glue("Running CV {i} for {outcome}..."))
      
      params <- list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = parameter_list[i],
        max_depth = max_child_bestparams[[outcome]]$max_depth,
        min_child_weight = max_child_bestparams[[outcome]]$min_child_weight,
        subsample = col_sub_bestparams[[outcome]]$subsample,
        colsample_bytree = col_sub_bestparams[[outcome]]$colsample_bytree
      )
      
      cv_model <- xgb.cv(
        params = params,
        data = train_matrix,
        nrounds = 1000,
        nfold = 5,
        early_stopping_rounds = 50,
        verbose = 1,
      )
      
      best_iteration_index <- which.max(cv_model$evaluation_log$test_auc_mean)
      best_iteration_auc <- cv_model$evaluation_log$test_auc_mean[best_iteration_index]
      cv_model$evaluation_log$test_logloss_mean
      if (best_iteration_auc > best_auc) {
        best_auc <- best_iteration_auc
        best_params <- params
        best_nrounds <- best_iteration_index
      }
      
    }
    
    final_bestparams[[outcome]] <- best_params
    
  }
}

for (i in 1:length(final_bestparams)) {
  
  param <- data.frame(final_bestparams[i])
  
  write_csv(param,glue("final_params_{combined_antimicrobial_map[i]}.csv"))
  
}


###Tuning time to last resistance
urines_ref <- read_csv("urines_ref.csv")
ur_util_ref <- ur_util
micro3 <- micro %>% rename(admittime = "charttime")

time_list <- c(30,180,720,1e4)

for (i in 1:length(time_list)) {
  
  antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                   "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                   "CLI", "LVX", "AMK", "TOB")
  urines_ref <- prev_AST_applier(urines_ref,micro3,glue("r_{as.character(time_list[i])}"),"R",timeframe=time_list[i])
  ur_util <- prev_AST_applier(ur_util,micro3,glue("r_{as.character(time_list[i])}"),"R",timeframe=time_list[i])
  
}

new_r_cols <- urines_ref %>% select(pAMPr_30:pTOBr_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_r_cols))

time_list <- c(7,30,180,720,1e4)

###Tuning time to last susceptibility
for (i in 1:length(time_list)) {
  
  antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                   "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                   "CLI", "LVX", "AMK", "TOB")
  urines_ref <- prev_AST_applier(urines_ref,micro3,glue("s_{as.character(time_list[i])}"),"S",timeframe=time_list[i])
  ur_util <- prev_AST_applier(ur_util,micro3,glue("s_{as.character(time_list[i])}"),"S",timeframe=time_list[i])
  
}

new_s_cols <- urines_ref %>% select(pAMPs_30:pTOBs_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_s_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###Tuning time to last treatment
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin", "Cefalexin", "Cefpodoxime proxetil",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- c("AMPrx", "AMXrx", "AMCrx", "SAMrx", "TZPrx", "CZOrx", "CZOrx", "CZOrx",
              "CROrx", "CAZrx", "FEPrx", "MEMrx", "ETPrx", "ATMrx", "CIPrx", "CIPrx",
              "GENrx", "TOBrx", "AMKrx", "RIFrx", "SXTrx", "NITrx", "ERYrx", "CLRrx",
              "AZMrx", "CLIrx", "VANrx", "MTRrx", "LNZrx", "DAPrx", "DOXrx")

apply_prev_rx <- function(df, suffix, antibiotic,time_to_event=365) {
  param_name <- paste0("p", suffix,"_",time_list[j])
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, time_to_event, 1)
}

time_list <- c(30,180,720,1e4)
drugs$ab_name <- drugs$abx_name

for (j in 1:length(time_list)) {
urines_ref <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],time_list[j])
}, .init = urines_ref) %>%
  ungroup()
ur_util <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],time_list[j])
}, .init = ur_util) %>%
  ungroup()
}

new_rx_cols <- urines_ref %>% select(pAMPrx_30:pDOXrx_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_rx_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###2 episodes of resistance in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                 "CLI", "LVX", "AMK", "TOB")
urines_ref <- prev_AST_applier(urines_ref,micro3,"r2","R",timeframe=365,n_events=2)
ur_util <- prev_AST_applier(ur_util,micro3,"r2","R",timeframe=365,n_events = 2)

new_2r_cols <- urines_ref %>% select(pAMPr2:pTOBr2)
urines5_combined <- data.frame(cbind(urines5_combined,new_2r_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###2 antimicrobial courses in the last year
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin", "Cefalexin", "Cefpodoxime proxetil",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- c("AMPrx", "AMXrx", "AMCrx", "SAMrx", "TZPrx", "CZOrx", "CZOrx", "CZOrx",
              "CROrx", "CAZrx", "FEPrx", "MEMrx", "ETPrx", "ATMrx", "CIPrx", "CIPrx",
              "GENrx", "TOBrx", "AMKrx", "RIFrx", "SXTrx", "NITrx", "ERYrx", "CLRrx",
              "AZMrx", "CLIrx", "VANrx", "MTRrx", "LNZrx", "DAPrx", "DOXrx")

apply_prev_rx <- function(df, suffix, antibiotic,time_to_event=365,no_events=2) {
  param_name <- paste0("p", suffix,"2")
  df %>%
    prev_rx_assign(!!sym(param_name), drugs, antibiotic, ab_name, time_to_event, no_events)
}

urines_ref <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],365,2)
}, .init = urines_ref) %>%
  ungroup()
ur_util <- reduce(seq_along(antibiotics), function(df, i) {
  apply_prev_rx(df, suffixes[i], antibiotics[i],365,2)
}, .init = ur_util) %>%
  ungroup()

new_rx_cols <- urines_ref %>% select(pAMPrx2:pDOXrx2)
urines5_combined <- data.frame(cbind(urines5_combined,new_rx_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###Tuning hospital admission
time_list <- c(30,180,720,1e4)

for (i in 1:length(time_list)) {
urines_ref <- urines_ref %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,time_list[i],1) %>%
  ungroup()
colnames(urines_ref)[ncol(urines_ref)] <- glue("pHADM_{time_list[i]}")
ur_util <- ur_util %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,time_list[i],1) %>%
  ungroup()
colnames(ur_util)[ncol(urines_ref)] <- glue("pHADM_{time_list[i]}")
}

new_hadm_cols <- urines_ref %>% select(pHADM_10000:pHADM_720)
urines5_combined <- data.frame(cbind(urines5_combined,new_hadm_cols))
urines_ref <- urines_ref %>% select(-pHADM_10000)
write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###2 hospital admissions in the last year
ur_util <- ur_util %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,365,2) %>%
  ungroup()
urines_ref <- urines_ref %>% 
  prev_event_assign(pHADM2,hadm,hadm_id,365,2) %>%
  ungroup()

hadm2 <- urines_ref %>% select(pHADM2)
urines5_combined <- data.frame(cbind(urines5_combined,hadm2))

###Tuning coded ICD-10 diagnosis
diag_codes <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", 
                "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
                "V", "W", "X", "Y", "Z")

for (i in 1:length(time_list)) {
urines_ref <- urines_ref %>% prev_ICD_applier(diagnoses,glue("pDIAG{time_list[i]}_"),diag_codes,time_frame = time_list[i])
ur_util <- ur_util %>% prev_ICD_applier(diagnoses,glue("pDIAG{time_list[i]}_"),diag_codes,time_frame = time_list[i])
}

new_icd_cols <- urines_ref %>% select(pDIAG30_A:pDIAG10000_Z)
urines5_combined <- data.frame(cbind(urines5_combined,new_icd_cols))

###Tuning coded ICD-10 procedure
proc_codes <- c("0", "3", "8", "5", "T", "4", "S", "A", "9", 
                "H", "I", "B", "7", "G", "1", "R", "J", "Q", 
                "K", "6", "M", "P", "L", "D", "F", "2", "N", 
                "C", "E", "X", "O")

for (i in 1:length(time_list)) {
urines_ref <- urines_ref %>% prev_ICD_applier(procedures,glue("pPROC{time_list[i]}_"),proc_codes,time_frame = time_list[i])
ur_util <- ur_util %>% prev_ICD_applier(procedures,glue("pPROC{time_list[i]}_"),proc_codes,time_frame = time_list[i])
}

new_proc_cols <- urines_ref %>% select(pPROC30_0:pPROC10000_O)
urines5_combined <- data.frame(cbind(urines5_combined,new_proc_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###Tuning recent organism growth
urine_df <- urine_df %>% mutate(admittime=charttime)
organisms <- urine_df %>% count(org_fullname) %>% arrange(desc(n)) %>% 
  dplyr::slice(1:10) %>% pull(org_fullname)

for (j in 1:length(time_list)) {
urines_ref <- reduce(seq_along(organisms), function(df, i) {
  apply_prev_event(df, paste0("pG", organisms[i],"Urine",as.character(time_list[j])), organisms[i],time_frame=time_list[j])
}, .init = urines_ref) %>%
  ungroup()
ur_util <- reduce(seq_along(organisms), function(df, i) {
  apply_prev_event(df, paste0("pG", organisms[i],"Urine",as.character(time_list[j])), organisms[i],time_frame=time_list[j])
}, .init = ur_util) %>%
  ungroup()
}

new_org_cols <- urines_ref %>% select(`pGEscherichia coliUrine30`:`pGMorganella morganiiUrine10000`)
urines5_combined <- data.frame(cbind(urines5_combined,new_org_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")

###Tuning other previous care events
time_list1 <- c(7,365,720,1e4)
time_list2 <- c(7,30,720,1e4)

for (i in 1:length(time_list1)) {

urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"DNR",field_value,pDNR2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pDNR_{time_list2[i]}")
urines_ref <- urines_ref %>%   
  care_event_assigner(poe,"Psychiatry",field_value,pPsych2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pPsych_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Nephrostomy",field_value,pNeph2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pNeph_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Surgery",field_value,pSURG2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pSURG_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pNUTR_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pPhysio_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Restraints",order_subtype,pRestr2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pRestr_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pOT_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Central TPN",order_subtype,pTPN2,"ordertime",time_list2[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pTPN_{time_list2[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"cath",field_value,pCATH2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pCATH_{time_list1[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Discharge",field_value,pDISC2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pDISC_{time_list1[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"ICU",field_value,pICU2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pICU_{time_list1[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"NGT",field_value,pNGT2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pNGT_{time_list1[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Hydration",field_value,pHyd2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pHyd_{time_list1[i]}")
urines_ref <- urines_ref %>% 
  care_event_assigner(poe,"Chemo",field_value,pChemo2,"ordertime",time_list1[i])
colnames(urines_ref)[ncol(urines_ref)] <- glue("pChemo_{time_list1[i]}")

ur_util <- ur_util %>% 
  care_event_assigner(poe,"DNR",field_value,pDNR2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pDNR_{time_list2[i]}")
ur_util <- ur_util %>%   
  care_event_assigner(poe,"Psychiatry",field_value,pPsych2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pPsych_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Nephrostomy",field_value,pNeph2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pNeph_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Surgery",field_value,pSURG2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pSURG_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Nutrition consult",order_subtype,pNUTR2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pNUTR_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Physical Therapy",order_subtype,pPhysio2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pPhysio_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Restraints",order_subtype,pRestr2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pRestr_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Occupational Therapy",order_subtype,pOT2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pOT_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Central TPN",order_subtype,pTPN2,"ordertime",time_list2[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pTPN_{time_list2[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"cath",field_value,pCATH2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pCATH_{time_list1[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Discharge",field_value,pDISC2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pDISC_{time_list1[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"ICU",field_value,pICU2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pICU_{time_list1[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"NGT",field_value,pNGT2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pNGT_{time_list1[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Hydration",field_value,pHyd2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pHyd_{time_list1[i]}")
ur_util <- ur_util %>% 
  care_event_assigner(poe,"Chemo",field_value,pChemo2,"ordertime",time_list1[i])
colnames(ur_util)[ncol(ur_util)] <- glue("pChemo_{time_list1[i]}")

}

new_poe_cols <- urines_ref %>% select(pDNR_7:pChemo_10000)
urines5_combined <- data.frame(cbind(urines5_combined,new_poe_cols))

write_csv(urines5_combined,"urines5_combined.csv")
write_csv(ur_util,"ur_util_combined.csv")
write_csv(urines_ref,"urines_ref_combined.csv")

###Final model

test_probs_df <- data.frame(matrix(nrow=floor(nrow(urines5_combined)*0.2),ncol=0))
micro_probs_df <- data.frame(matrix(nrow=nrow(ur_xg_combined),ncol=0))
aucs <- data.frame(matrix(nrow=1,ncol=0))
shap_summary_tables <- list()
metrics_list <- list() 

for (outcome in colnames(urines5_outcomes)) {
  
  if (sum(!is.na(urines5_combined[[outcome]])) > 0) {
    
    trainIndex <- createDataPartition(urines5_combined[[outcome]], p = .8, list = FALSE, times = 1)
    urines5Train <- urines5_combined[trainIndex, ]
    urines5Test <- urines5_combined[-trainIndex, ]
    
    predictor_columns <- colnames(urines5_predictors)
    selected_columns <- intersect(predictor_columns, colnames(urines5Train))
    missing_cols <- setdiff(selected_columns, colnames(ur_xg_combined))
    ur_xg_combined[missing_cols] <- 0
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(all_of(selected_columns))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(all_of(selected_columns))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(all_of(selected_columns))), 
                                label = ur_xg_combined[[outcome]])

    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = final_bestparams[[outcome]]$eta,
      max_depth = final_bestparams[[outcome]]$max_depth,
      min_child_weight = final_bestparams[[outcome]]$min_child_weight,
      subsample = final_bestparams[[outcome]]$subsample,
      colsample_bytree = final_bestparams[[outcome]]$colsample_bytree
    )
    
    print("Running CV...")
    
    cv_model <- xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = 1000,
      nfold = 5,
      early_stopping_rounds = 50,
      verbose = 1,
    )
    
    optimal_nrounds <- cv_model$best_iteration
    
    print("Training...")
    
    xgb_model <- xgb.train(
      params = params,
      data = train_matrix,
      nrounds = optimal_nrounds,
    )
    
    print("Shapping...")
    
    shap_values <- predict(xgb_model, newdata = train_matrix, predcontrib = TRUE)
    shap_df <- as.data.frame(shap_values)
    shap_df <- shap_df[, -ncol(shap_df)]
    shap_summary <- data.frame(
      Feature = colnames(shap_df),
      MeanAbsSHAP = colMeans(abs(shap_df))
    )
    shap_summary <- shap_summary %>% filter(MeanAbsSHAP!=0)
    
    #Feature selection
    train_matrix <- xgb.DMatrix(data = as.matrix(urines5Train %>% select(shap_summary %>% pull(Feature))), 
                                label = urines5Train[[outcome]])
    test_matrix <- xgb.DMatrix(data = as.matrix(urines5Test %>% select(shap_summary %>% pull(Feature))), 
                               label = urines5Test[[outcome]])
    micro_matrix <- xgb.DMatrix(data = as.matrix(ur_xg_combined %>% select(shap_summary %>% pull(Feature))), 
                                label = ur_xg_combined[[outcome]])
    
    #Run again with selected features
    xgb_model <- xgb.train(
      params = params,
      data = train_matrix,
      nrounds = optimal_nrounds,
    )
    
    pred_prob_test <- predict(xgb_model, newdata = test_matrix)
    roc_result <- roc(urines5Test[[outcome]], pred_prob_test)
    auc_value <- auc(roc_result)
    print(paste("AUC-ROC:", auc_value))
    aucs[[outcome]] <- auc_value
    pdf(glue("{outcome}_xg_roc.pdf"), width = 10, height = 10)
    plot(roc_result, main = glue("{outcome} ROC Curve"), col = "blue")
    dev.off()
    pred_prob_micro <- predict(xgb_model, newdata = micro_matrix)
    
    test_probs_df[[outcome]] <- pred_prob_test
    micro_probs_df[[outcome]] <- pred_prob_micro
    
    shap_summary <- shap_summary[order(-shap_summary$MeanAbsSHAP), ]
    
    shap_summary_tables[[outcome]] <- shap_summary
    
    pred_test_class <- ifelse(pred_prob_test > 0.5, 1, 0)
    actual_test_class <- urines5Test[[outcome]]
    
    confusion <- confusionMatrix(factor(pred_test_class), factor(actual_test_class))
    accuracy <- confusion$overall['Accuracy']
    precision <- confusion$byClass['Precision']
    recall <- confusion$byClass['Recall']
    f1_score <- 2 * (precision * recall) / (precision + recall) # F1 calculation
    
    metrics_list[[outcome]] <- list(
      AUC = auc_value,
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1_Score = f1_score
    )
    
  }
}

for (i in 1:length(shap_summary_tables)) {
  
  shappy <- data.frame(shap_summary_tables[i]) %>% mutate(across(2, ~ round(., 3))) %>% filter(if_any(2, ~ . != 0))
  
  colnames(shappy) <- c("Feature","Variable")
  
  write_csv(shappy,glue("SHAP_{combined_antimicrobial_map[i]}.csv"))
  
}

for (i in 1:length(metrics_list)) {
  
  metricky <- data.frame(metrics_list[i])
  
  write_csv(metricky,glue("metrics_{combined_antimicrobial_map[i]}.csv"))
  
}

write_csv(micro_probs_df,"micro_xg_probs_df.csv")
write_csv(test_probs_df,"test_xg_probs_df.csv")
write_csv(aucs,"xg_aucs.csv")

problist <- micro_probs_df %>% melt()

util_xg_df <- util_probs_df %>% 
  mutate(S=problist$value, R=1-problist$value)

weightseq <- seq(0,20,1)
iv_perclist <- c()
po_perclist <- c()
overall_perclist <- c()
ivac_list <- c()
poac_list <- c()
ovac_list <- c()
ovor_list <- c()
oviv_list <- c()

for(weight in seq_along(weightseq)) {
  
  res_sens_analysis(ur_util,util_xg_df,NEWS_variable=weightseq[weight],R_value=1)
  
  iv_perclist <- c(iv_perclist,iv_perc)
  po_perclist <- c(po_perclist,po_perc)
  overall_perclist <- c(overall_perclist,overall_perc)
  ivac_list <- c(ivac_list,iv_s_access)
  poac_list <- c(poac_list,po_s_access)
  ovac_list <- c(ovac_list,overall_s_access)
  ovor_list <- c(ovor_list,overall_s_oral)
  oviv_list <- c(oviv_list,overall_s_iv)
  
}

label_binder <- function(vec,label) {
  
  data.frame(vec,Metric=label,Weight=weightseq)
  
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
iv_perclist <- iv_perclist %>% mutate(Percentage=100-Percentage)
po_perclist <- po_perclist %>% mutate(Percentage=100-Percentage)
overall_perclist <- overall_perclist %>% mutate(Percentage=100-Percentage)

iv_xg_plot_df <- data.frame(rbind(
  iv_perclist,ivac_list
))
po_xg_plot_df <- data.frame(rbind(
  po_perclist,poac_list
))
overall_xg_plot_df <- data.frame(rbind(
  overall_perclist,ovac_list,ovor_list,oviv_list
))

write_csv(iv_xg_plot_df,"iv_xg_plot_df.csv")
write_csv(po_xg_plot_df,"po_xg_plot_df.csv")
write_csv(overall_xg_plot_df,"overall_xg_plot_df.csv")

iv_xg_plot_df %>% susc_plotter("IV ","NEWS",agent_col1=TZP,agent_name1="Piperacillin-tazobactam")
po_xg_plot_df %>% susc_plotter("oral ","NEWS",agent_col1=NIT,agent_name1="Nitrofurantoin")
overall_xg_plot_nopo <- overall_xg_plot_df %>% filter(Metric!="Oral agents")
overall_xg_plot_nopo %>% susc_plotter_overall("overall ", "NEWS",agent_col1=NIT,agent_name1="Nitrofurantoin",
                                    agent_col2=TZP,agent_name2="Piperacillin-tazobactam")

###Combination R weighting sensitivity analysis
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

weightseq <- seq(0,30,5)
iv_perclist_combo <- c()
po_perclist_combo <- c()
overall_perclist_combo <- c()
ivac_list_combo <- c()
poac_list_combo <- c()
ovac_list_combo <- c()
oviv_list_combo <- c()
ovor_list_combo <- c()


for(weight in seq_along(weightseq)) {
  
  combo_res_sens_analysis(ur_util,util_probs_df,NEWS_variable=weightseq[weight],R_value=1)
  
  iv_perclist_combo <- c(iv_perclist_combo,iv_perc)
  po_perclist_combo <- c(po_perclist_combo,po_perc)
  overall_perclist_combo <- c(overall_perclist_combo,overall_perc)
  ivac_list_combo <- c(ivac_list_combo,iv_s_access)
  poac_list_combo <- c(poac_list_combo,po_s_access)
  ovac_list_combo <- c(ovac_list_combo,overall_s_access)
  ovor_list_combo <- c(ovor_list_combo,overall_s_oral)
  oviv_list_combo <- c(oviv_list_combo,overall_s_iv)
  
}

iv_perclist_combo <- iv_perclist_combo %>% label_binder("All agents")
po_perclist_combo <- po_perclist_combo %>% label_binder("All agents")
overall_perclist_combo <- overall_perclist_combo %>% label_binder("All agents")
ivac_list_combo <- ivac_list_combo %>% label_binder("Access agents")
poac_list_combo <- poac_list_combo %>% label_binder("Access agents")
ovac_list_combo <- ovac_list_combo %>% label_binder("Access agents")
ovor_list_combo <- ovor_list_combo %>% label_binder("Oral agents")
oviv_list_combo <- oviv_list_combo %>% label_binder("IV agents")
iv_perclist_combo <- iv_perclist_combo %>% rename(Percentage = "vec")
po_perclist_combo <- po_perclist_combo %>% rename(Percentage = "vec")
overall_perclist_combo <- overall_perclist_combo %>% rename(Percentage = "vec")
ivac_list_combo <- ivac_list_combo %>% rename(Percentage = "vec")
poac_list_combo <- poac_list_combo %>% rename(Percentage = "vec")
ovac_list_combo <- ovac_list_combo %>% rename(Percentage = "vec")
ovor_list_combo <- ovor_list_combo %>% rename(Percentage = "vec")
oviv_list_combo <- oviv_list_combo %>% rename(Percentage = "vec")
iv_perclist_combo <- iv_perclist_combo %>% mutate(Percentage=100-Percentage)
po_perclist_combo <- po_perclist_combo %>% mutate(Percentage=100-Percentage)
overall_perclist_combo <- overall_perclist_combo %>% mutate(Percentage=100-Percentage)

iv_plot_df_combo <- data.frame(rbind(
  iv_perclist_combo,ivac_list_combo
))
po_plot_df_combo <- data.frame(rbind(
  po_perclist_combo,poac_list_combo
))
overall_plot_df_combo <- data.frame(rbind(
  overall_perclist_combo,ovac_list_combo
))

write_csv(iv_plot_df_combo,"iv_plot_df_combo.csv")
write_csv(po_plot_df_combo,"po_plot_df_combo.csv")
write_csv(overall_plot_df_combo,"overall_plot_df_combo.csv")

iv_plot_df_combo %>% susc_plotter("IV ","Resistance","(including combinations)",
                                     agent_col1=TZP,agent_name1="Piperacillin-tazobactam")
po_plot_df_combo %>% susc_plotter("oral ","Resistance","(including combinations)",
                                  agent_col1=NIT,agent_name1="Nitrofurantoin")
overall_plot_df_combo_nopo <- overall_xg_plot_df %>% filter(Metric!="Oral agents")
overall_plot_df_combo_nopo %>% susc_plotter_overall("overall (no oral) ","Resistance","(including combinations)",
                                               agent_col1=NIT,agent_name1="Nitrofurantoin",
                                               agent_col2=TZP,agent_name2="Piperacillin-tazobactam")


###R weighting analysis with improved probability predictions
result_key <- ur_util %>% select(micro_specimen_id,AMP:VAN,AMP_SAM:NIT_VAN)
util_probs_df <- util_probs_df %>% left_join(result_key)
util_probs_df <- util_probs_df %>% mutate(better_R=case_when(Antimicrobial=="Ampicillin"&AMP=="S"~R/2,
                                            Antimicrobial=="Ampicillin"&AMP=="R"~(R+1)/2,Antimicrobial=="Ampicillin"&AMP=="S"~R/2,
                                            Antimicrobial=="Ampicillin-sulbactam"&SAM=="S"~R/2,
                                            Antimicrobial=="Ampicillin-sulbactam"&SAM=="R"~(R+1)/2,
                                            Antimicrobial=="Piperacillin-tazobactam"&TZP=="S"~R/2,
                                            Antimicrobial=="Piperacillin-tazobactam"&TZP=="R"~(R+1)/2,
                                            Antimicrobial=="Cefazolin"&CZO=="S"~R/2,
                                            Antimicrobial=="Cefazolin"&CZO=="R"~(R+1)/2,
                                            Antimicrobial=="Ceftriaxone"&CRO=="S"~R/2,
                                            Antimicrobial=="Ceftriaxone"&CRO=="R"~(R+1)/2,
                                            Antimicrobial=="Ceftazidime"&CAZ=="S"~R/2,
                                            Antimicrobial=="Ceftazidime"&CAZ=="R"~(R+1)/2,
                                            Antimicrobial=="Cefepime"&FEP=="S"~R/2,
                                            Antimicrobial=="Cefepime"&FEP=="R"~(R+1)/2,
                                            Antimicrobial=="Meropenem"&MEM=="S"~R/2,
                                            Antimicrobial=="Meropenem"&MEM=="R"~(R+1)/2,
                                            Antimicrobial=="Ciprofloxacin"&CIP=="S"~R/2,
                                            Antimicrobial=="Ciprofloxacin"&CIP=="R"~(R+1)/2,
                                            Antimicrobial=="Gentamicin"&GEN=="S"~R/2,
                                            Antimicrobial=="Gentamicin"&GEN=="R"~(R+1)/2,
                                            Antimicrobial=="Trimethoprim-sulfamethoxazole"&SXT=="S"~R/2,
                                            Antimicrobial=="Trimethoprim-sulfamethoxazole"&SXT=="R"~(R+1)/2,
                                            Antimicrobial=="Nitrofurantoin"&NIT=="S"~R/2,
                                            Antimicrobial=="Nitrofurantoin"&NIT=="R"~(R+1)/2,
                                            Antimicrobial=="Vancomycin"&VAN=="S"~R/2,
                                            Antimicrobial=="Vancomycin"&VAN=="R"~(R+1)/2,
                                            TRUE~NA))

weightseq <- seq(0,30,5)
iv_perclist_imp <- c()
po_perclist_imp <- c()
overall_perclist_imp <- c()
ivac_list_imp <- c()
poac_list_imp <- c()
ovac_list_imp <- c()

for(weight in seq_along(weightseq)) {
  
  res_sens_analysis_2(ur_util,util_probs_df,weightseq[weight])
  
  iv_perclist_imp <- c(iv_perclist_imp,iv_perc)
  po_perclist_imp <- c(po_perclist_imp,po_perc)
  overall_perclist_imp <- c(overall_perclist_imp,overall_perc)
  ivac_list_imp <- c(ivac_list_imp,iv_s_access)
  poac_list_imp <- c(poac_list_imp,po_s_access)
  ovac_list_imp <- c(ovac_list_imp,overall_s_access)
  
}

label_binder <- function(vec,label) {
  
  data.frame(vec,Metric=label,Weight=weightseq)
  
}

iv_perclist_imp <- iv_perclist_imp %>% label_binder("All agents")
po_perclist_imp <- po_perclist_imp %>% label_binder("All agents")
overall_perclist_imp <- overall_perclist_imp %>% label_binder("All agents")
ivac_list_imp <- ivac_list_imp %>% label_binder("Access agents")
poac_list_imp <- poac_list_imp %>% label_binder("Access agents")
ovac_list_imp <- ovac_list_imp %>% label_binder("Access agents")
iv_perclist_imp <- iv_perclist_imp %>% rename(Percentage = "vec")
po_perclist_imp <- po_perclist_imp %>% rename(Percentage = "vec")
overall_perclist_imp <- overall_perclist_imp %>% rename(Percentage = "vec")
ivac_list_imp <- ivac_list_imp %>% rename(Percentage = "vec")
poac_list_imp <- poac_list_imp %>% rename(Percentage = "vec")
ovac_list_imp <- ovac_list_imp %>% rename(Percentage = "vec")
iv_perclist_imp <- iv_perclist_imp %>% mutate(Percentage=100-Percentage)
po_perclist_imp <- po_perclist_imp %>% mutate(Percentage=100-Percentage)
overall_perclist_imp <- overall_perclist_imp %>% mutate(Percentage=100-Percentage)

iv_plot_df_imp <- data.frame(rbind(
  iv_perclist_imp,ivac_list_imp
))
po_plot_df_imp <- data.frame(rbind(
  po_perclist_imp,poac_list_imp
))
overall_plot_df_imp <- data.frame(rbind(
  overall_perclist_imp,ovac_list_imp
))

write_csv(iv_plot_df_imp,"iv_plot_df_imp.csv")
write_csv(po_plot_df_imp,"po_plot_df_imp.csv")
write_csv(overall_plot_df_imp,"overall_plot_df_imp.csv")

iv_plot_df_imp %>% susc_plotter_iv("IV ","Resistance","(prediction accuracy improved)",agent_col1=GEN,agent_name1="Gentamicin",
                                agent_col2=TZP,agent_name2="Piperacillin-tazobactam")
po_plot_df_imp %>% susc_plotter("oral ","Resistance","(prediction accuracy improved)",agent_col1=NIT,agent_name1="Nitrofurantoin")
overall_plot_df_imp %>% susc_plotter("overall ", "Resistance","(prediction accuracy improved)",agent_col1=NIT,agent_name1="Nitrofurantoin")

###2nd-line oral R weight sensitivity analysis
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
                 "MEM","CIP","GEN","SXT")
ab_singles <- all_singles
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
                "GEN","SXT")
iv_ab_singles <- iv_singles
iv_combos <- combn(iv_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_ivs <- c(iv_singles, iv_combos)
oral_singles <- c("AMP","SAM","CIP",
                  "SXT")
oral_ab_singles <- oral_singles
oral_combos <- combn(oral_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_orals <- c(oral_singles, oral_combos)
access_singles <- c("AMP","SAM","GEN",
                    "SXT","CZO")
access_combos <- combn(access_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_access <- c(access_singles, access_combos)
watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
watch_combos <- combn(watch_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_watch <- c(watch_singles, watch_combos)

weightseq <- seq(0,30,5)
iv_perclist_2nd <- c()
po_perclist_2nd <- c()
overall_perclist_2nd <- c()
ivac_list_2nd <- c()
poac_list_2nd <- c()
ovac_list_2nd <- c()

for(weight in seq_along(weightseq)) {
  
  res_sens_analysis(ur_util,util_probs_df,weightseq[weight])
  
  iv_perclist_2nd <- c(iv_perclist_2nd,iv_perc)
  po_perclist_2nd <- c(po_perclist_2nd,po_perc)
  overall_perclist_2nd <- c(overall_perclist_2nd,overall_perc)
  ivac_list_2nd <- c(ivac_list_2nd,iv_s_access)
  poac_list_2nd <- c(poac_list_2nd,po_s_access)
  ovac_list_2nd <- c(ovac_list_2nd,overall_s_access)
  
}

iv_perclist_2nd <- iv_perclist_2nd %>% label_binder("All agents")
po_perclist_2nd <- po_perclist_2nd %>% label_binder("All agents")
overall_perclist_2nd <- overall_perclist_2nd %>% label_binder("All agents")
ivac_list_2nd <- ivac_list_2nd %>% label_binder("Access agents")
poac_list_2nd <- poac_list_2nd %>% label_binder("Access agents")
ovac_list_2nd <- ovac_list_2nd %>% label_binder("Access agents")
iv_perclist_2nd <- iv_perclist_2nd %>% rename(Percentage = "vec")
po_perclist_2nd <- po_perclist_2nd %>% rename(Percentage = "vec")
overall_perclist_2nd <- overall_perclist_2nd %>% rename(Percentage = "vec")
ivac_list_2nd <- ivac_list_2nd %>% rename(Percentage = "vec")
poac_list_2nd <- poac_list_2nd %>% rename(Percentage = "vec")
ovac_list_2nd <- ovac_list_2nd %>% rename(Percentage = "vec")
iv_perclist_2nd <- iv_perclist_2nd %>% mutate(Percentage=100-Percentage)
po_perclist_2nd <- po_perclist_2nd %>% mutate(Percentage=100-Percentage)
overall_perclist_2nd <- overall_perclist_2nd %>% mutate(Percentage=100-Percentage)

iv_plot_df_2nd <- data.frame(rbind(
  iv_perclist_2nd,ivac_list_2nd
))
po_plot_df_2nd <- data.frame(rbind(
  po_perclist_2nd,poac_list_2nd
))
overall_plot_df_2nd <- data.frame(rbind(
  overall_perclist_2nd,ovac_list_2nd
))

write_csv(iv_plot_df_2nd,"iv_plot_df_2nd.csv")
write_csv(po_plot_df_2nd,"po_plot_df_2nd.csv")
write_csv(overall_plot_df_2nd,"overall_plot_df_2nd.csv")

po_plot_df_2nd %>% susc_plotter("oral ","Resistance"," (Second-line)",agent_col1=SAM,agent_name1="Ampicillin-sulbactam")

