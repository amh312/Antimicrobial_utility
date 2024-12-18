#MICROSIMULATION ANALYSIS

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

##Uploads and reference lists

###Uploads
ur_util <- read_csv("ur_util_final.csv")
urines_abx <- read_csv("ur_util_final.csv")
abx <- read_csv("interim_abx.csv")
urines5_desc <- read_csv("urines5_ref.csv")
pats <- read_csv("patients.csv")

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

##Sample size calculation for chi-squared test

###Pooled proportion
prop_1 <- 0.7 #expected standard
prop_2 <- 0.75 #expected intervention
pooled_prop <- (prop_1+prop_2)/2

###Effect size
ef_size <- sqrt((prop_1-prop_2)^2/(pooled_prop*(1-pooled_prop)))

###Z values
Za <- 1.96 #5% sig level
Zb <- 0.84 #80% power

###Sample size
samp_size <- (Za+Zb)^2/ef_size^2
glue("sample size requirement for chi-squared test = {samp_size}")

##Descriptive data
abx %>% distinct(subject_id,.keep_all = T) %>% count(CDI)

###Table 1: Patient/specimen characteristics
yearkey <- pats %>% distinct(subject_id,.keep_all = T) %>% select(subject_id,anchor_year_group)
agekey <- pats %>% distinct(subject_id,.keep_all = T) %>% select(subject_id,anchor_age) %>% 
  mutate(standard_age = case_when(
    anchor_age<30~18,
    anchor_age>=30&anchor_age<40~30,
    anchor_age>=40&anchor_age<50~40,
    anchor_age>=50&anchor_age<60~50,
    anchor_age>=60&anchor_age<70~60,
    anchor_age>=70&anchor_age<80~70,
    anchor_age>=80&anchor_age<90~80,
    anchor_age>=90~90,
  )) %>% select(-anchor_age)
hadmkey <- hadm %>% distinct(subject_id,.keep_all = T) %>% select(subject_id,race,marital_status,language,insurance)

urines5_desc %>% tab_mutater("Urine model") %>% print(n=100)
tab_mutater <- function(df,dataset,patient_level="Y") {
  
  supercounter_1 <- function(df) {
    
    counter <- function(df,charac) {
      
      charac <- enquo(charac)
      
      df %>% count(!!charac) %>% arrange(desc(n)) %>% 
        mutate(Characteristic=names(.)[1],
               `n (%)`=glue("{n} ({round((n/nrow(df))*100,1)})")) %>% 
        rename(Subtype = 1) %>% select(-n)
      
    }
    
    tibble(rbind(
      df %>% counter(Gender),
      df %>% counter(Race),
      df %>% counter(`Age group`),
      df %>% counter(`Marital status`),
      df %>% counter(`Language spoken`),
      df %>% counter(Insurance),
      df %>% counter(`Year group`),
      df %>% counter(`Ampicillin result`),
      df %>% counter(`Ampicillin/sulbactam result`),
      df %>% counter(`Piperacillin/tazobactam result`),
      df %>% counter(`Cefazolin result`),
      df %>% counter(`Ceftriaxone result`),
      df %>% counter(`Ceftazidime result`),
      df %>% counter(`Cefepime result`),
      df %>% counter(`Meropenem result`),
      df %>% counter(`Ciprofloxacin result`),
      df %>% counter(`Gentamicin result`),
      df %>% counter(`Trimethoprim/sulfamethoxazole result`),
      df %>% counter(`Nitrofurantoin result`),
      df %>% counter(`Vancomycin result`)
    )) %>% relocate(Characteristic,.before = "Subtype")
    
  }

  supercounter_2 <- function(df) {
    
    counter <- function(df,charac) {
      
      charac <- enquo(charac)
      
      df %>% count(!!charac) %>% arrange(desc(n)) %>% 
        mutate(Characteristic=names(.)[1],
               `n (%)`=glue("{n} ({round((n/nrow(df))*100,1)})")) %>% 
        rename(Subtype = 1) %>% select(-n)
      
    }
    
    tibble(rbind(
      df %>% counter(Gender),
      df %>% counter(Race),
      df %>% counter(`Age group`),
      df %>% counter(`Marital status`),
      df %>% counter(`Language spoken`),
      df %>% counter(Insurance),
      df %>% counter(`Year group`),
      df %>% counter(CDI),
      df %>% counter(Toxicity)
    )) %>% relocate(Characteristic,.before = "Subtype")
    
  }
  
  if (patient_level=="Y") {
    
    df <- df
    
  } else {
    
    df <- df %>% distinct(subject_id,.keep_all = T)
    
  }
  
  if ("micro_specimen_id" %in% colnames(df)) {
  
df1 <- df %>% mutate(
  Gender = case_when(MALE ~ "M", TRUE~"F"),
  Race = case_when(
    grepl("WHITE",race)~"White",
    grepl("BLACK",race)~"Black",
    grepl("ASIAN",race)~"Asian",
    grepl("HISPANIC",race)~"Hispanic",
    grepl("UNKNOWN",race)~"Unknown",
    TRUE~"Other"
  ),
  `Marital status` = case_when(
    grepl("^MARRIED$",marital_status)~"Married",
    TRUE~"Unmarried/unknown"
  ),
  `Age group` = case_when(
    grepl("18",standard_age)~"18-29",
    grepl("30",standard_age)~"30-39",
    grepl("40",standard_age)~"40-49",
    grepl("50",standard_age)~"50-59",
    grepl("60",standard_age)~"60-69",
    grepl("70",standard_age)~"70-79",
    grepl("80",standard_age)~"80-89",
    grepl("90",standard_age)~"≥ 90",
  ),
  `Language spoken`=case_when(
    grepl("ENGLISH",language)~"English",
    TRUE~"Other/unknown",
  ),
  Insurance = case_when(
    grepl("(UNKNOWN|Other)",insurance)~"Other/unknown",
    TRUE~insurance
  ),
  `Genus grown`=case_when(grepl("(Raoul|Haf|Staph)",org_genus)~"Other",
                          TRUE~org_genus),
  `Ampicillin result`=AMP,
  `Ampicillin/sulbactam result`=SAM,
  `Piperacillin/tazobactam result`=TZP,
 `Cefazolin result`=CZO,
 `Ceftriaxone result`=CRO,
 `Ceftazidime result`=CAZ,
 `Cefepime result`=FEP,
 `Meropenem result`=MEM,
 `Ciprofloxacin result`=CIP,
 `Gentamicin result`=GEN,
  `Trimethoprim/sulfamethoxazole result`=SXT,
 `Nitrofurantoin result`=NIT,
 `Vancomycin result`=VAN
) %>% left_join(yearkey) %>% 
  rename(`Year group`="anchor_year_group") %>% 
  filter(`Year group`!="2020 - 2022") %>% 
  supercounter_1()

  } else {
  
    df1 <- df %>% left_join(hadmkey) %>% left_join(agekey) %>% 
      mutate(
      Gender = case_when(MALE ~ "M", TRUE~"F"),
      Race = case_when(
        grepl("WHITE",race)~"White",
        grepl("BLACK",race)~"Black",
        grepl("ASIAN",race)~"Asian",
        grepl("HISPANIC",race)~"Hispanic",
        grepl("UNKNOWN",race)~"Unknown",
        TRUE~"Other"
      ),
      `Marital status` = case_when(
        grepl("^MARRIED$",marital_status)~"Married",
        TRUE~"Unmarried/unknown"
      ),
      `Age group` = case_when(
        grepl("18",standard_age)~"18-29",
        grepl("30",standard_age)~"30-39",
        grepl("40",standard_age)~"40-49",
        grepl("50",standard_age)~"50-59",
        grepl("60",standard_age)~"60-69",
        grepl("70",standard_age)~"70-79",
        grepl("80",standard_age)~"80-89",
        grepl("90",standard_age)~"≥ 90",
      ),
      `Language spoken`=case_when(
        grepl("ENGLISH",language)~"English",
        TRUE~"Other/unknown",
      ),
      Insurance = case_when(
        grepl("(UNKNOWN|Other)",insurance)~"Other/unknown",
        TRUE~insurance
      ),
      CDI = case_when(CDI~"CDI",TRUE~NA),
      Toxicity = case_when(overall_tox~"Toxicity",
                           TRUE~NA)
    ) %>% left_join(yearkey) %>% 
      rename(`Year group`="anchor_year_group") %>% 
      filter(`Year group`!="2020 - 2022") %>% 
      supercounter_2() %>% filter(!is.na(Subtype)) %>% 
      mutate(Characteristic=case_when(grepl("(CDI|Toxicity)",Characteristic)~"Antibiotic outcome",
                                      TRUE~Characteristic)
             )
    
}

colnames(df1)[3] <- glue("{dataset} {colnames(df1)[3]}")

return(df1)

}

desc_tab <-  abx %>% tab_mutater("Prescription model","N") %>% 
  full_join(urines5_desc %>% tab_mutater("Urine model")) %>% 
  left_join(ur_util %>% tab_mutater("Urine microsimulation"))
totals <- list("Total","Patients",
               glue("{as.character(nrow(abx %>% distinct(subject_id)))} (100)"),
               glue("{as.character(nrow(urines5_desc))} (100)"),
               glue("{as.character(nrow(ur_util))} (100)")
            )
desc_tab[nrow(desc_tab)+1,] <- totals
desc_tab <- desc_tab %>% mutate(Characteristic=case_when(
  Characteristic==lag(Characteristic)~"",
  TRUE~Characteristic
))

write_csv(desc_tab,"uf_desctab.csv")

###Table 2: Prescription characteristics

ab_tab <- abx %>% mutate(abx_name=str_replace_all(abx_name,"_"," & ")) %>% count(abx_name) %>% 
  arrange(desc(n)) %>% mutate(abx_name=case_when(
    n/nrow(abx)<0.0025~"Other",TRUE~abx_name
  )) %>% group_by(abx_name) %>%
  summarise(n=sum(n)) %>% ungroup() %>% print(n=100) %>% 
  arrange(desc(n)) %>% 
  mutate(n=glue("{n} ({round((n/nrow(abx))*100,1)})")) %>%
  rename(`n (%)`="n",`Antibiotic(s)`="abx_name")

write_csv(ab_tab,"ab_tab.csv")

##AST performance analysis - number of results per panel

###Number of AST results per panel
ur_util <- ur_util %>% rpp_ast()

###Assemble data frame for dot plot data visualisation
acs_df <- ur_util %>% assemble_dotplot_df()

write_csv(acs_df,"uf_ast_sourcedata_aware_dotplot.csv")

###Dot plot of number of all S results and Access S results per panel
main_aware_plot <- acs_df %>% main_dotplotter("PDAST\nall S","Standard\nall S","PDAST\nall Access S","Standard\nall Access S",
                                              "All agents","WHO access agents")

###Dot plot of number of all R results and Access R results per panel
accessr_aware_plot <- acs_df %>% main_dotplotter("PDAST\nall R","Standard\nall R","PDAST\nall Access R","Standard\nall Access R",
                                                 "All agents (R)","WHO access agents (R)","(Access agent resistance)")

##AST performance analysis - number of S results per antimicrobial agent

###Count S results per antimicrobial for personalised approach
pdast_all_singles <- ur_util %>% number_ab_results(PDAST_1,PDAST_6,all_singles,"S","I")

###Count S results per antimicrobial for standard approach
standard_all_singles <- ur_util %>% number_ab_results(STANDARD_AST_1,STANDARD_AST_6,all_singles,"S","I")

###Count R results per antimicrobial for personalised approach
pdast_all_singles_r <- ur_util %>% number_ab_results(PDAST_1,PDAST_6,all_singles,"R","NT")

###Count R results per antimicrobial for standard approach
standard_all_singles_r <- ur_util %>% number_ab_results(STANDARD_AST_1,STANDARD_AST_6,all_singles,"R","NT")

###Cleveland dot plot of number of S results per antimicrobial agent
abs_df <- abs_df_assemble(pdast_all_singles,standard_all_singles)
write_csv(abs_df,"s_uf_sourcedata_abs_cleveplot.csv")
s_results_by_ab <- abs_df %>% cleveland_ab_plot("susceptible")

###Join current antimicrobial prescription to urines dataframe
ab_key <- abx %>% filter(is_abx) %>% 
  select(subject_id,starttime,stoptime)
urines_abx <- urines_abx %>% left_join(ab_key,by="subject_id") %>% 
  mutate(on_ab = case_when(
    storetime > starttime & storetime < stoptime ~ TRUE,
    TRUE ~ FALSE )) %>% filter(on_ab)

###Convert I to S and NT to R
urref <- urines_abx %>% select(AMP:VAN)
urref[urref=="NT"] <- "R"
urref[urref=="I"] <- "S"
urines_abx[,18:30] <- urref

##Analysis of AST results for the antimicrobial agent the patient is prescribed

###Find prescribed antimicrobial-result matches
urines_abx <- urines_abx %>% 
  mutate(on_standard = case_when(as.ab(abx_name)==STANDARD_AST_1 ~ NIT,
                                 as.ab(abx_name)==STANDARD_AST_2 ~ SXT,
                                 as.ab(abx_name)==STANDARD_AST_3 ~ CIP,
                                 as.ab(abx_name)==STANDARD_AST_4 ~ TZP,
                                 as.ab(abx_name)==STANDARD_AST_5 ~ GEN,
                                 as.ab(abx_name)==STANDARD_AST_6 ~ CRO,
                                 TRUE ~ "NT")) %>%
  rowwise() %>%
  mutate(on_PDAST = case_when(
    as.ab(abx_name) == PDAST_1 ~ get(PDAST_1),
    as.ab(abx_name) == PDAST_2 ~ get(PDAST_2),
    as.ab(abx_name) == PDAST_3 ~ get(PDAST_3),
    as.ab(abx_name) == PDAST_4 ~ get(PDAST_4),
    as.ab(abx_name) == PDAST_5 ~ get(PDAST_5),
    as.ab(abx_name) == PDAST_6 ~ get(PDAST_6),
    TRUE ~ "NT"
  )) %>%
  ungroup()

###Prepare data frame for bar plot
ablist <- all_singles %>% ab_name() %>% str_replace("-","/")
urines_abx <- urines_abx %>% filter(abx_name %in% ablist)
abx_graph1 <- urines_abx %>% count(abx_name,on_standard,on_PDAST) %>% 
  filter(!(on_standard=="NT" & on_PDAST=="NT")) %>% 
  mutate(Panel=case_when(on_standard!="NT"~"Standard",TRUE~"PDAST"),
         Result=case_when(Panel=="PDAST"~on_PDAST,TRUE~on_standard))
abx_graph2 <- abx_graph1 %>% filter(on_standard!="NT" & on_PDAST!="NT") %>% 
  mutate(Panel="PDAST",Result=on_PDAST)
abx_graph <- abx_graph1 %>% rbind(abx_graph2)
access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")
watch_abs <- c("TZP","CRO","CAZ",
               "FEP","MEM","CIP","VAN")
axiscols <- ifelse(abx_graph %>% group_by(abx_name) %>% 
                     mutate(n=sum(n)) %>% 
                     arrange(n) %>% ungroup() %>% 
                     distinct(abx_name) %>% unlist() %in% ab_name(access_abs),"seagreen",
                   "darkorange")
abx_graph$abx_name <- factor(abx_graph$abx_name,
                             levels=abx_graph %>% group_by(abx_name) %>% 
                               mutate(n=sum(n)) %>% 
                               arrange(n) %>% ungroup() %>% 
                               distinct(abx_name) %>% unlist())
max_count <- ceiling(abx_graph1 %>% select(-Panel,-Result) %>% 
                       group_by(abx_name) %>% mutate(n=sum(n)) %>% arrange(desc(n)) %>%
                       ungroup() %>% dplyr::slice(1) %>% select(n) %>% unlist() /25) * 25

write_csv(abx_graph,"uf_sourcedata_prescr_abx.csv")

###Data visualisation of antimicrobial-result matches
abx_prescribed <- ggplot(abx_graph, aes(x = abx_name, y = if_else(Panel == "Standard", -n, n),
                                        fill = Result,color=Result)) +
  geom_bar(stat = "identity", width=0.85,position = "stack") +
  scale_y_continuous(
    limits = c(-max_count, max_count), 
    breaks = seq(-max_count, max_count, by = 25), 
    labels = abs(seq(-max_count, max_count, by = 25)) 
  ) +
  coord_flip() +
  labs(y = "Number of results", x = "Antimicrobial agent prescribed",
       fill = "Result") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  geom_hline(yintercept = 0,linetype="solid",color="black") +
  ggtitle("Results provided for the antimicrobial agent prescribed (inpatients)") +
  geom_text(x=13.5,y=-100,label="Standard approach",color="#3C3C3C",size=4) +
  geom_text(x=13.5,y=100,label="Personalised approach",color="#3C3C3C",size=4) +
  theme(plot.title = element_text(size = 16, margin = margin(b = 20)),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(
          colour = axiscols))

ggsave("uf_abx_prescribed.pdf", plot = abx_prescribed, device = "pdf", width = 10, height = 4,
       path="/Users/alexhoward/Documents/Projects/UDAST_code")

###Counting the total number of antimicrobial-result matches per paneL
PDAST_total <- totaller("PDAST")
Standard_total <- totaller("Standard")
glue("The personalised approach provided {PDAST_total} results for antimicrobial
     agents patients were currently prescribed, versus {Standard_total} results
     provided by the standard panel")
ur_util %>% filter(!is.na(abx_name)) %>% nrow()

###Counting the total number of S and R results for prescribed agents per panel by AWaRe category
abx_graph <- abx_graph %>% 
  mutate(AWaRe = case_when(abx_name %in% ab_name(access_abs) ~ "Access",
                           TRUE ~ "Watch"))
awaresum <- abx_graph %>% group_by(Panel,Result,AWaRe) %>% 
  summarise(n_results=sum(n)) %>% ungroup()
print(awaresum)
PDAST_S <- awaretotaller("PDAST")
Standard_S <- awaretotaller(("Standard"))
glue("The personalised approach provided {PDAST_S} 'S' results for Access category
      antimicrobial agents patients were currently prescribed, versus {Standard_S} 
      comparable results provided by the standard panel")

###Counting cases where on Watch abx and an Access S alternative / IVOST was provided
urines_abx <- urines_abx %>% rowwise() %>%
  mutate(on_PDAST1 = case_when(
    PDAST_1 %in% access_abs ~ get(PDAST_1),TRUE~"NT"),
    on_PDAST2 = case_when(
      PDAST_2 %in% access_abs ~ get(PDAST_2),TRUE~"NT"),
    on_PDAST3 = case_when(
      PDAST_3 %in% access_abs ~ get(PDAST_3),TRUE~"NT"),
    on_PDAST4 = case_when(
      PDAST_4 %in% access_abs ~ get(PDAST_4),TRUE~"NT"),
    on_PDAST5 = case_when(
      PDAST_5 %in% access_abs ~ get(PDAST_5),TRUE~"NT"),
    on_PDAST6 = case_when(
      PDAST_6 %in% access_abs ~ get(PDAST_6),TRUE~"NT")
  ) %>%
  ungroup() %>% mutate(
    PDAST_S_in_range = apply(select(., "on_PDAST1":"on_PDAST6"), 1, function(row) {
      any(row == "S")
    }),
    Standard_S_in_range = case_when(GEN=="S"|NIT=="S"|
                                      SXT=="S" ~ TRUE,
                                    TRUE ~ FALSE),
    PDAST_ac_alt = case_when(PDAST_S_in_range & abx_name %in% ab_name(watch_abs) ~ TRUE,
                             TRUE ~ FALSE),
    Standard_ac_alt = case_when(Standard_S_in_range & abx_name %in% ab_name(watch_abs) ~ TRUE,
                                TRUE ~ FALSE))
n_ac_switch_pdast <- sum(urines_abx$PDAST_ac_alt)
n_ac_switch_standard <- sum(urines_abx$Standard_ac_alt)
ivs_only <- c("CRO","TZP","FEP","MEM","CAZ","VAN")
urines_abx <- urines_abx %>% mutate(PDAST_gen_out = case_when(
  PDAST_ac_alt & NIT!="S"&SXT!="S"&AMP!="S"&SAM!="S"&CZO!="S" ~FALSE,
  TRUE~TRUE
),
PDAST_ivost = case_when(abx_name %in% ab_name(ivs_only) &
                          PDAST_ac_alt & NIT=="S"|SXT=="S"|AMP=="S"|SAM=="S" ~TRUE,
                        TRUE~FALSE
),
Standard_ivost = case_when( abx_name %in% ab_name(ivs_only) &
                              Standard_ac_alt & NIT=="S"|SXT=="S"|AMP=="S"|SAM=="S" ~TRUE,
                            TRUE~FALSE
))

n_pdast_ivost <- sum(urines_abx$PDAST_ivost)
n_standard_ivost <- sum(urines_abx$Standard_ivost)

glue("The personalised approach provided an opportunity to switch from a
     prescribed Watch category agent to an Access category agent in {n_ac_switch_pdast}
     instances, compared to {n_ac_switch_standard} instances using the standard approach.
     
     The personalised approach provided an step-down from an IV agent to an Access
     category oral agent in {n_pdast_ivost} cases compared to {n_standard_ivost} cases
     using the standard approach
     
     ")

###Susceptible and resistant results for the agent recommended for treatment
urines_abx <- ur_util
urref <- urines_abx %>% select(AMP:VAN)
urref[urref=="NT"] <- "R"
urref[urref=="I"] <- "S"
urines_abx[,18:30] <- urref
urines_abx <- urines_abx %>%
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

###Data visualisation of antimicrobial-result matches

abx_graph <- urines_abx %>% ab_rec_df(Intravenous_1,"Intravenous_1",PDIVRx_1_result,STIVRx_1_result,
                                      urines_abx %>% distinct(STANDARD_IV_1) %>% ab_name())
write_csv(abx_graph,"uf_sourcedata_prescr_abx.csv")
abx_graph %>% ab_result_graph("first","intravenous",9,900)

abx_graph_2 <- urines_abx %>% ab_rec_df(Intravenous_2,"Intravenous_2",PDIVRx_2_result,STIVRx_2_result,
                                        urines_abx %>% distinct(STANDARD_IV_2) %>% ab_name())
write_csv(abx_graph_2,"uf_sourcedata_prescr_abx_2.csv")
abx_graph_2 %>% ab_result_graph("second","intravenous",11,900)

abx_graph_3 <- urines_abx %>% ab_rec_df(Oral_1,"Oral_1",PDPORx_1_result,STPORx_1_result,
                                        urines_abx %>% distinct(STANDARD_PO_1) %>% ab_name())
write_csv(abx_graph_3,"uf_sourcedata_prescr_abx_3.csv")
abx_graph_3 %>% ab_result_graph("first","oral",5.65,900)

abx_graph_4 <- urines_abx %>% ab_rec_df(Oral_2,"Oral_2",PDIVRx_2_result,STPORx_2_result,
                                        urines_abx %>% distinct(STANDARD_PO_2) %>% ab_name())
write_csv(abx_graph_4,"uf_sourcedata_prescr_abx_4.csv")
abx_graph_4 %>% ab_result_graph("second","oral",5.65,1000)

###Comparison of antimicrobial agents prescribed and recommended
urines_abx <- urines_abx %>% filter(!grepl("_",abx_name))
urines_abx <- urines_abx %>% mutate(
  Wa_Ac_S_IV = case_when(
  PDIVRx_1_result=="S"&Intravenous_1%in%(access_abs)&
    abx_name%in%(watch_abs %>% ab_name()) ~ TRUE,TRUE~FALSE),
  Wa_Ac_S_PO = case_when(
    PDPORx_1_result=="S"&Oral_1%in%(access_abs)&
      abx_name%in%(watch_abs %>% ab_name()) ~ TRUE,TRUE~FALSE),
  Ac_Wa_IV = case_when(
    Intravenous_1%in%(watch_abs)&
      abx_name%in%(access_abs %>% ab_name()) ~ TRUE,TRUE~FALSE),
  Ac_Wa_PO = case_when(
    Oral_1%in%(watch_abs)&
      abx_name%in%(access_abs %>% ab_name()) ~ TRUE,TRUE~FALSE)
)

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
         STPORx_2_result = get(STANDARD_PO_2)) %>%
  ungroup()

all_res <- nrow(ur_util %>% filter(PDRx_1_result=='R'))
iv_res <- nrow(ur_util %>% filter(PDIVRx_1_result=='R'))
po_res <- nrow(ur_util %>% filter(PDPORx_1_result=='R'))

glue("Where an IV agent was required, the personalised approach recommended a switch from a
     prescribed Watch category agent to an effective Access category agent in 
     {round((sum(urines_abx$Wa_Ac_S_IV)/nrow(urines_abx))*100,1)}% (n={sum(urines_abx$Wa_Ac_S_IV)}) of prescriptions.
     
     Where an oral agent was required, the personalised approach recommended a switch from a
     prescribed Watch category agent to an effective Access category agent in 
     {round((sum(urines_abx$Wa_Ac_S_PO)/nrow(urines_abx))*100,1)}% (n={sum(urines_abx$Wa_Ac_S_PO)}) of prescriptions.
     
     Overall, the personalised approach recommended an agent with a subsequent resistant result in
     {round(all_res/nrow(ur_util)*100,1)}% (n={all_res}) of cases where urine culture was sent.
     
     Where an IV agent was required, the personalised approach recommended an agent with a subsequent resistant result in
     {round(iv_res/nrow(ur_util)*100,1)}% (n={iv_res}) of cases where urine culture was sent.
     
     Where an oral agent was required, the personalised approach recommended an agent with a subsequent resistant result in
     {round(po_res/nrow(ur_util)*100,1)}% (n={po_res}) of cases where urine culture was sent.
     
     ")

ur_util %>% filter(PDRx_1_result=="S") %>% nrow()/nrow(ur_util)
urines_abx %>% filter(PDPORx_1_result=="R") %>% filter(org_fullname=="Escherichia coli") %>% count(Oral_1) %>% arrange(desc(n))
urines_abx %>% filter(PDIVRx_1_result=="R") %>% count(Intravenous_1) %>% arrange(desc(n))
ur_util %>% filter(PDRx_1%in%access_singles) %>% nrow()/nrow(ur_util)
write_csv(ur_util,"ur_util_microsim.csv")

