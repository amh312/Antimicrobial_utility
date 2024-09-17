#MICROSIMULATION ANALYSIS

###Counting the number of S or I results per antimicrobial per panel
number_results_per_panel <- function(df,start_col,end_col,which_abs,result_1,result_2) {

  start_col <- enquo(start_col)
  end_col <- enquo(end_col)
  
  n_all_s <- c()
  
  for(i in 1:nrow(df)) {
    
    n_ac_s <- sum((df %>%
                     select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%slice(i) %>% unlist(),which_abs))) %>% 
                     slice(i)) == result_1 |
                    (df %>%
                       select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%slice(i) %>% unlist(),which_abs))) %>% 
                       slice(i)) == result_2)
    
    n_all_s <- n_all_s %>% append(n_ac_s)
    
  }
  
  n_all_s
  
}

###Counting the number of S or I results per antimicrobial (total)
number_ab_results <- function(df,start_col,end_col,which_abs,result_1,result_2) {
  
  all_si <- c()
  start_col <- enquo(start_col)
  end_col <- enquo(end_col)
  
  for(i in 1:nrow(df)) {
    
    all_s <- df %>%
      select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. ==result_1) %>% rownames()
    
    all_i <- df %>%
      select(all_of(intersect(df %>% select(!!start_col:!!end_col) %>%slice(i) %>% unlist(),all_abs))) %>% 
      slice(i) %>% t() %>% data.frame() %>% filter(. ==result_2) %>% rownames()
    
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

##Uploads and reference lists

###Uploads
ur_util <- read_csv("ur_util_final.csv")

###Reference lists
all_singles <- c("AMP","SAM","TZP","CZO","CRO","CAZ","FEP",
             "MEM","CIP","GEN","SXT","NIT","VAN")
all_combos <- combn(all_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_abs <- c(all_singles,all_combos)
iv_singles <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
            "GEN","SXT","VAN")
iv_combos <- combn(iv_abs, 2, FUN = function(x) paste(x, collapse = "_"))
all_ivs <- c(iv_singles, iv_combos)
oral_singles <- c("AMP","SAM","CIP",
              "SXT","NIT")
oral_combos <- combn(oral_abs, 2, FUN = function(x) paste(x, collapse = "_"))
all_orals <- c(oral_singles, oral_combos)
access_singles <- c("AMP","SAM","GEN",
                  "SXT","NIT","CZO")
access_combos <- combn(access_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_access <- c(access_singles, access_combos)
watch_singles <- c("CRO","CAZ","FEP","MEM","TZP","CIP","VAN")
watch_combos <- combn(watch_singles, 2, FUN = function(x) paste(x, collapse = "_"))
all_watch <- c(watch_singles, watch_combos)

##Calculations

###PDAST results per panel
ur_util$PDAST_rpp_ass <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_singles,"S","I")
ur_util$PDAST_rpp_acs <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_combos,"S","I")
ur_util$PDAST_rpp_aas <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_abs,"S","I")
ur_util$PDAST_rpp_iss <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,iv_singles,"S","I")
ur_util$PDAST_rpp_ics <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,iv_combos,"S","I")
ur_util$PDAST_rpp_ais <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_ivs,"S","I")
ur_util$PDAST_rpp_oss <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,oral_singles,"S","I")
ur_util$PDAST_rpp_ocs <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,oral_combos,"S","I")
ur_util$PDAST_rpp_aos <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_orals,"S","I")
ur_util$PDAST_rpp_acss <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,access_singles,"S","I")
ur_util$PDAST_rpp_accs <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,access_combos,"S","I")
ur_util$PDAST_rpp_aacs <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_access,"S","I")
ur_util$PDAST_rpp_wass <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,watch_singles,"S","I")
ur_util$PDAST_rpp_wacs <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,watch_combos,"S","I")
ur_util$PDAST_rpp_awas <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_watch,"S","I")
ur_util$PDAST_rpp_asr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_singles,"R","NT")
ur_util$PDAST_rpp_acr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_combos,"R","NT")
ur_util$PDAST_rpp_aar <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_abs,"R","NT")
ur_util$PDAST_rpp_isr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,iv_singles,"R","NT")
ur_util$PDAST_rpp_icr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,iv_combos,"R","NT")
ur_util$PDAST_rpp_air <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_ivs,"R","NT")
ur_util$PDAST_rpp_osr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,oral_singles,"R","NT")
ur_util$PDAST_rpp_ocr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,oral_combos,"R","NT")
ur_util$PDAST_rpp_aor <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_orals,"R","NT")
ur_util$PDAST_rpp_acsr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,access_singles,"R","NT")
ur_util$PDAST_rpp_accr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,access_combos,"R","NT")
ur_util$PDAST_rpp_aacr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_access,"R","NT")
ur_util$PDAST_rpp_wasr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,watch_singles,"R","NT")
ur_util$PDAST_rpp_wacr <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,watch_combos,"R","NT")
ur_util$PDAST_rpp_awar <- ur_util %>% number_results_per_panel(PDAST_1,PDAST_6,all_watch,"R","NT")

###Standard results per panel
ur_util$STANDARD_AST_rpp_ass <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_singles,"S","I")
ur_util$STANDARD_AST_rpp_acs <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_combos,"S","I")
ur_util$STANDARD_AST_rpp_aas <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_abs,"S","I")
ur_util$STANDARD_AST_rpp_iss <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_singles,"S","I")
ur_util$STANDARD_AST_rpp_ics <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_combos,"S","I")
ur_util$STANDARD_AST_rpp_ais <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_ivs,"S","I")
ur_util$STANDARD_AST_rpp_oss <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_singles,"S","I")
ur_util$STANDARD_AST_rpp_ocs <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_combos,"S","I")
ur_util$STANDARD_AST_rpp_aos <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_orals,"S","I")
ur_util$STANDARD_AST_rpp_acss <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_singles,"S","I")
ur_util$STANDARD_AST_rpp_accs <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_combos,"S","I")
ur_util$STANDARD_AST_rpp_aacs <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_access,"S","I")
ur_util$STANDARD_AST_rpp_wass <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_singles,"S","I")
ur_util$STANDARD_AST_rpp_wacs <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_combos,"S","I")
ur_util$STANDARD_AST_rpp_awas <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_watch,"S","I")
ur_util$STANDARD_AST_rpp_asr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_singles,"R","NT")
ur_util$STANDARD_AST_rpp_acr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_combos,"R","NT")
ur_util$STANDARD_AST_rpp_aar <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_abs,"R","NT")
ur_util$STANDARD_AST_rpp_isr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_singles,"R","NT")
ur_util$STANDARD_AST_rpp_icr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,iv_combos,"R","NT")
ur_util$STANDARD_AST_rpp_air <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_ivs,"R","NT")
ur_util$STANDARD_AST_rpp_osr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_singles,"R","NT")
ur_util$STANDARD_AST_rpp_ocr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,oral_combos,"R","NT")
ur_util$STANDARD_AST_rpp_aor <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_orals,"R","NT")
ur_util$STANDARD_AST_rpp_acsr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_singles,"R","NT")
ur_util$STANDARD_AST_rpp_accr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,access_combos,"R","NT")
ur_util$STANDARD_AST_rpp_aacr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_access,"R","NT")
ur_util$STANDARD_AST_rpp_wasr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_singles,"R","NT")
ur_util$STANDARD_AST_rpp_wacr <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,watch_combos,"R","NT")
ur_util$STANDARD_AST_rpp_awar <- ur_util %>% number_results_per_panel(STANDARD_AST_1,STANDARD_AST_6,all_watch,"R","NT")

###Assemble data frame for dot plot data visualisation
create_df <- function(df, col_name, aware_result, panel,type,result) {
  result <- df %>%
    select({{ col_name }}) %>%
    cbind(AWaRe_results = aware_result, Panel = panel,Agents=type,Result=result) %>%
    as.data.frame()
  colnames(result) <- c("n", "AWaRe_results", "Panel","Agents","Result")
  return(result)
}

PDAST_rpp_ass <- ur_util %>% create_df(PDAST_rpp_ass,"PDAST\nsingle S","PDAST","Single","S")
PDAST_rpp_acs <- ur_util %>% create_df(PDAST_rpp_acs,"PDAST\ncombo S","PDAST","Combo","S")
PDAST_rpp_aas <- ur_util %>% create_df(PDAST_rpp_aas,"PDAST\nall S","PDAST","All","S")
PDAST_rpp_iss <- ur_util %>% create_df(PDAST_rpp_iss,"PDAST\nsingle IV S","PDAST","Single IV","S")
PDAST_rpp_ics <- ur_util %>% create_df(PDAST_rpp_ics,"PDAST\ncombo IV S","PDAST","Combo IV","S")
PDAST_rpp_ais <- ur_util %>% create_df(PDAST_rpp_ais,"PDAST\nall IV S","PDAST","All IV","S")
PDAST_rpp_oss <- ur_util %>% create_df(PDAST_rpp_oss,"PDAST\nsingle oral S","PDAST","Single Oral","S")
PDAST_rpp_ocs <- ur_util %>% create_df(PDAST_rpp_ocs,"PDAST\ncombo oral S","PDAST","Combo Oral","S")
PDAST_rpp_aos <- ur_util %>% create_df(PDAST_rpp_aos,"PDAST\nall oral S","PDAST","All Oral","S")
PDAST_rpp_acss <- ur_util %>% create_df(PDAST_rpp_acss,"PDAST\nsingle Access S","PDAST","Single Access","S")
PDAST_rpp_accs <- ur_util %>% create_df(PDAST_rpp_accs,"PDAST\ncombo Access S","PDAST","Combo Access","S")
PDAST_rpp_aacs <- ur_util %>% create_df(PDAST_rpp_aacs,"PDAST\nall Access S","PDAST","All Access","S")
PDAST_rpp_wass <- ur_util %>% create_df(PDAST_rpp_wass,"PDAST\nsingle Watch S","PDAST","Single Watch","S")
PDAST_rpp_wacs <- ur_util %>% create_df(PDAST_rpp_wacs,"PDAST\ncombo watch S","PDAST","Combo Watch","S")
PDAST_rpp_awas <- ur_util %>% create_df(PDAST_rpp_awas,"PDAST\nall watch S","PDAST","All Watch","S")
PDAST_rpp_asr <- ur_util %>% create_df(PDAST_rpp_asr,"PDAST\nsingle R","PDAST","Single","R")
PDAST_rpp_acr <- ur_util %>% create_df(PDAST_rpp_acr,"PDAST\ncombo R","PDAST","Combo","R")
PDAST_rpp_aar <- ur_util %>% create_df(PDAST_rpp_aar,"PDAST\nall R","PDAST","All","R")
PDAST_rpp_isr <- ur_util %>% create_df(PDAST_rpp_isr,"PDAST\nsingle IV R","PDAST","Single IV","R")
PDAST_rpp_icr <- ur_util %>% create_df(PDAST_rpp_icr,"PDAST\ncombo IV R","PDAST","Combo IV","R")
PDAST_rpp_air <- ur_util %>% create_df(PDAST_rpp_air,"PDAST\nall IV R","PDAST","All IV","R")
PDAST_rpp_osr <- ur_util %>% create_df(PDAST_rpp_osr,"PDAST\nsingle oral R","PDAST","Single Oral","R")
PDAST_rpp_ocr <- ur_util %>% create_df(PDAST_rpp_ocr,"PDAST\ncombo oral R","PDAST","Combo Oral","R")
PDAST_rpp_aor <- ur_util %>% create_df(PDAST_rpp_aor,"PDAST\nall oral R","PDAST","All Oral","R")
PDAST_rpp_acsr <- ur_util %>% create_df(PDAST_rpp_acsr,"PDAST\nsingle Access R","PDAST","Single Access","R")
PDAST_rpp_accr <- ur_util %>% create_df(PDAST_rpp_accr,"PDAST\ncombo Access R","PDAST","Combo Access","R")
PDAST_rpp_aacr <- ur_util %>% create_df(PDAST_rpp_aacr,"PDAST\nall Access R","PDAST","All Access","R")
PDAST_rpp_wasr <- ur_util %>% create_df(PDAST_rpp_wasr,"PDAST\nsingle Watch R","PDAST","Single Watch","R")
PDAST_rpp_wacr <- ur_util %>% create_df(PDAST_rpp_wacr,"PDAST\ncombo Watch R","PDAST","Combo Watch","R")
PDAST_rpp_awar <- ur_util %>% create_df(PDAST_rpp_awar,"PDAST\nall Watch R","PDAST","All Watch","R")
STANDARD_AST_rpp_ass <- ur_util %>% create_df(STANDARD_AST_rpp_ass,"Standard\nsingle S","Standard","Single","S")
STANDARD_AST_rpp_acs <- ur_util %>% create_df(STANDARD_AST_rpp_acs,"Standard\ncombo S","Standard","Combo","S")
STANDARD_AST_rpp_aas <- ur_util %>% create_df(STANDARD_AST_rpp_aas,"Standard\nall S","Standard","All","S")
STANDARD_AST_rpp_iss <- ur_util %>% create_df(STANDARD_AST_rpp_iss,"Standard\nsingle IV S","Standard","Single IV","S")
STANDARD_AST_rpp_ics <- ur_util %>% create_df(STANDARD_AST_rpp_ics,"Standard\ncombo IV S","Standard","Combo IV","S")
STANDARD_AST_rpp_ais <- ur_util %>% create_df(STANDARD_AST_rpp_ais,"Standard\nall IV S","Standard","All IV","S")
STANDARD_AST_rpp_oss <- ur_util %>% create_df(STANDARD_AST_rpp_oss,"Standard\nsingle oral S","Standard","Single Oral","S")
STANDARD_AST_rpp_ocs <- ur_util %>% create_df(STANDARD_AST_rpp_ocs,"Standard\ncombo oral S","Standard","Combo Oral","S")
STANDARD_AST_rpp_aos <- ur_util %>% create_df(STANDARD_AST_rpp_aos,"Standard\nall oral S","Standard","All Oral","S")
STANDARD_AST_rpp_acss <- ur_util %>% create_df(STANDARD_AST_rpp_acss,"Standard\nsingle Access S","Standard","Single Access","S")
STANDARD_AST_rpp_accs <- ur_util %>% create_df(STANDARD_AST_rpp_accs,"Standard\ncombo Access S","Standard","Combo Access","S")
STANDARD_AST_rpp_aacs <- ur_util %>% create_df(STANDARD_AST_rpp_aacs,"Standard\nall Access S","Standard","All Access","S")
STANDARD_AST_rpp_wass <- ur_util %>% create_df(STANDARD_AST_rpp_wass,"Standard\nsingle Watch S","Standard","Single Watch","S")
STANDARD_AST_rpp_wacs <- ur_util %>% create_df(STANDARD_AST_rpp_wacs,"Standard\ncombo Watch S","Standard","Combo Watch","S")
STANDARD_AST_rpp_awas <- ur_util %>% create_df(STANDARD_AST_rpp_awas,"Standard\nall Watch S","Standard","All Watch","S")
STANDARD_AST_rpp_asr <- ur_util %>% create_df(STANDARD_AST_rpp_asr,"Standard\nsingle R","Standard","Single","R")
STANDARD_AST_rpp_acr <- ur_util %>% create_df(STANDARD_AST_rpp_acr,"Standard\ncombo R","Standard","Combo","R")
STANDARD_AST_rpp_aar <- ur_util %>% create_df(STANDARD_AST_rpp_aar,"Standard\nall R","Standard","All","R")
STANDARD_AST_rpp_isr <- ur_util %>% create_df(STANDARD_AST_rpp_isr,"Standard\nsingle IV R","Standard","Single IV","R")
STANDARD_AST_rpp_icr <- ur_util %>% create_df(STANDARD_AST_rpp_icr,"Standard\ncombo IV R","Standard","Combo IV","R")
STANDARD_AST_rpp_air <- ur_util %>% create_df(STANDARD_AST_rpp_air,"Standard\nall IV R","Standard","All IV","R")
STANDARD_AST_rpp_osr <- ur_util %>% create_df(STANDARD_AST_rpp_osr,"Standard\nsingle oral R","Standard","Single Oral","R")
STANDARD_AST_rpp_ocr <- ur_util %>% create_df(STANDARD_AST_rpp_ocr,"Standard\ncombo oral R","Standard","Combo Oral","R")
STANDARD_AST_rpp_aor <- ur_util %>% create_df(STANDARD_AST_rpp_aor,"Standard\nall oral R","Standard","All Oral","R")
STANDARD_AST_rpp_acsr <- ur_util %>% create_df(STANDARD_AST_rpp_acsr,"Standard\nsingle Access R","Standard","Single Access","R")
STANDARD_AST_rpp_accr <- ur_util %>% create_df(STANDARD_AST_rpp_accr,"Standard\ncombo Access R","Standard","Combo Access","R")
STANDARD_AST_rpp_aacr <- ur_util %>% create_df(STANDARD_AST_rpp_aacr,"Standard\nall Access R","Standard","All Access","R")
STANDARD_AST_rpp_wasr <- ur_util %>% create_df(STANDARD_AST_rpp_wasr,"Standard\nsingle Watch R","Standard","Single Watch","R")
STANDARD_AST_rpp_wacr <- ur_util %>% create_df(STANDARD_AST_rpp_wacr,"Standard\ncombo Watch R","Standard","Combo Watch","R")
STANDARD_AST_rpp_awar <- ur_util %>% create_df(STANDARD_AST_rpp_awar,"Standard\nall Watch R","Standard","All Watch","R")

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

write_csv(summdf,"uf_ast_sourcedata_aware_dotplot.csv")

###Dot plot of number of all S results and Access S results per panel
main_aware_plot <- acs_df %>% main_dotplotter("PDAST\nall S","Standard\nall S","PDAST\nall Access S","Standard\nall Access S",
                                              "All agents","WHO access agents")

###Dot plot of number of all R results and Access R results per panel
accessr_aware_plot <- acs_df %>% main_dotplotter("PDAST\nall R","Standard\nall R","PDAST\nall Access R","Standard\nall Access R",
                                                 "All agents (R)","WHO access agents (R)","(Access agent resistance)")

###Dot plot of number of Watch S results and Watch R results per panel
watch_plot <- acs_df %>% main_dotplotter("PDAST\nall Watch S","Standard\nall Watch S","PDAST\nall Watch R","Standard\nall Watch R",
                                         "WHO watch agents (S)","WHO watch agents (R)","(Watch agent results)")


write_csv(ur_util,"ur_util_microsim.csv")
