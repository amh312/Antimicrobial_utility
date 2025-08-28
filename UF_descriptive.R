#DESCRIPTIVE ANALYSIS

##Functions

  ###Table mutation
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

  ###Antibiogram
  antibio_maker <- function(df,min_no,dataset){
    
    urines5_antibio <- df %>% mutate(
      MEM=case_when(MEM=="NT"~"R",TRUE~MEM),
      SXT=case_when(SXT=="NT"~"R",TRUE~SXT)
    ) %>% mutate(across(AMP:VAN, as.sir)) %>% 
      select(org_fullname,AMP:VAN)
    
    ##NB: those with < 30 results ignored
    u5_antibio <- antibiogram(urines5_antibio,col_mo = "org_fullname",
                              minimum = min_no) %>% 
      arrange(rowSums(across(2:ncol(.))))
    
    abnames <- ab_name(colnames(u5_antibio)[2:ncol(u5_antibio)])
    colnames(u5_antibio) <- c("Organism",abnames)
    
    u5_antibio <- u5_antibio %>% mutate(Organism=str_trim(str_replace(Organism, "\\s*\\(.*", "")))
    
    u5_plotdf <- u5_antibio %>% select(-Organism) %>% mutate(Row=row_number()) %>% 
      pivot_longer(-Row,names_to="Column",values_to="% S")
    
    Orglist <- u5_antibio %>% pull(Organism)
    
    antiplot <- ggplot(u5_plotdf, aes(x = Column, y = factor(Row), fill = `% S`)) +
      geom_tile(color = "white") +
      geom_text(aes(label = `% S`), color = "black") +
      scale_fill_gradient(low = "red", high = "green") +
      theme_minimal() +
      labs(x = "", y = "Row") +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_y_discrete(labels=Orglist)+
      ggtitle(glue("Susceptibility antibiogram of {dataset} dataset"))+
      ylab("Organism")
    
    ggsave(glue("antibio_{dataset}.pdf"), plot = antiplot, device = "pdf", width = 10, height = 8,
           path="/Users/alexhoward/Documents/Projects/UDAST_code")
    
    print(antiplot)
    
  }
  
  ###Organism data for table
  orgcounter <- function(df,dataset,other=TRUE) {
    
    if (other==TRUE) {
    
    df %>% count(org_genus) %>% arrange(desc(n)) %>% 
      mutate(org_genus=case_when(n<10~"Other",TRUE~org_genus)) %>% 
      group_by(org_genus) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      distinct(org_genus, .keep_all = TRUE) %>% ungroup() %>% 
      arrange(desc(n)) %>% mutate(
        Characteristic="Genus grown"
      ) %>% rename(Subtype="org_genus") %>% 
      mutate(perc=round((n/sum(n))*100,1)) %>% 
      mutate(perc=glue("{n} ({perc})")) %>% 
      relocate(Characteristic,.before="Subtype") %>% 
      select(-n)
      
    } else {
      
      df %>% count(org_genus) %>% arrange(desc(n)) %>%
        mutate(
          Characteristic="Genus grown"
        ) %>% rename(Subtype="org_genus") %>% 
        mutate(perc=round((n/sum(n))*100,1)) %>% 
        mutate(perc=glue("{n} ({perc})")) %>% 
        relocate(Characteristic,.before="Subtype") %>% 
        select(-n)
      
    }
    
  }
  
  ###Get diagnoses and pcs for ED stay ID
  rel_edstays <- function(df,ed_df,ed_diag,ed_triag) {
    
    urmod_edstays <- ed_df %>% semi_join(df,by="subject_id")
    urmod_key <- df %>% select(subject_id,charttime)
    urmod_comb <- urmod_edstays %>% left_join(urmod_key,by="subject_id") %>% 
      mutate(ed_urine=case_when(charttime>intime&charttime<outtime~TRUE,
                                TRUE~FALSE)) %>% filter(ed_urine) %>% select(subject_id,stay_id)
    df2 <- df %>% left_join(urmod_comb)
    triage_key <- ed_triag %>% select(stay_id,chiefcomplaint) %>% 
      group_by(stay_id) %>%
      summarise(chiefcomplaint = paste(chiefcomplaint, collapse = ", "), .groups = "drop") %>% 
      ungroup()
    diag_key <- ed_diag %>% select(stay_id,icd_title) %>% 
      group_by(stay_id) %>%
      summarise(icd_title = paste(icd_title, collapse = ", "), .groups = "drop") %>% 
      ungroup()
    df2 %>% left_join(triage_key) %>% left_join(diag_key) %>% 
      mutate(
        chiefcomplaint=case_when(is.na(stay_id)~NA,TRUE~chiefcomplaint),
        icd_title=case_when(is.na(stay_id)~NA,TRUE~icd_title)
      )
      
    
  }
  
##Uploads and reference lists

  ###Uploads
  ur_util <- read_csv("ur_util_final.csv")
  urines_abx <- read_csv("ur_util_final.csv")
  abx <- read_csv("interim_abx.csv")
  urines5_desc <- read_csv("urines5_ref.csv")
  pats <- read_csv("patients.csv")
  hadm <- read_csv("admissions.csv")
  ed_diagnosis <- read_csv("diagnosis.csv")
  edstays <- read_csv("edstays.csv")
  triage <- read_csv("triage.csv")
  
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

##Table 1: Patient/specimen characteristics
  
  ###Year key
  yearkey <- pats %>% distinct(subject_id,.keep_all = T) %>% select(subject_id,anchor_year_group)
  
  ###Age key
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
  
  ###Hospital admission key
  hadmkey <- hadm %>% distinct(subject_id,.keep_all = T) %>% select(subject_id,race,marital_status,language,insurance)
  
  ###Compile descriptive table
  desc_tab <-  abx %>% tab_mutater("Prescription model","N") %>% 
    full_join(urines5_desc %>% tab_mutater("Urine model")) %>% 
    left_join(ur_util %>% tab_mutater("Urine simulation"))
  
  ###Add organism data to table
  model_orgs <- urines5_desc %>% orgcounter() %>% rename(`Urine model n (%)`="perc")
  sim_orgs <- ur_util %>% orgcounter(other=FALSE) %>% rename(`Urine simulation n (%)`="perc")
  all_orgs <- model_orgs %>% left_join(sim_orgs) 
  all_orgs[is.na(all_orgs)] <- "0 (0)"
  all_orgs <- all_orgs %>% mutate(`Prescription model n (%)`=NA) %>% 
    relocate(`Prescription model n (%)`,.after="Urine model n (%)")
  desc_tab <- desc_tab %>% rbind(all_orgs)
 
  ###Illness severity key
  severity <- ur_util %>% count(acuity) %>% mutate(Characteristic="Illness severity score",
                                                   `Prescription model n (%)`=NA,
                                                   `Urine model n (%)`=NA,
                                                   `Urine simulation n (%)`=
                                                     glue("{n} ({round((n/nrow(ur_util))*100,1)})")) %>% 
    rename(Subtype="acuity") %>% relocate(Characteristic,.before="Subtype") %>% 
    select(-n)
  
  ###Bind illness severity to table
  desc_tab <- desc_tab %>% rbind(severity) %>% tibble()
  
  ###Add ED presenting complaints and diagnoses (urine dfs only)
  staykey <- edstays %>% select(stay_id,intime) %>% distinct(stay_id,.keep_all = T) %>% 
    rename(charttime="intime")
  triagekey <- triage %>% left_join(staykey) %>% relocate(charttime,.before = "temperature") %>% 
    mutate(chartdate=as.Date(charttime)) %>% select(subject_id,chartdate,acuity,stay_id) %>% 
    filter(!is.na(acuity))
  compkey <- triage %>% select(stay_id,chiefcomplaint) %>% 
    group_by(stay_id) %>%
    summarise(chiefcomplaint = paste(chiefcomplaint, collapse = ", "), .groups = "drop") %>% 
    ungroup()
  diag_key <- ed_diagnosis %>% select(stay_id,icd_title) %>% 
    group_by(stay_id) %>%
    summarise(icd_title = paste(icd_title, collapse = ", "), .groups = "drop") %>% 
    ungroup()
  
  ed_util <- ur_util %>% select(-c(acuity,chiefcomplaint)) %>% semi_join(triagekey,by=c("subject_id","chartdate")) %>% 
    left_join(triagekey,by=c("subject_id","chartdate"))
  ur_util <- ed_util %>% left_join(compkey) %>% left_join(diag_key)

  ed_urines5 <- urines5_desc %>% semi_join(triagekey,by=c("subject_id","chartdate")) %>% 
    left_join(triagekey,by=c("subject_id","chartdate"))
  ed_urines5 <- ed_urines5 %>% left_join(compkey) %>% left_join(diag_key) %>% 
    select(subject_id,stay_id,icd_title,chiefcomplaint)
  urines5_desc <- urines5_desc %>% left_join(ed_urines5,by="subject_id")
  
  ###Add inpatient diagnoses (all 3 dataframes)
  
  
  
  
  
  
  
  
  
  ###Totals
  totals <- list("Total","Patients",
                 glue("{as.character(nrow(abx %>% distinct(subject_id)))} (100)"),
                 glue("{as.character(nrow(urines5_desc))} (100)"),
                 glue("{as.character(nrow(ur_util))} (100)")
              )
  desc_tab[nrow(desc_tab)+1,] <- totals
  
  ###Remove repetition in characteristic
  desc_tab <- desc_tab %>% mutate(Characteristic=case_when(
    Characteristic==lag(Characteristic)~"",
    TRUE~Characteristic
  ))

  ###Write descriptive table to csv
  write_csv(desc_tab,"uf_desctab.csv")

##Table 2: Prescription characteristics

  ###Clean antimicrobial names
  ab_tab <- abx %>% mutate(ab_name=str_replace_all(ab_name,"_"," & ")) %>% 
    
    ###Count antimicrobials
    count(ab_name) %>% 
    
    ###Arrange descending
    arrange(desc(n)) %>% mutate(ab_name=case_when(
      n/nrow(abx)<0.0025~"Other",TRUE~ab_name
    )) %>% group_by(ab_name) %>%
    
    ###Summarise n by antibiotic
    summarise(n=sum(n)) %>% ungroup() %>% print(n=100) %>% 
    arrange(desc(n)) %>% 
    
    ###Add percentage
    mutate(n=glue("{n} ({round((n/nrow(abx))*100,1)})")) %>%
    rename(`n (%)`="n",`Antibiotic(s)`="ab_name")
  
  ###Write table to csv
  write_csv(ab_tab,"ab_tab.csv")
  
  ##Table 3: Antibiograms
  
  ##Model development data
  antibio_maker(urines5_desc,10,"model development")
  antibio_maker(ur_util,10,"simulation study")
    
