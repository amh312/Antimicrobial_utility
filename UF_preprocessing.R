#FEATURE ENGINEERING

##Functions

###Assigning previous event feature variable
prev_event_assign <- function(df,B_var,event_df,event_var,no_days,no_events) {
  
  #convert charttime to posix
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  
  event_df %>%
    
    #syncing var types and getting rid of nas
    mutate(event = {{event_var}}) %>% 
    select(subject_id, event, charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(!is.na(event)) %>% 
    
    #bind event df of interest
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #check for prev events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #ad prev event var
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter to just urines
    filter(grepl('URINE', spec_type_desc))
  
}

###Dummy variables for main presenting complains
pc_dummies <- function(df) {
  
  df %>% mutate(
    pc_dyspnea = case_when(grepl("Dyspnea",chiefcomplaint,ignore.case=T)|
                             grepl("shortness",chiefcomplaint,ignore.case=T)|
                             grepl("sob",chiefcomplaint,ignore.case=T)|
                             grepl("hypoxia",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_abdopain = case_when(grepl("abd",chiefcomplaint,ignore.case=T)&
                              grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_confusion = case_when((grepl("altered",chiefcomplaint,ignore.case=T)&
                                grepl("mental",chiefcomplaint,ignore.case=T)|
                                grepl("confus",chiefcomplaint,ignore.case=T))~ TRUE,
                             TRUE~FALSE),
    pc_chestpain = case_when(grepl("chest",chiefcomplaint,ignore.case=T)&
                               grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                             TRUE~FALSE),
    pc_weakness = case_when(grepl("weakness",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_dyspnea = case_when(grepl("fever",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_wound = case_when(grepl("wound",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_fall = case_when(grepl("fall",chiefcomplaint,ignore.case=T)~ TRUE,
                        TRUE~FALSE),
    pc_prbleed = case_when(grepl("brbpr",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_vomiting = case_when(grepl("N/V",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_backpain = case_when(grepl("back",chiefcomplaint,ignore.case=T)&
                              grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_lethargy = case_when(grepl("lethargy",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_diarvom = case_when(grepl("N/V/D",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_diarrhea = case_when(grepl("diarrhea",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_headache = case_when(grepl("headache",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_syncope = case_when(grepl("syncope",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_seizure = case_when(grepl("seizure",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_flankpain = case_when(grepl("flank",chiefcomplaint,ignore.case=T)~ TRUE,
                             TRUE~FALSE),
    pc_lowbp = case_when(grepl("hypotension",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_anemia = case_when(grepl("anemia",chiefcomplaint,ignore.case=T)~ TRUE,
                          TRUE~FALSE),
    pc_pain = case_when(grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                        TRUE~FALSE),
    pc_swelling = case_when(grepl("swelling",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_cough = case_when(grepl("cough",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_sepsis = case_when(grepl("cough",chiefcomplaint,ignore.case=T)~ TRUE,
                          TRUE~FALSE),
    pc_fever = case_when(grepl("fever",chiefcomplaint,ignore.case=T)|
                           grepl("pyrexia",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE)
  )
  
}

###Assigning previous event type feature variable
prev_event_type_assign <- function(df,B_var,event_df,event_var,event_type,no_days,no_events) {
  
  df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  event_df %>%
    
    #sync up var types
    mutate(event = {{event_var}}) %>% 
    select(subject_id, event, charttime = 'admittime') %>% 
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(event_type, event)) %>%
    
    #bind event of interest df
    bind_rows(df) %>% 
    mutate(event = case_when(!is.na(event) ~ "Yes",
                             TRUE ~ "No")) %>% 
    
    #check prev events
    MIMER::check_previous_events(cols="event", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make var
    mutate({{B_var}} := case_when(pr_event==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(event = NULL, pr_event=NULL) %>% 
    
    #filter to urines
    filter(grepl('URINE', spec_type_desc))
  
  
}

###Search for previous event across multiple variables
apply_prev_event <- function(df, event,organism) {
  df %>%
    prev_event_type_assign(!!sym(event), urine_df, org_fullname,organism, 365, 1)
}

###Assigning "NT" variable to relevant NAs in microbiology dataframe
NT_assigner <- function(df) {
  
  #split dfs into non_missing and missing
  micaborgs <- df %>% filter(!is.na(org_name))
  micabnas <- df %>% filter(is.na(org_name))
  
  #select ab columns
  micaborgab <- micaborgs %>% select(PEN:MTR)
  
  #replace nas with not tested
  micaborgab[is.na(micaborgab)] <- "NT"
  
  #select corresponding columns
  micaborgs[,17:81] <- micaborgab
  df2 <- tibble(rbind(micaborgs,micabnas))
  df2 %>% rename(admittime = "charttime")
}

###Applying previous AST result search across multiple result types
prev_AST_applier <- function(df1,micro_data,suffix,result,timeframe=365,n_events=1) {
  
  #add p prefix and result type suffix to abx list
  ast_featnames <- paste0("p", antibiotics, suffix)
  
  #redefine prev event func to allow arg passage from wrapping func
  apply_prev_event <- function(df, newvar, antibiotic) {
    df %>%
      prev_event_type_assign(!!sym(newvar), micro_data, !!sym(antibiotic), result, timeframe, n_events)
  }
  
  #use reduce to loop along abx list and iteratively update urines df
  df1 <- reduce(seq_along(antibiotics),
                
                #loop over abx and check for previous specified result
                function(df, i) {
                  apply_prev_event(df, ast_featnames[i], antibiotics[i])
                  }
                
                #define urines df as starting point
                , .init = df1) %>%
    
    #remove leftover grouping from MIMER function
    ungroup()
  
}

###Assigning previous antimicrobial treatment variable
prev_rx_assign <- function(df, B_var, drug_df, abx, abx_groupvar,no_days,no_events) {
  
  #convert charttime to poxix
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  #quosure for var
  abx_groupvar <- enquo(abx_groupvar)
  
  #sync up var types in dfs before bind
  drug_df %>%
    select('subject_id', ab_name,charttime='starttime') %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    filter(grepl(glue("{abx}"), !!abx_groupvar)) %>% 
    
    #bnd df with vars of interest
    bind_rows(ur_df) %>% 
    mutate(abx_treatment = case_when(!is.na(ab_name) ~ "Yes",
                                     TRUE ~ "No")) %>% 
    
    #look for prev vars
    MIMER::check_previous_events(cols="abx_treatment", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_rx_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    
    #make new var
    mutate({{B_var}} := case_when(pr_rx_abx_treatment==TRUE ~ TRUE,
                                  TRUE ~ FALSE)) %>% 
    mutate(abx_treatment=NULL,pr_rx_abx_treatment=NULL) %>% 
    
    #filter back to urines only
    filter(grepl('URINE', spec_type_desc))
  
}

###Finding abnormal inflammatory markers on day prior to urine test
labevent_search <- function(df,search_term,feature_name) {
  
  #quosure for feature var
  feature_name <- enquo(feature_name)
  
  #find result name to filter by
  filter_term <- d_labitems %>%
    filter(grepl(search_term,label,ignore.case=T)) %>% 
    count(itemid) %>% arrange(n) %>% dplyr::slice(1) %>% select(itemid) %>% unlist()
  
  #filter by result name
  filtered_df <- labevents %>% filter(itemid==filter_term) %>% 
    filter(!is.na(valuenum)) %>% rename(admittime="charttime")
  
  #check for prev abnormal result
  df %>% 
    prev_event_type_assign(!!feature_name,filtered_df,flag,"abnormal",1,1) %>%
    ungroup()
  
}

###Assigning gender feature variable
gender_assign <- function(df,B_var,gender_df) {
  
  #make gender key df
  gender_df %>%
    select('subject_id', 'gender') %>%
    
    #join to urine df by key
    right_join(df) %>%
    
    #add new var
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    
    #get rid of residual gender var
    mutate(gender=NULL)
  
}

###Finding patient demographic characeristics
demographic_assign <- function(df,demographic) {
  
  #quosure for demo var
  demographic <- enquo(demographic)
  
  #demographic key df
  hadm_demographic <- hadm %>%
    select(subject_id,!!demographic) %>%
    distinct(subject_id,.keep_all = T)
  
  #join demo key to urine df
  df %>% left_join(hadm_demographic,by="subject_id") %>%
    mutate(!!demographic:=case_when(is.na(!!demographic) ~ "UNKNOWN",
                                    TRUE~!!demographic))
  
}

###Check for hospital admission from outpatient location
outpatient_check <- function(df) {
  
  #if no admission location, label as opd
  df %>% mutate(admission_location=case_when(is.na(admission_location) ~ "OUTPATIENT",
                                             TRUE~admission_location))
  
}

###Applying ICD-1O code search across multiple ICD-10 code prefixes
prev_ICD_applier <- function(df,icd_df,prefix,codes) {
  
  #define prev event applier func
  apply_prev_event_assignments <- function(df, code) {
    
    #add prefix to ICD code
    icd_featname <- paste0(prefix, code)
    
    #check for specified prev icd code using sym to convert str to varname
    df %>%
      prev_event_type_assign(!!sym(icd_featname), icd_df, icd_group, code, 365, 1)
  }
  
  #iteratively update pos_urines with reduce
  pos_urines <- reduce(codes,
                       
                       #loop prev event search across ICD list
                       function(df, code) {
                         apply_prev_event_assignments(df, code)
                         }
                       
                       #specify urines df as starting point for reduce
                       , .init = pos_urines) %>%
    
    #no U codes so remove
    mutate(pDIAG_U = FALSE) %>%
    
    #get rid of MIMER grouping
    ungroup()
  
}

###Checking for previous care events
care_event_assigner <- function(df,search_df,search_term,search_column,feature_name,event_date_col,timeframe,n_events=1) {
  
  #quosure feat name and search col
  feature_name <- enquo(feature_name)
  search_column <- enquo(search_column)
  
  #make key for care event of interest
  care_event <- search_df %>% filter(grepl(search_term,!!search_column,ignore.case=T)) %>% mutate(
    !!search_column:=search_term) %>% rename(admittime=event_date_col)
  
  #check for that event previously
  df %>% 
    prev_event_type_assign(!!feature_name,care_event,!!search_column,search_term,timeframe,n_events) %>%
    ungroup()
  
}

###Applying BMI category search across multiple categories
assign_bmi_events <- function(df, bmi_df, categories, days, min_events) {
  
  #use reduce to apply event search across bmi cats
  reduce(categories,
         
         #loop prev event search across bmi cats with preceding p
         function(acc, category) {
           
           #make features name for bmi
           bmi_featname <- paste0("p", category)
           
           prev_event_type_assign(acc, !!sym(bmi_featname), bmi_df, BMI_cat, category, days, min_events)
           
           }
         
         #specify urines df as starting point
         , .init = df)
}

###Attaching CDI label
CDI_label <- function(df,filter_term) {
  
  #quosure for filter term
  filter_term <- enquo(filter_term)
  
  #ensure spec types all urine
  df <- df %>% mutate(spec_type_desc="URINE") %>% 
    
    #push dates ahead 3m then look back for CDI (to check subsequent 3m)
    prev_event_type_assign(CDI,cdi_ref,org_fullname,"Clostridioides difficile",
                           28*3,1) %>% 
    
    #get rid of MIMER residual grouping
    ungroup()  %>% 
    
    #use filter term to get back to original df
    filter(!is.na(!!filter_term))
  
  #factorise new CDI var
  df$CDI <- factor(df$CDI)
  
  #return df
  df
  
}

###Attaching previous CDI label
pCDI_label <- function(df,filter_term) {
  
  #quosure filter term
  filter_term <- enquo(filter_term)
  
  #look for c diff in the last 10,000 days
  df %>% 
    prev_event_assign(pCDI,cdi_ref,org_fullname,1e4,1) %>% ungroup() %>% 
    
    #filter back to original df
    filter(!is.na(!!filter_term))
  
}

###Attaching previous hospital admission label
hadm_label <- function(df,filter_term) {
  
  #quosure filter term
  filter_term <- enquo(filter_term)
  
  #check for prev hospital admission in last year
  df %>% 
    prev_event_assign(pHADM,hadm,hadm_id,365,1) %>% 
    
    #get rid of MIMER grouping
    ungroup() %>% 
    
    #filter back to original df
    filter(!is.na(!!filter_term))
  
}

###Attaching AKI label
AKI_label <- function(df,filter_term) {
  
  #quosure filter term
  filter_term <- enquo(filter_term)
  
  #check for aki in last 5d
  df %>% 
    prev_event_type_assign(AKI,creats,AKI3,TRUE,
                           5,1) %>% 
    
    #get rid of mimer grouping
    ungroup() %>% 
    
    #factorise aki variable
    mutate(AKI = factor(AKI)) %>% 
    
    #filter df back to original
    filter(!is.na(!!filter_term))
  
}

###Attaching nephrotoxic drugs
nephrotoxic_join <- function(df) {
  
  #join to key with nephrotoxic drugs
  df %>% left_join(nephrotoxics_key) %>% mutate(
    
    #make neph agent variable
    Nephrotoxic_agent = case_when(is.na(Nephrotoxic_agent) ~ FALSE, TRUE ~ Nephrotoxic_agent)
  )
}

###Attaching contrast
contrast_join <- function(df) {
  
  #join to contrast key df
  df %>% left_join(contrast) %>% mutate(
    
    #make new contrast var
    Contrast = case_when(is.na(Contrast) ~ FALSE, TRUE ~ Contrast)
  )
}

###AKI adjustment for nephrotoxins/contrast
AKI_adjusted_check <- function(df) {
  
  df %>% mutate(AKI = case_when(
    
    #qualify AKI by excluding cases with nephrotoxic agents and contrast
    AKI==TRUE & Nephrotoxic_agent==FALSE & Contrast==FALSE ~ TRUE,
    TRUE ~ FALSE
  ))
  
}

###Attaching previous AKI label
prAKI_label <- function(df,filter_term) {
  
  #quosure for filter term
  filter_term <- enquo(filter_term)
  
  df %>% 
    
    #check for previous aki
    prev_event_assign(prAKI,creats,AKI,1e4,1) %>% 
    
    #get rid of mimer grouping
    ungroup() %>% 
    
    #filter back to original df
    filter(!is.na(!!filter_term))
  
}

###Attaching previous diagnosis label
diag_label <- function(df,filter_term,label,timeframe,newvar) {
  
  #quosure filter term
  filter_term <- enquo(filter_term)
  newvar <- enquo(newvar)
  
  #check for previous diagnosis
  df <- df %>% 
    prev_event_type_assign(!!newvar,hadm,description,
                           label,
                           timeframe,1) %>% 
    
    #get rid of mimer grouping
    ungroup() %>% 
    
    #filter back to original df
    filter(!is.na(!!filter_term))
  
}

###Detection of new low blood test value
new_lowvalue <- function(df,new_colname) {
  
  #quosure new column name
  new_colname <- enquo(new_colname)
  
  #group by patient
  df %>% group_by(subject_id) %>% mutate(
    
    #look for value below bottom end of ref range
    !!new_colname := case_when(
      (valuenum < as.numeric(ref_range_lower)) &
        
        #ensure low value is new by looking at row above
        !(lag(valuenum) < lag(as.numeric(ref_range_lower))) ~ TRUE,
      TRUE~FALSE)
  ) %>% ungroup()
  
}

###Detection of new high blood test value
new_highvalue <- function(df,new_colname) {
  
  #quosure new column name
  new_colname <- enquo(new_colname)
  
  #group by patient
  df %>% group_by(subject_id) %>% mutate(
    
    #look for value above upper end of ref range
    !!new_colname := case_when(
      (valuenum > as.numeric(ref_range_upper)) &
        
        #ensure high value is new by looking at row above
        !(lag(valuenum) > lag(as.numeric(ref_range_upper))) ~ TRUE,
      TRUE~FALSE)
  ) %>% ungroup()
  
}

###Attaching specified abnormal blood test value label
abnormal_label <- function(df,df2,new_column,search_term,filter_term) {
  
  #quosure column names
  filter_term <- enquo(filter_term)
  new_column <- enquo(new_column)
  search_term <- enquo(search_term)
  
  #jump event df backwards a week to adjust to look at that period after
  df2 <- df2 %>% mutate(admittime =
                          charttime - (60*60*24*7))
  
  df %>% 
    
    #look at the 48h-7d period after the urine
    prev_event_type_assign(!!new_column,df2,!!search_term,TRUE,
                           5,1) %>% ungroup() %>% 
    
    #factorise the new column
    mutate(!!new_column := factor(!!new_column)) %>% 
    
    #remove nas
    filter(!is.na(!!filter_term))
  
}

###Attaching bleeding diagnosis
bleed_join <- function(df) {
  
  #join to bleeding diagnosis key
  df %>% left_join(bleeding) %>% mutate(
    
    #make bleeding diagnosis variable
    Bleeding_diagnosis = case_when(is.na(Bleeding_diagnosis) ~ FALSE, TRUE ~ Bleeding_diagnosis)
  )
}

###Attaching cytotoxins
cytotoxic_join <- function(df) {
  
  #join to cytotoxic key df
  df %>% left_join(cytotoxics_key) %>% mutate(
    
    #make cytotoxic var
    Cytotoxic_agent = case_when(is.na(Cytotoxic_agent) ~ FALSE, TRUE ~ Cytotoxic_agent)
  )
}

###Check for marrow suppression adjusted by bleeding/cytotoxins
marrow_check <- function(df) {
  
  #look for cytopenias
  df %>% mutate(marrow_suppress = case_when(
    (leukopenia==TRUE|anaemia==TRUE|thrombocytopenia==TRUE)
    
    #qualify by excluding bleeding and cytotoxics
    & Bleeding_diagnosis==FALSE & Cytotoxic_agent==FALSE~TRUE, TRUE~FALSE
  ))
  
}

###Attaching recent biliary procedure(s)
bil_proc_join <- function(df) {
  
  #join to biliary procedure key df
  df %>% left_join(biliary) %>% mutate(
    
    #make biliary procedure var for same admission
    Biliary_procedure = case_when(is.na(Biliary_procedure) ~ FALSE, TRUE ~ Biliary_procedure)
  )
}

###Adjusted check for hepatotoxicity
lft_check <- function(df) {
  
  #look for deranged lfts
  df %>% mutate(deranged_lfts = case_when(
    (high_alp==TRUE|high_alt==TRUE|high_ast==TRUE) &
      
      #exclude liver failure and biliary procedure from same admission
      pLIVER==FALSE&Biliary_procedure==FALSE~TRUE, TRUE~FALSE
  ))
  
}

###Overall toxicity check
toxicity_check <- function(df) {
  
  #look for renal, liver, haem derangement or coded adverse event
  df %>% mutate(overall_tox = case_when(
    AKI==TRUE|marrow_suppress==TRUE|deranged_lfts==TRUE|abx_ae==TRUE ~TRUE, TRUE~FALSE
  ))
  
}

###Check for sepsis or high obs frequency on admission
update_sepsis <- function(df) {
  
  df %>% mutate(admission_sepsis = case_when(((admission_infection==TRUE |
                                               pc_fever==TRUE) &
                                               (SIRS==TRUE |
                                               ob_freq > median(df$ob_freq)
                                               )) |
                                               admission_sepsis==TRUE |
                                               pc_sepsis==TRUE~TRUE,
                                             TRUE~FALSE))
}

###Check for death in the following 28 days
death_check <- function(df,df2,new_column,filter_term) {
  
  #var quosures
  new_column <- enquo(new_column)
  filter_term <- enquo(filter_term)
  
  #jump event df back 28days to capture time after
  df2 <- df2 %>% mutate(admittime = deathtime-(60*60*24*28))
  
  df %>% 
    
    #look back at the 28d jumped forwards (i.e., looking forwards 28d)
    prev_event_assign(!!new_column,df2,deathtime,
                      28,1) %>% ungroup() %>% 
    
    #factorise new column
    mutate(!!new_column := factor(!!new_column)) %>% 
    filter(!is.na(!!filter_term))
  
}

###Check for prior abnormal BMI categorisations
categorise_bmi <- function(df) {
  df %>%
    
    #filter events df down to bmi checks
    filter(grepl("BMI", result_name)) %>%
    mutate(
      
      #categorise bmis
      BMI_cat = case_when(
        as.numeric(result_value) >= 30 ~ "Obese",
        as.numeric(result_value) >= 25 & as.numeric(result_value) < 30 ~ "Overweight",
        as.numeric(result_value) >= 18.5 & as.numeric(result_value) < 25 ~ "Normal weight",
        as.numeric(result_value) < 18.5 ~ "Underweight"
      ),
      
      #ensure admittime is posix
      admittime = as.POSIXct(chartdate, format = '%Y-%m-%d %H:%M:%S')
    )
}

###Check for ICU in the following 28 days
icu_check <- function(df,df2,new_column,filter_term) {
  
  #var quosures
  new_column <- enquo(new_column)
  filter_term <- enquo(filter_term)
  
  #push event df dates back 28 days to allow looking forwards
  df2 <- df2 %>% mutate(admittime = admittime-(60*60*24*28))
  
  df %>% 
    
    #look back 28 days at jumped back df (effectively looking at 28d after)
    prev_event_type_assign(!!new_column,df2,field_value,"ICU",28,1) %>%
    ungroup() %>% 
    
    #factorise new variable
    mutate(!!new_column := factor(!!new_column)) %>% 
    
    #filter out nas to get back to same df
    filter(!is.na(!!filter_term))
  
}

###Check if still an inpatient 7 days later
inpt7d_check <- function(df){
  
  #make discharge key df
  disch_key <- hadm %>% select(hadm_id,dischtime) %>% distinct(hadm_id,.keep_all = T)
  
  #join to discharge key df
  df %>% left_join(disch_key) %>% 
    
    #make boolean var for inpatient in last 7d
    mutate(inpatient_7d = case_when(dischtime>as_datetime(charttime) 
                                    + (60*60*24*7) ~ TRUE,
                                    TRUE~FALSE))
}

###Check overall sepsis adverse outcomes
sepsis_ae_check <- function(df) {
  
  #look for cases that met sepsis criteria and died, admit icu or in a week
  df %>% mutate(sepsis_ae = case_when(admission_sepsis==TRUE &
                                        (death_28d==TRUE |
                                           icu_28d==TRUE |
                                           inpatient_7d==TRUE) ~ TRUE,
                                      TRUE~FALSE),
                
                #factorise new var
                sepsis_ae=factor(sepsis_ae))
}

###Collapse multiple organisms to single resistance result per sample
res_collapse <- function(df,col_name) { 
  
  #quosure var name
  col_name <- enquo(col_name)
  
  #group by specimen
  df %>% group_by(micro_specimen_id) %>% 
    
    #convert ast results to integers with desired priority
    mutate(!!col_name := case_when(!!col_name=="R"~3,
                                   !!col_name=="I"~2,
                                   !!col_name=="S"~1,
                                   !!col_name=="NT"~0)) %>% 
    
    #revert to the highest value for that specimen
    mutate(!!col_name := max(!!col_name)) %>%
    
    #convert that highest number back to the corresponding result
    mutate(!!col_name := case_when(!!col_name==3~"R",
                                   !!col_name==2~"I",
                                   !!col_name==1~"S",
                                   !!col_name==0~"NT")) %>% 
    ungroup()
  
}

###Wrapping function to apply collapse across multiple resistances
big_res_collapse <- function(df) { 
  
  df %>% 
    res_collapse(AMP) %>% 
    res_collapse(SAM) %>%
    res_collapse(TZP) %>%
    res_collapse(CZO) %>%
    res_collapse(CRO) %>%
    res_collapse(CAZ) %>%
    res_collapse(FEP) %>%
    res_collapse(MEM) %>%
    res_collapse(CIP) %>%
    res_collapse(GEN) %>%
    res_collapse(SXT) %>%
    res_collapse(NIT) %>%
    res_collapse(VAN) %>%
    distinct(micro_specimen_id,.keep_all = T)
  
}

###Assigning splitting index
split_indexer <- function(df,size,seed_no) {
  
  #round down to nearest integer to desired sample proportion
  smp_size <- floor(size * nrow(df))
  
  #set specified seed
  set.seed(seed_no)
  
  #set training index based on sample size
  train_ind <- sample(seq_len(nrow(df)), size = smp_size)
  
}

###Converting multinomial resistance variables to binary variables (full df)
binarise_full_df <- function(df,NT_val,I_val) {
  
  #select abx columns
  urref <- df %>% select(AMP:VAN)
  
  #convert not testeds to desired value
  urref[urref=="NT"] <- NT_val
  
  #convert Is to desired value
  urref[urref=="I"] <- I_val
  
  #replace abx values with new dataframe
  df[,1:13] <- urref
  
  #return df as output
  df
  
}

###Dummy variables for main presenting complains
pc_dummies <- function(df) {
  
  df %>% mutate(
    pc_dyspnea = case_when(grepl("Dyspnea",chiefcomplaint,ignore.case=T)|
                             grepl("shortness",chiefcomplaint,ignore.case=T)|
                             grepl("sob",chiefcomplaint,ignore.case=T)|
                             grepl("hypoxia",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_abdopain = case_when(grepl("abd",chiefcomplaint,ignore.case=T)&
                              grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_confusion = case_when((grepl("altered",chiefcomplaint,ignore.case=T)&
                                grepl("mental",chiefcomplaint,ignore.case=T)|
                                grepl("confus",chiefcomplaint,ignore.case=T))~ TRUE,
                             TRUE~FALSE),
    pc_chestpain = case_when(grepl("chest",chiefcomplaint,ignore.case=T)&
                               grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                             TRUE~FALSE),
    pc_weakness = case_when(grepl("weakness",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_dyspnea = case_when(grepl("fever",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_wound = case_when(grepl("wound",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_fall = case_when(grepl("fall",chiefcomplaint,ignore.case=T)~ TRUE,
                        TRUE~FALSE),
    pc_prbleed = case_when(grepl("brbpr",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_vomiting = case_when(grepl("N/V",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_backpain = case_when(grepl("back",chiefcomplaint,ignore.case=T)&
                              grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_lethargy = case_when(grepl("lethargy",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_diarvom = case_when(grepl("N/V/D",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_diarrhea = case_when(grepl("diarrhea",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_headache = case_when(grepl("headache",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_syncope = case_when(grepl("syncope",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_seizure = case_when(grepl("seizure",chiefcomplaint,ignore.case=T)~ TRUE,
                           TRUE~FALSE),
    pc_flankpain = case_when(grepl("flank",chiefcomplaint,ignore.case=T)~ TRUE,
                             TRUE~FALSE),
    pc_lowbp = case_when(grepl("hypotension",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_anemia = case_when(grepl("anemia",chiefcomplaint,ignore.case=T)~ TRUE,
                          TRUE~FALSE),
    pc_pain = case_when(grepl("pain",chiefcomplaint,ignore.case=T)~ TRUE,
                        TRUE~FALSE),
    pc_swelling = case_when(grepl("swelling",chiefcomplaint,ignore.case=T)~ TRUE,
                            TRUE~FALSE),
    pc_cough = case_when(grepl("cough",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE),
    pc_sepsis = case_when(grepl("cough",chiefcomplaint,ignore.case=T)~ TRUE,
                          TRUE~FALSE),
    pc_fever = case_when(grepl("fever",chiefcomplaint,ignore.case=T)|
                           grepl("pyrexia",chiefcomplaint,ignore.case=T)~ TRUE,
                         TRUE~FALSE)
  )
  
}

###Check for combinations
combocheck <- function(df) {
  
  df %>%
    
    #arranging and grouping
    arrange(subject_id, starttime) %>%
    group_by(subject_id) %>%
    
    #look back 36 rows for drug stoptimes that overlap this drug's starttime
    mutate(combination_drug2 = case_when(
      starttime < lag(stoptime) &
        ab_name!=lag(ab_name)~ lag(ab_name),
      TRUE~NA
    ),
    combination_drug3 = case_when(
      starttime < lag(stoptime,n=2) &
        ab_name!=lag(ab_name,n=2)~ lag(ab_name,n=2),
      TRUE~NA
    ),
    combination_drug4 = case_when(
      starttime < lag(stoptime,n=3) &
        ab_name!=lag(ab_name,n=3)~ lag(ab_name,n=3),
      TRUE~NA
    ),
    combination_drug5 = case_when(
      starttime < lag(stoptime,n=4) &
        ab_name!=lag(ab_name,n=4)~ lag(ab_name,n=4),
      TRUE~NA
    ),
    combination_drug6 = case_when(
      starttime < lag(stoptime,n=5) &
        ab_name!=lag(ab_name,n=5)~ lag(ab_name,n=5),
      TRUE~NA
    ),
    combination_drug7 = case_when(
      starttime < lag(stoptime,n=6) &
        ab_name!=lag(ab_name,n=6)~ lag(ab_name,n=6),
      TRUE~NA
    ),
    combination_drug8 = case_when(
      starttime < lag(stoptime,n=7) &
        ab_name!=lag(ab_name,n=7)~ lag(ab_name,n=7),
      TRUE~NA
    ),
    combination_drug9 = case_when(
      starttime < lag(stoptime,n=8) &
        ab_name!=lag(ab_name,n=8)~ lag(ab_name,n=8),
      TRUE~NA
    ),
    combination_drug10 = case_when(
      starttime < lag(stoptime,n=9) &
        ab_name!=lag(ab_name,n=9)~ lag(ab_name,n=9),
      TRUE~NA
    ),
    combination_drug11 = case_when(
      starttime < lag(stoptime,n=10) &
        ab_name!=lag(ab_name,n=10)~ lag(ab_name,n=10),
      TRUE~NA
    ),
    combination_drug12 = case_when(
      starttime < lag(stoptime,n=11) &
        ab_name!=lag(ab_name,n=11)~ lag(ab_name,n=11),
      TRUE~NA
    ),
    combination_drug13 = case_when(
      starttime < lag(stoptime,n=12) &
        ab_name!=lag(ab_name,n=12)~ lag(ab_name,n=12),
      TRUE~NA
    ),
    combination_drug14 = case_when(
      starttime < lag(stoptime,n=13) &
        ab_name!=lag(ab_name,n=13)~ lag(ab_name,n=13),
      TRUE~NA
    ),
    combination_drug15 = case_when(
      starttime < lag(stoptime,n=14) &
        ab_name!=lag(ab_name,n=14)~ lag(ab_name,n=14),
      TRUE~NA
    ),
    combination_drug16 = case_when(
      starttime < lag(stoptime,n=15) &
        ab_name!=lag(ab_name,n=15)~ lag(ab_name,n=15),
      TRUE~NA
    ),
    combination_drug17 = case_when(
      starttime < lag(stoptime,n=16) &
        ab_name!=lag(ab_name,n=16)~ lag(ab_name,n=16),
      TRUE~NA
    ),
    combination_drug18 = case_when(
      starttime < lag(stoptime,n=17) &
        ab_name!=lag(ab_name,n=17)~ lag(ab_name,n=17),
      TRUE~NA
    ),
    combination_drug19 = case_when(
      starttime < lag(stoptime,n=18) &
        ab_name!=lag(ab_name,n=18)~ lag(ab_name,n=18),
      TRUE~NA
    ),
    combination_drug20 = case_when(
      starttime < lag(stoptime,n=19) &
        ab_name!=lag(ab_name,n=19)~ lag(ab_name,n=19),
      TRUE~NA
    ),
    combination_drug21 = case_when(
      starttime < lag(stoptime,n=20) &
        ab_name!=lag(ab_name,n=20)~ lag(ab_name,n=20),
      TRUE~NA
    ),
    combination_drug22 = case_when(
      starttime < lag(stoptime,n=21) &
        ab_name!=lag(ab_name,n=21)~ lag(ab_name,n=21),
      TRUE~NA
    ),
    combination_drug23 = case_when(
      starttime < lag(stoptime,n=22) &
        ab_name!=lag(ab_name,n=22)~ lag(ab_name,n=22),
      TRUE~NA
    ),
    combination_drug24 = case_when(
      starttime < lag(stoptime,n=23) &
        ab_name!=lag(ab_name,n=23)~ lag(ab_name,n=23),
      TRUE~NA
    ),
    combination_drug25 = case_when(
      starttime < lag(stoptime,n=24) &
        ab_name!=lag(ab_name,n=24)~ lag(ab_name,n=24),
      TRUE~NA
    ),
    combination_drug26 = case_when(
      starttime < lag(stoptime,n=25) &
        ab_name!=lag(ab_name,n=25)~ lag(ab_name,n=25),
      TRUE~NA
    ),
    combination_drug27 = case_when(
      starttime < lag(stoptime,n=26) &
        ab_name!=lag(ab_name,n=26)~ lag(ab_name,n=26),
      TRUE~NA
    ),
    combination_drug28 = case_when(
      starttime < lag(stoptime,n=27) &
        ab_name!=lag(ab_name,n=27)~ lag(ab_name,n=27),
      TRUE~NA
    ),
    combination_drug29 = case_when(
      starttime < lag(stoptime,n=28) &
        ab_name!=lag(ab_name,n=28)~ lag(ab_name,n=28),
      TRUE~NA
    ),
    combination_drug30 = case_when(
      starttime < lag(stoptime,n=29) &
        ab_name!=lag(ab_name,n=29)~ lag(ab_name,n=29),
      TRUE~NA
    ),
    combination_drug31 = case_when(
      starttime < lag(stoptime,n=30) &
        ab_name!=lag(ab_name,n=30)~ lag(ab_name,n=30),
      TRUE~NA
    ),
    combination_drug32 = case_when(
      starttime < lag(stoptime,n=31) &
        ab_name!=lag(ab_name,n=31)~ lag(ab_name,n=31),
      TRUE~NA
    ),
    combination_drug33 = case_when(
      starttime < lag(stoptime,n=32) &
        ab_name!=lag(ab_name,n=32)~ lag(ab_name,n=32),
      TRUE~NA
    ),
    combination_drug34 = case_when(
      starttime < lag(stoptime,n=33) &
        ab_name!=lag(ab_name,n=33)~ lag(ab_name,n=33),
      TRUE~NA
    ),
    combination_drug35 = case_when(
      starttime < lag(stoptime,n=34) &
        ab_name!=lag(ab_name,n=34)~ lag(ab_name,n=34),
      TRUE~NA
    ),
    combination_drug36 = case_when(
      starttime < lag(stoptime,n=35) &
        ab_name!=lag(ab_name,n=35)~ lag(ab_name,n=35),
      TRUE~NA
    ),
    combination_drug37 = case_when(
      starttime < lag(stoptime,n=36) &
        ab_name!=lag(ab_name,n=36)~ lag(ab_name,n=36),
      TRUE~NA
    )) %>%
    ungroup()
  
}
combodrugcheck <- function(df) {
  
  #make boolean for the presence of any combinations of abx
  df %>% group_by(subject_id) %>% 
    mutate(
      combo_agent = case_when(
        any(!is.na(combination_drug2))|
          any(!is.na(combination_drug3))|
          any(!is.na(combination_drug4))|
          any(!is.na(combination_drug5))|
          any(!is.na(combination_drug6))|
          any(!is.na(combination_drug7))|
          any(!is.na(combination_drug8))|
          any(!is.na(combination_drug9))|
          any(!is.na(combination_drug10))|
          any(!is.na(combination_drug11))|
          any(!is.na(combination_drug12))|
          any(!is.na(combination_drug13))|
          any(!is.na(combination_drug14))|
          any(!is.na(combination_drug15))|
          any(!is.na(combination_drug16))|
          any(!is.na(combination_drug17))|
          any(!is.na(combination_drug18))|
          any(!is.na(combination_drug19))|
          any(!is.na(combination_drug20))|
          any(!is.na(combination_drug21))|
          any(!is.na(combination_drug22))|
          any(!is.na(combination_drug23))|
          any(!is.na(combination_drug24))|
          any(!is.na(combination_drug25))|
          any(!is.na(combination_drug26))|
          any(!is.na(combination_drug27))|
          any(!is.na(combination_drug28))|
          any(!is.na(combination_drug29))|
          any(!is.na(combination_drug30))|
          any(!is.na(combination_drug31))|
          any(!is.na(combination_drug32))|
          any(!is.na(combination_drug33))|
          any(!is.na(combination_drug34))|
          any(!is.na(combination_drug35))|
          any(!is.na(combination_drug36))|
          any(!is.na(combination_drug37))~ TRUE,
        TRUE~FALSE
      )
    ) %>% ungroup()
  
}

###Prior agent start and stop times for combinations
prior_agent_stoptimes <- function(df) {
  
  df %>%
    
    #arranging and grouping
    arrange(subject_id, starttime) %>%
    group_by(subject_id) %>%
    
    #get stop times of any agents given combined with this drug
    mutate(prioragent_stoptime2 = case_when(
      !is.na(combination_drug2)~lag(stoptime),
      TRUE~NA
    ),
    prioragent_stoptime3 = case_when(
      !is.na(combination_drug3)~lag(stoptime,n=2),
      TRUE~NA
    ),
    prioragent_stoptime4 = case_when(
      !is.na(combination_drug4)~lag(stoptime,n=3),
      TRUE~NA
    ),
    prioragent_stoptime5 = case_when(
      !is.na(combination_drug5)~lag(stoptime,n=4),
      TRUE~NA
    ),
    prioragent_stoptime6 = case_when(
      !is.na(combination_drug6)~lag(stoptime,n=5),
      TRUE~NA
    ),
    prioragent_stoptime7 = case_when(
      !is.na(combination_drug7)~lag(stoptime,n=6),
      TRUE~NA
    ),
    prioragent_stoptime8 = case_when(
      !is.na(combination_drug8)~lag(stoptime,n=7),
      TRUE~NA
    ),
    prioragent_stoptime9 = case_when(
      !is.na(combination_drug9)~lag(stoptime,n=8),
      TRUE~NA
    ),
    prioragent_stoptime10 = case_when(
      !is.na(combination_drug10)~lag(stoptime,n=9),
      TRUE~NA
    ),
    prioragent_stoptime11 = case_when(
      !is.na(combination_drug11)~lag(stoptime,n=10),
      TRUE~NA
    ),
    prioragent_stoptime12 = case_when(
      !is.na(combination_drug12)~lag(stoptime,n=11),
      TRUE~NA
    ),
    prioragent_stoptime13 = case_when(
      !is.na(combination_drug13)~lag(stoptime,n=12),
      TRUE~NA
    ),
    prioragent_stoptime14 = case_when(
      !is.na(combination_drug14)~lag(stoptime,n=13),
      TRUE~NA
    ),
    prioragent_stoptime15 = case_when(
      !is.na(combination_drug15)~lag(stoptime,n=14),
      TRUE~NA
    ),
    prioragent_stoptime16 = case_when(
      !is.na(combination_drug16)~lag(stoptime,n=15),
      TRUE~NA
    ),
    prioragent_stoptime17 = case_when(
      !is.na(combination_drug17)~lag(stoptime,n=16),
      TRUE~NA
    ),
    prioragent_stoptime18 = case_when(
      !is.na(combination_drug18)~lag(stoptime,n=17),
      TRUE~NA
    ),
    prioragent_stoptime19 = case_when(
      !is.na(combination_drug19)~lag(stoptime,n=18),
      TRUE~NA
    ),
    prioragent_stoptime20 = case_when(
      !is.na(combination_drug20)~lag(stoptime,n=19),
      TRUE~NA
    ),
    prioragent_stoptime21 = case_when(
      !is.na(combination_drug21)~lag(stoptime,n=20),
      TRUE~NA
    ),
    prioragent_stoptime22 = case_when(
      !is.na(combination_drug22)~lag(stoptime,n=21),
      TRUE~NA
    ),
    prioragent_stoptime23 = case_when(
      !is.na(combination_drug23)~lag(stoptime,n=22),
      TRUE~NA
    ),
    prioragent_stoptime24 = case_when(
      !is.na(combination_drug24)~lag(stoptime,n=23),
      TRUE~NA
    ),
    prioragent_stoptime25 = case_when(
      !is.na(combination_drug25)~lag(stoptime,n=24),
      TRUE~NA
    ),
    prioragent_stoptime26 = case_when(
      !is.na(combination_drug26)~lag(stoptime,n=25),
      TRUE~NA
    ),
    prioragent_stoptime27 = case_when(
      !is.na(combination_drug27)~lag(stoptime,n=26),
      TRUE~NA
    ),
    prioragent_stoptime28 = case_when(
      !is.na(combination_drug28)~lag(stoptime,n=27),
      TRUE~NA
    ),
    prioragent_stoptime29 = case_when(
      !is.na(combination_drug29)~lag(stoptime,n=28),
      TRUE~NA
    ),
    prioragent_stoptime30 = case_when(
      !is.na(combination_drug30)~lag(stoptime,n=29),
      TRUE~NA
    ),
    prioragent_stoptime31 = case_when(
      !is.na(combination_drug31)~lag(stoptime,n=30),
      TRUE~NA
    ),
    prioragent_stoptime32 = case_when(
      !is.na(combination_drug32)~lag(stoptime,n=31),
      TRUE~NA
    ),
    prioragent_stoptime33 = case_when(
      !is.na(combination_drug33)~lag(stoptime,n=32),
      TRUE~NA
    ),
    prioragent_stoptime34 = case_when(
      !is.na(combination_drug34)~lag(stoptime,n=33),
      TRUE~NA
    ),
    prioragent_stoptime35 = case_when(
      !is.na(combination_drug35)~lag(stoptime,n=34),
      TRUE~NA
    ),
    prioragent_stoptime36 = case_when(
      !is.na(combination_drug36)~lag(stoptime,n=35),
      TRUE~NA
    ),
    prioragent_stoptime37 = case_when(
      !is.na(combination_drug37)~lag(stoptime,n=36),
      TRUE~NA
    )) %>%
    ungroup()
  
}
starttime_check <- function(df) {
  
  
  df %>%
    
    #use map to apply function rowwise
    mutate(combo_starttime = map_dbl(
      
      #map function over list of columns in parallel
      pmap(list(
      starttime,prioragent_stoptime2,prioragent_stoptime3,
      prioragent_stoptime4,prioragent_stoptime5,prioragent_stoptime6,
      prioragent_stoptime7,prioragent_stoptime8,prioragent_stoptime9,
      prioragent_stoptime10,prioragent_stoptime11,prioragent_stoptime12,
      prioragent_stoptime13,prioragent_stoptime14,prioragent_stoptime15,
      prioragent_stoptime16,prioragent_stoptime17,prioragent_stoptime18,
      prioragent_stoptime19,prioragent_stoptime20,prioragent_stoptime21,
      prioragent_stoptime22,prioragent_stoptime23,prioragent_stoptime24,
      prioragent_stoptime25,prioragent_stoptime26,prioragent_stoptime27,
      prioragent_stoptime28,prioragent_stoptime29,prioragent_stoptime30,
      prioragent_stoptime31,prioragent_stoptime32,prioragent_stoptime33,
      prioragent_stoptime34,prioragent_stoptime35,prioragent_stoptime36,
      prioragent_stoptime37
      
      #use min to retrieve earliest starttime from the above columns
    ), ~ min(c(...), na.rm = TRUE)),
    
    #ensure result returned as posix date
    as.POSIXct))
}
curr_regime_stoptimes <- function(df) {
  
  df %>%
    
    #arranging and grouping
    arrange(subject_id, starttime) %>%
    group_by(subject_id) %>%
    
    #check ahead for when current combination ends
    mutate(curr_regime_stoptime2 = case_when(
      !is.na(lead(combination_drug2))~lead(starttime),
      TRUE~NA
    ),
    curr_regime_stoptime3 = case_when(
      !is.na(lead(combination_drug3))~lead(starttime,n=2),
      TRUE~NA
    ),
    curr_regime_stoptime4 = case_when(
      !is.na(lead(combination_drug4))~lead(starttime,n=3),
      TRUE~NA
    ),
    curr_regime_stoptime5 = case_when(
      !is.na(lead(combination_drug5))~lead(starttime,n=4),
      TRUE~NA
    ),
    curr_regime_stoptime6 = case_when(
      !is.na(lead(combination_drug6))~lead(starttime,n=5),
      TRUE~NA
    ),
    curr_regime_stoptime7 = case_when(
      !is.na(lead(combination_drug7))~lead(starttime,n=6),
      TRUE~NA
    ),
    curr_regime_stoptime8 = case_when(
      !is.na(lead(combination_drug8))~lead(starttime,n=7),
      TRUE~NA
    ),
    curr_regime_stoptime9 = case_when(
      !is.na(lead(combination_drug9))~lead(starttime,n=8),
      TRUE~NA
    ),
    curr_regime_stoptime10 = case_when(
      !is.na(lead(combination_drug10))~lead(starttime,n=9),
      TRUE~NA
    ),
    curr_regime_stoptime11 = case_when(
      !is.na(lead(combination_drug11))~lead(starttime,n=10),
      TRUE~NA
    ),
    curr_regime_stoptime12 = case_when(
      !is.na(lead(combination_drug12))~lead(starttime,n=11),
      TRUE~NA
    ),
    curr_regime_stoptime13 = case_when(
      !is.na(lead(combination_drug13))~lead(starttime,n=12),
      TRUE~NA
    ),
    curr_regime_stoptime14 = case_when(
      !is.na(lead(combination_drug14))~lead(starttime,n=13),
      TRUE~NA
    ),
    curr_regime_stoptime15 = case_when(
      !is.na(lead(combination_drug15))~lead(starttime,n=14),
      TRUE~NA
    ),
    curr_regime_stoptime16 = case_when(
      !is.na(lead(combination_drug16))~lead(starttime,n=15),
      TRUE~NA
    ),
    curr_regime_stoptime17 = case_when(
      !is.na(lead(combination_drug17))~lead(starttime,n=16),
      TRUE~NA
    ),
    curr_regime_stoptime18 = case_when(
      !is.na(lead(combination_drug18))~lead(starttime,n=17),
      TRUE~NA
    ),
    curr_regime_stoptime19 = case_when(
      !is.na(lead(combination_drug19))~lead(starttime,n=18),
      TRUE~NA
    ),
    curr_regime_stoptime20 = case_when(
      !is.na(lead(combination_drug20))~lead(starttime,n=19),
      TRUE~NA
    ),
    curr_regime_stoptime21 = case_when(
      !is.na(lead(combination_drug21))~lead(starttime,n=20),
      TRUE~NA
    ),
    curr_regime_stoptime22 = case_when(
      !is.na(lead(combination_drug22))~lead(starttime,n=21),
      TRUE~NA
    ),
    curr_regime_stoptime23 = case_when(
      !is.na(lead(combination_drug23))~lead(starttime,n=22),
      TRUE~NA
    ),
    curr_regime_stoptime24 = case_when(
      !is.na(lead(combination_drug24))~lead(starttime,n=23),
      TRUE~NA
    ),
    curr_regime_stoptime25 = case_when(
      !is.na(lead(combination_drug25))~lead(starttime,n=24),
      TRUE~NA
    ),
    curr_regime_stoptime26 = case_when(
      !is.na(lead(combination_drug26))~lead(starttime,n=25),
      TRUE~NA
    ),
    curr_regime_stoptime27 = case_when(
      !is.na(lead(combination_drug27))~lead(starttime,n=26),
      TRUE~NA
    ),
    curr_regime_stoptime28 = case_when(
      !is.na(lead(combination_drug28))~lead(starttime,n=27),
      TRUE~NA
    ),
    curr_regime_stoptime29 = case_when(
      !is.na(lead(combination_drug29))~lead(starttime,n=28),
      TRUE~NA
    ),
    curr_regime_stoptime30 = case_when(
      !is.na(lead(combination_drug30))~lead(starttime,n=29),
      TRUE~NA
    ),
    curr_regime_stoptime31 = case_when(
      !is.na(lead(combination_drug31))~lead(starttime,n=30),
      TRUE~NA
    ),
    curr_regime_stoptime32 = case_when(
      !is.na(lead(combination_drug32))~lead(starttime,n=31),
      TRUE~NA
    ),
    curr_regime_stoptime33 = case_when(
      !is.na(lead(combination_drug33))~lead(starttime,n=32),
      TRUE~NA
    ),
    curr_regime_stoptime34 = case_when(
      !is.na(lead(combination_drug34))~lead(starttime,n=33),
      TRUE~NA
    ),
    curr_regime_stoptime35 = case_when(
      !is.na(lead(combination_drug35))~lead(starttime,n=34),
      TRUE~NA
    ),
    curr_regime_stoptime36 = case_when(
      !is.na(lead(combination_drug36))~lead(starttime,n=35),
      TRUE~NA
    ),
    curr_regime_stoptime37 = case_when(
      !is.na(lead(combination_drug37))~lead(starttime,n=36),
      TRUE~NA
    )) %>%
    ungroup()
  
}
stoptime_check <- function(df) {
  df %>%
    
    #map_dbl to allow application of posix reformat over output vector
    mutate(combo_stoptime = map_dbl(
      
      #pmap to allow application of min value check across list of vectors
      pmap(list(
      stoptime,curr_regime_stoptime2,curr_regime_stoptime3,
      curr_regime_stoptime4,curr_regime_stoptime5,curr_regime_stoptime6,
      curr_regime_stoptime7,curr_regime_stoptime8,curr_regime_stoptime9,
      curr_regime_stoptime10,curr_regime_stoptime11,curr_regime_stoptime12,
      curr_regime_stoptime13,curr_regime_stoptime14,curr_regime_stoptime15,
      curr_regime_stoptime16,curr_regime_stoptime17,curr_regime_stoptime18,
      curr_regime_stoptime19,curr_regime_stoptime20,curr_regime_stoptime21,
      curr_regime_stoptime22,curr_regime_stoptime23,curr_regime_stoptime24,
      curr_regime_stoptime25,curr_regime_stoptime26,curr_regime_stoptime27,
      curr_regime_stoptime28,curr_regime_stoptime29,curr_regime_stoptime30,
      curr_regime_stoptime31,curr_regime_stoptime32,curr_regime_stoptime33,
      curr_regime_stoptime34,curr_regime_stoptime35,curr_regime_stoptime36,
      curr_regime_stoptime37
      
      #apply min across list of stoptime vectors to produce vector of earliest stoptimes
    ), ~ min(c(...), na.rm = TRUE)),
    
    #apply posix across the output vector to convert back to datetime
    as.POSIXct))
}

###Check if second agent of combination
combo_booleans <- function(df) {
  
  df <- df %>%
    
    #combine drugs in combination into single column with underscore seperator
    unite("all_combo_abs", ab_name,combination_drug2:combination_drug37, sep = "_", remove = FALSE,na.rm=TRUE)
  
  #add suffix to ampicillin to differentiate from ampicillin-sulbactam
  df <- df %>% mutate(
    all_combo_abs = paste0(all_combo_abs,"_end")
  )
  
  #booleans for presence of agents in new united variable
  df %>% mutate(
    Ampicillin = case_when(
      grepl("Ampicillin_",all_combo_abs)|grepl("Ampicillin_end",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    `Ampicillin/sulbactam` = case_when(
      grepl("Ampicillin/sulbactam",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    `Piperacillin/tazobactam` = case_when(
      grepl("Piperacillin/tazobactam",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Cefazolin = case_when(
      grepl("Cefazolin",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Ceftriaxone = case_when(
      grepl("Ceftriaxone",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Ceftazidime = case_when(
      grepl("Ceftazidime",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Cefepime = case_when(
      grepl("Cefepime",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Meropenem = case_when(
      grepl("Meropenem",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Ciprofloxacin = case_when(
      grepl("Ciprofloxacin",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Gentamicin = case_when(
      grepl("Gentamicin",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    `Trimethoprim/sulfamethoxazole` = case_when(
      grepl("Trimethoprim/sulfamethoxazole",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Nitrofurantoin = case_when(
      grepl("Nitrofurantoin",all_combo_abs) ~ TRUE,
      TRUE~FALSE),
    Vancomycin = case_when(
      grepl("Vancomycin",all_combo_abs) ~ TRUE,
      TRUE~FALSE)
  )
  
}

###Count number of drugs in combination
count_drugs <- function(df) {
  
  #count number of drugs in combination variable
  df %>% mutate(drug_count=Ampicillin+`Ampicillin/sulbactam`+
                  `Piperacillin/tazobactam`+Cefazolin+Ceftriaxone+
                  Ceftazidime+Cefepime+
                  Meropenem+Ciprofloxacin+
                  Gentamicin+`Trimethoprim/sulfamethoxazole`+
                  Nitrofurantoin+Vancomycin)
}

###Replace name, starttimes, and stoptimes with combination values
combo_mutate <- function(df) {
  df %>%
    mutate(ab_combination = 
             
             #select specified range of columns
             apply(select(., start_col:end_col), 1, 
                                  function(row) {
                                    
                                    #get names in that row
                                    true_names <- names(row)[row]
                                    
                                    #concatenate abx column names where true
                                    str_c(true_names, collapse = "_")
                                    }
    ))
}
times_amend <- function(df) {
  
  #ensure there is either a single or combination starttime
  df %>% mutate(ab_name=ab_combination,ab_name=ab_combination,
                starttime=case_when(
                  !is.na(combo_starttime)~as_datetime(combo_starttime),
                  TRUE~as_datetime(starttime)
                ),
                
                #do the same for stoptime
                stoptime=case_when(
                  !is.na(combo_stoptime)~as_datetime(combo_stoptime),
                  TRUE~as_datetime(stoptime)
                )) %>% 
    select(-(startdate:ab_combination)) %>% select(-all_combo_abs)
}

###Filtering and binding combinations to full dataset
dup_remove <- function(df) {
  
  #distinct by patient, antibiotic and starttime
  df %>% distinct(subject_id,ab_name,as.Date(starttime),.keep_all = T) %>% 
    
    #remove starttimes
    select(-`as.Date(starttime)`)
}
bind_combos <- function(abx_df,combo_df) {
  
  #filter abx df to just single abx
  abx_df <- abx_df %>% anti_join(combo_df,by="poe_id")
  
  #bind single and combination abs togather
  tibble(rbind(abx_df,combo_df)) %>% arrange(subject_id,starttime)
}
rowids <- function(df) {
  
  #get row ids
  df %>% mutate(row_id = seq(1,nrow(df))) %>% 
    
    #move row ids to first column
    relocate(row_id,.before = "subject_id")
}

###Make 2-agent combination susceptibilities
double_sens_columns <- function(df) {
  
  ##
  for (pair in ab_pairs) {
    
    #make vectors with two AST result columns of interest
    col1 <- df[[pair[1]]]
    col2 <- df[[pair[2]]]
    
    #paste abx together to make new columb name
    new_col_name <- paste(pair, collapse = "_")
    
    #if any r results, return r, otherwise s
    df[[new_col_name]] <- ifelse(col1 == 'R' & col2 == 'R', 'R', 'S')
  }
  
  df
  
}

###Amending results for chromosomal AmpCs
AmpC_variable <- function(df) {
  
  #check if organism is contitutive ampc producer
  df %>% mutate(AmpC=case_when(org_fullname %in% ampc_species~TRUE,TRUE~FALSE))
}
AmpC_converter <- function(df) {
  
  #populate resistances conferred by ampc
  df %>% mutate(
    AMP=case_when(AmpC~"R",TRUE~AMP),SAM=case_when(AmpC~"R",TRUE~SAM),
    TZP=case_when(AmpC~"R",TRUE~TZP),CZO=case_when(AmpC~"R",TRUE~CZO),
    CRO=case_when(AmpC~"R",TRUE~CRO),CAZ=case_when(AmpC~"R",TRUE~CAZ),
    AMP=case_when(AmpC~"R",TRUE~AMP),AMP_SAM=case_when(AmpC~"R",TRUE~AMP_SAM),
    AMP_TZP=case_when(AmpC~"R",TRUE~AMP_TZP),AMP_CZO=case_when(AmpC~"R",TRUE~AMP_CZO),
    AMP_CRO=case_when(AmpC~"R",TRUE~AMP_CRO),AMP_CAZ=case_when(AmpC~"R",TRUE~AMP_CAZ),
    SAM_TZP=case_when(AmpC~"R",TRUE~SAM_TZP),SAM_CZO=case_when(AmpC~"R",TRUE~SAM_CZO),
    SAM_CRO=case_when(AmpC~"R",TRUE~SAM_CRO),SAM_CAZ=case_when(AmpC~"R",TRUE~SAM_CAZ),
    TZP_CZO=case_when(AmpC~"R",TRUE~TZP_CZO),TZP_CRO=case_when(AmpC~"R",TRUE~TZP_CRO),
    CZO_CAZ=case_when(AmpC~"R",TRUE~CZO_CAZ),CZO_CRO=case_when(AmpC~"R",TRUE~CZO_CRO),
    CZO_CAZ=case_when(AmpC~"R",TRUE~CZO_CAZ),CRO_CAZ=case_when(AmpC~"R",TRUE~CRO_CAZ)
  )
}

###Key generation for triage data and observations
keymaker <- function(df) {
  
  #make key dataframe with observations
  df %>% select(hadm_id,temperature,heartrate,
                resprate,o2sat,sbp,dbp,
                chiefcomplaint) %>% 
    distinct(hadm_id,.keep_all = T)
  
}

###Binding blood test events to urines df
blood_binder <- function(df,df2,test_name) {
  
  #quosure
  test_name <- enquo(test_name)
  
  #change charttime to character to match target df
  df2 <- df2 %>% mutate(charttime=as.character(charttime))
  
  #group by patient and arrange by charttime
  df %>% 
    bind_rows(df2) %>%
    group_by(subject_id) %>% 
    arrange(charttime) %>% 
    
    #check for blood test in previous 5 rows
    mutate(!!test_name := case_when(is.na(lag(micro_specimen_id))~lag(valuenum),
                                    !is.na(lag(micro_specimen_id))&is.na(lag(micro_specimen_id,n=2))~lag(valuenum,n=2),
                                    !is.na(lag(micro_specimen_id))&!is.na(lag(micro_specimen_id,n=2))&is.na(lag(micro_specimen_id,n=3))~lag(valuenum,n=3),
                                    !is.na(lag(micro_specimen_id))&!is.na(lag(micro_specimen_id,n=2))&!is.na(lag(micro_specimen_id,n=3))&is.na(lag(micro_specimen_id,n=4))~lag(valuenum,n=4),
                                    !is.na(lag(micro_specimen_id))&!is.na(lag(micro_specimen_id,n=2))&!is.na(lag(micro_specimen_id,n=3))&!is.na(lag(micro_specimen_id,n=4))&is.na(lag(micro_specimen_id,n=5))~lag(valuenum,n=5),
                                    TRUE~NA
    )) %>%
    ungroup() %>%
    
    #filter back to target dataframe
    filter(!is.na(micro_specimen_id)) %>% 
    
    #remove linking variables
    select(-(labevent_id:priority))
  
}

###Check for antibiotic prescription in the last year according to vector element
apply_prev_rx365 <- function(df, suffix, antibiotic) {
  
  #make previous antibiotic feature name
  abx_featname <- paste0("p", suffix)
  
  df %>%
    
    #apply previous treatment check to that antibiotic
    prev_rx_assign(!!sym(abx_featname), drugs, antibiotic, ab_name, 365, 1)
}

###Check for antibiotic prescription in the last week according to vector element
apply_prev_rx7 <- function(df, suffix, antibiotic) {
  
  #paste day 7 onto antibiotic name for feature name
  abx_featname <- paste0("d7", suffix)
  
  df %>%
    
    #check for specified previous antibiotic for that antibiotic name
    prev_rx_assign(!!sym(abx_featname), drugs, antibiotic, ab_name, 7, 1)
  
}


##Data upload (CSV files accessible at https://physionet.org/content/mimiciv/2.2/)

omr <- read_csv("omr.csv") #Measurements e.g., height, weight
hadm <- read_csv("admissions.csv") #Admission data
labevents <- read_csv("labevents.csv") #Laboratory tests (non-micro)
d_labitems <- read_csv("d_labitems.csv") #Laboratory test codes
pats <- read_csv("patients.csv") #Patient demographics
services <- read_csv("services.csv") #Service providers
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv") #icd codes
diagnoses_raw <- read_csv("diagnoses_icd.csv") #icd epi
diagnoses <- read_csv("diagnoses_clean.csv") #diagnoses cleaned in UF_clean
procedures <- read_csv("procedures_clean.csv") #procedures cleaned in UF_clean
poe <- read_csv("poe_clean.csv") #care events cleaned in UF_clean
micro <- read_csv("micro_clean2.csv") #micro cleaned in UF_clean
drugs <- read_csv("drugs_clean.csv") #prescriptions cleaned in UF_clean
vitalsign <- read_csv("vitalsign.csv") #vita slign data for observations
triage <- read_csv("triage.csv") #triage info for illness severity score
edstays <- read_csv("edstays.csv") #edstays to link severity scores to admissions
pos_urines <- read_csv("pos_urines_pre_features.csv") #urines cleaned in UF_clean

##Finding previous AST results

###Assigning modified microbiology dataframes to enable prev_event_type_assign
micro3 <- micro %>% rename(admittime = "charttime")
micro2 <- micro %>% NT_assigner()

###At least one resistant isolate in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "GEN", "SXT", "NIT", "VAN", "AMPC", "TCY", "PEN", 
                 "CLI", "LVX", "AMK", "TOB")
pos_urines <- prev_AST_applier(pos_urines,micro3,"r","R")

###At least one previous resistant isolate in the last week
pos_urines <- prev_AST_applier(pos_urines,micro2,"7dr","R",7,1)

###At least one susceptible isolate in the last year
pos_urines <- prev_AST_applier(pos_urines,micro3,"s","S")

###At least one 'I' isolate in the last year
antibiotics <- antibiotics[antibiotics != "AMPC" & 
                             antibiotics != "SXT" & 
                             antibiotics != "VAN" ] #No 'I' results
pos_urines <- prev_AST_applier(pos_urines,micro3,"i","I")

###At least one isolate with the antimicrobial not tested for in the last year
antibiotics <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", 
                 "CIP", "SXT", "VAN", "PEN")
pos_urines <- prev_AST_applier(pos_urines,micro2,"nt","NT")

###At least one growth of top ten common specified organisms in urine in the last year
urine_df <- micro %>% filter(test_name=="URINE CULTURE" & !is.na(org_fullname)) %>% 
  mutate(admittime=charttime)
organisms <- urine_df %>% count(org_fullname) %>% arrange(desc(n)) %>% 
  dplyr::slice(1:10) %>% pull(org_fullname)
org_featnames <- paste0("pG", organisms,"Urine")
pos_urines <- reduce(seq_along(organisms), function(df, i) {
  
  #apply previous event check across organism list
  apply_prev_event(df, org_featnames[i], organisms[i])
  
  #feed in urines df at each iteration
}, .init = pos_urines) %>%
  ungroup()

##Finding previous antimicrobial treatment

###Assigning reference lists of antimicrobials and new feature suffixes
antibiotics <- c("Ampicillin", "Amoxicillin", "Amoxicillin/clavulanic acid", "Ampicillin/sulbactam",
                 "Piperacillin/tazobactam", "Cefazolin",
                 "Ceftriaxone", "Ceftazidime", "Cefepime", "Meropenem", "Ertapenem",
                 "Aztreonam", "Ciprofloxacin", "Levofloxacin", "Gentamicin", "Tobramycin",
                 "Amikacin", "Rifampicin", "Trimethoprim/sulfamethoxazole", "Nitrofurantoin",
                 "Erythromycin", "Clarithromycin", "Azithromycin", "Clindamycin", "Vancomycin",
                 "Metronidazole", "Linezolid", "Daptomycin", "Doxycycline")
suffixes <- antibiotics %>% as.ab() %>% paste0("rx")

###At least one inpatient antimicrobial prescription in the last year
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  
  #apply previous prescription check across antibiotic names
  apply_prev_rx365(df, suffixes[i], antibiotics[i])
  
  #feed urines df back in to update
}, .init = pos_urines) %>%
  ungroup()

###At least one inpatient antimicrobial prescription in the last week
pos_urines <- reduce(seq_along(antibiotics), function(df, i) {
  
  #apply 7d check across full antibiotic list
  apply_prev_rx7(df, suffixes[i], antibiotics[i])
  
  #feed urine df back in to update
}, .init = pos_urines) %>%
  ungroup()

##Find inflammatory marker results

###Elevated C-reactive protein in the 24 hours prior to urine specimen
pos_urines <- pos_urines %>% labevent_search("reactive",highCRP)

###Abnormal total peripheral white cell count in the 24 hours prior to urine specimen
pos_urines <- pos_urines %>% labevent_search("White",abnormalWCC)

##Find patient characteristic and history variables

###At least one hospital admission in the last year
pos_urines <- pos_urines %>% 
  prev_event_assign(pHADM,hadm,hadm_id,365,1) %>%
  ungroup()

###At least one discharge to a nursing home
pos_urines <- pos_urines %>% 
  prev_event_type_assign(pNH,hadm,discharge_location,"NURSING",1e4,1) %>%
  ungroup()

###Male patient
pos_urines <- pos_urines %>% 
  gender_assign(MALE,pats)

###Patient age group at time of admission
agelist <- seq(30,80,10)
pats <- pats %>% mutate(standard_age=anchor_age)
pats <- reduce(seq_along(agelist),function(df,i){
  
  #call anyone under 30 18
  df %>% mutate(standard_age=case_when(anchor_age < 30 ~ 18,
                                       
    #otherwise take the floor value of that decade
    anchor_age>=agelist[i] &
      anchor_age<agelist[i]+10 ~ agelist[i],
    
    #unless over 90, in which case call 90
    anchor_age >=90 ~ 90,TRUE~standard_age))
  
},.init=pats)
pats <- pats %>% group_by(subject_id) %>% summarise(standard_age=mean(standard_age,na.rm=TRUE))
pos_urines <- left_join(pos_urines,pats,by="subject_id")

###Other patient demographics
pos_urines <- pos_urines %>% 
  demographic_assign(race) %>% #Race
  demographic_assign(marital_status) %>% #Marital status
  demographic_assign(insurance) %>% #Insurance type
  demographic_assign(language) #Whether the patient speaks English

###Hospital admission from outpatient location
hadm_admission <- hadm %>%
  select(hadm_id,admission_location) %>% outpatient_check() %>% 
  mutate(hadm_id = case_when(is.na(hadm_id) ~ 0,
                             TRUE ~ hadm_id)) %>% 
  distinct(hadm_id,.keep_all = T)
pos_urines <- left_join(pos_urines,hadm_admission,by="hadm_id") %>% 
  outpatient_check()

###Coded ICD-10 diagnosis groups for admissions in the last year
diag_codes <- LETTERS[LETTERS!="U"]
pos_urines <- pos_urines %>% prev_ICD_applier(diagnoses,"pDIAG_",diag_codes)

###Coded ICD-10 procedure groups for admissions in the last year
proc_codes <- c("0", "3", "8", "5", "T", "4", "S", "A", "9", 
                "H", "I", "B", "7", "G", "1", "R", "J", "Q", 
                "K", "6", "M", "P", "L", "D", "F", "2", "N", 
                "C", "E", "X", "O")
pos_urines <- pos_urines %>% prev_ICD_applier(procedures,"pPROC_",proc_codes)

###At least one coded previous UTI diagnosis in the last year
uti_key <-d_icd_diagnoses %>% filter(grepl("urinary tract infection",long_title,ignore.case=T) |
                                       grepl("acute pyelon",long_title,ignore.case=T) |
                                       (grepl("urinary catheter",long_title,ignore.case=T) & 
                                          grepl("infec",long_title,ignore.case=T)))
hadm_key <- hadm %>% select(hadm_id,admittime)
uti_df <- diagnoses_raw %>% left_join(uti_key,by=c("icd_code","icd_version")) %>% 
  filter(!is.na(long_title)) %>% left_join(hadm_key,by="hadm_id")
pos_urines <- pos_urines %>% care_event_assigner(uti_df,"(urin|pyelo|cath)",long_title,pUTI,"admittime",365,1)

###Presence of an outpatient provider ID
pos_urines <- pos_urines %>% mutate(provider_id = case_when(order_provider_id!="" ~ TRUE,
                                                            TRUE ~ FALSE))

###Current inpatient specialty
serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
pos_urines <- pos_urines %>% left_join(serv_key,by="hadm_id") %>% mutate(
  curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                           TRUE ~ curr_service))

###At least one measured obese, underweight, or overweight BMI category in the last 3 years
bmi <- categorise_bmi(omr)
bmi_categories <- c("Obese", "Underweight", "Overweight")
pos_urines <- assign_bmi_events(pos_urines, bmi, bmi_categories, 1095, 1) %>%
  ungroup()

###Observation frequency on day of test
obs <- poe %>% filter(order_subtype=="Vitals/Monitoring") %>% 
  mutate(ordertime=as.Date(ordertime)) %>% 
  group_by(subject_id,ordertime) %>% count(order_subtype) %>% 
  arrange(desc(n)) %>% select(-order_subtype)
pos_urines <- pos_urines %>% mutate(ordertime=chartdate) %>%
  left_join(obs,by=c("subject_id","ordertime")) %>% 
  rename(ob_freq = "n") %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ 0,
                                                       TRUE ~ ob_freq),
                                   ob_freq = standardize(ob_freq)) %>% 
  select(-ordertime)

###Care events (fields)
fieldlist <- c("cath","DNR","Discharge","ICU","Psychiatry","Nephrostomy",         
                "Surgery","Hydration","NGT","Chemo")
fieldnames <- c("pCATH","pDNR","pDISC","pICU","pPsych","pNeph","pSURG",
                "pHyd","pNGT","pChemo")
fieldtfram <- c(28,365,28,28,365,365,365,28,28,28)
pos_urines <- reduce(seq_along(fieldlist),function(df,i){
  
  df %>%
    
    #check specified event
    care_event_assigner(poe,fieldlist[i],field_value,
                        
                        #populate specified feature name
                        !!sym(fieldnames[i]),"ordertime",
                        
                        #check specified timeframe
                        fieldtfram[i])
    
},.init=pos_urines)

#Care events (orders)
orderlist <- c("Nutrition consult","Physical Therapy","Restraints",
               "Occupational Therapy","Central TPN")
ordernames <- c("pNUTR","pPhysio","pRestr","pOT","pTPN")
pos_urines <- reduce(seq_along(orderlist),function(df,i){
  
  df %>% 
    
    #check for specified event
    care_event_assigner(poe,orderlist[i],order_subtype,
                        
                        #allocate specified feature name, check 365 days
                        !!sym(ordernames[i]),"ordertime",365)
  
},.init=pos_urines)

##Name of organism grown (dummy boolean features)
org_recipe <- recipe(~org_fullname,data=pos_urines)
orgdummyvars <- org_recipe %>% step_dummy(org_fullname) %>% prep(training = pos_urines)
orgvar_cols <- bake(orgdummyvars,new_data = NULL)
pos_urines <- pos_urines %>% cbind(orgvar_cols) %>% tibble()

##Standardised observations

###Observation key from ED datasets
staykey <- edstays %>% select(hadm_id,stay_id) %>% 
  distinct(hadm_id,.keep_all = T) %>% filter(!is.na(hadm_id))
triage <- triage %>% left_join(staykey)
triagekey <- triage %>% select(hadm_id,temperature,heartrate,
                               resprate,o2sat,sbp,dbp,acuity,
                               chiefcomplaint) %>% 
  distinct(hadm_id,.keep_all = T)
pos_urines <- pos_urines %>% left_join(triagekey)
pos_urines <- pos_urines %>% mutate(o2sat = case_when(o2sat>100 ~ 100, #remove implausible sats
                                                      TRUE~o2sat))
pos_urines <- pos_urines %>% mutate(dbp = case_when(dbp==775 ~ 75, #remove implausible hr
                                                      TRUE~dbp))
obslist <- c("heartrate","resprate","sbp","dbp","acuity","o2sat","temperature")
pos_urines <- reduce(seq_along(obslist),function(df,i){
  
  #impute missing obs values with mean of the rest of the column
  df %>% mutate(!!sym(obslist[i]):=case_when(is.na(!!sym(obslist[i]))~
                                               mean(!!sym(obslist[i]),na.rm=T),
                                            TRUE~!!sym(obslist[i]))) %>% 
    
    #standardise obs
    mutate(!!sym(obslist[i]):=standardize(!!sym(obslist[i])))
  
},.init=pos_urines)

###Adding presenting complaint variables
pos_urines %>% count(chiefcomplaint) %>% arrange(desc(n)) %>% print(n=200)
pos_urines <- pos_urines %>% pc_dummies() %>% select(-chiefcomplaint)

###Preceding blood tests
labevents_urine <- labevents %>% semi_join(pos_urines,by="hadm_id")
bltestvec <- c("Alanine Aminotransferase","Albumin","Alkaline Phosphatase",
             "Bicarbonate","Bilirubin, Total","^Creatinine$","Potassium",
             "Sodium","Hematocrit","Hemoglobin","Lymphocytes","Monocytes",
             "Neutrophils","Platelet Count","^RDW$","Red Blood Cells",
             "White Blood Cells","^Lactate$","Lactate Dehydrogenase")
pos_urines <- reduce(seq_along(bltestvec),
                     function(df,i){
                       
                       #filter labitems to test of interest
                       bltest <- d_labitems %>% filter(grepl(bltestvec[i],label))
                       
                       #match labevents to that label
                       labev_bltest <- labevents_urine %>% semi_join(bltest,by="itemid")
                       
                       #remove grepl search punctuation from names to get colnames
                       bltestnam <- gsub("[[:punct:]]+","",bltestvec)
                       
                       df %>% 
                         
                         #bind and check for tests in urine dataframe
                         blood_binder(labev_bltest,!!sym(bltestvec[i])) %>% 
                         
                         #if test is na then impute the mean
                         mutate(
                           !!sym(bltestnam[i]):=case_when(is.na(!!sym(bltestnam[i])) ~ mean(as.numeric(!!sym(bltestnam[i])), na.rm = TRUE),
                                               TRUE ~ as.numeric(!!sym(bltestnam[i])))
                         ) %>% 
                         
                         #standardise values
                         mutate(!!sym(bltestnam[i]):=standardize(!!sym(bltestnam[i])))
                       
                     },.init=pos_urines)

###Arrival transport and disposition (from edstays)
transp_key <- edstays %>% filter(!is.na(hadm_id)&!is.na(stay_id)) %>% 
  select(hadm_id,stay_id,arrival_transport,disposition) %>% distinct(hadm_id,.keep_all = T)
pos_urines <- pos_urines %>% left_join(transp_key,by="hadm_id")
pos_urines <- pos_urines %>% mutate(arrival_transport=case_when(is.na(arrival_transport)~"UNKNOWN",
                                                                TRUE~arrival_transport),
                                    disposition=case_when(is.na(disposition)~"UNKNOWN",
                                                                TRUE~disposition))

###Admission medications
medrecon <- read_csv("medrecon.csv")
medkey <- medrecon %>% count(etcdescription) %>% arrange(desc(n)) %>% filter(n>1000)
medrecon <- medrecon %>% semi_join(medkey,by="etcdescription")
medrecon <- medrecon %>% distinct(stay_id,etcdescription)
medrecon <- medrecon %>% mutate(value = TRUE) %>%
  pivot_wider(
    id_cols = stay_id,                    
    names_from = etcdescription,          
    values_from = value,                  
    values_fill = list(value = FALSE)     
  )
pos_urines <- pos_urines %>% left_join(medrecon,by="stay_id")
pos_urines <- pos_urines %>% select(-stay_id)
pos_urines <- pos_urines %>% mutate(across(`Asthma/COPD Therapy - Beta 2-Adrenergic Agents, Inhaled, Short Acting`:
                               `Hepatitis B Treatment- Nucleoside Analogs (Antiviral)`,
                      ~replace_na(.x,FALSE))) 

write_csv(pos_urines,"pos_urines_w_features.csv")

##Preprocessing for prediction model

###Filter to the last urine for each subject and collapse AST results
pos_urines <- pos_urines %>% group_by(subject_id) %>% 
  arrange(chartdate) %>% summarise_all(last) %>% 
  big_res_collapse()

###Split into stem model development and microsimulation datasets
train_ind <- pos_urines %>% split_indexer(0.9,123)
urines <- pos_urines[train_ind,]
urines_assess <- pos_urines[-train_ind,]

###Modify to include effect of AmpCs
ampc_species <- c("Citrobacter braakii","Citrobacter freundii","Citrobacter gillenii",
                  "Citrobacter murliniae","Citrobacter rodenticum",
                  "Citrobacter sedlakii","Citrobacter werkmanii",
                  "Citrobacter youngae","Enterobacter",
                  "Hafnia alvei","Klebsiella aerogenes",
                  "Morganella morganii","Providencia",
                  "Serratia marcescens")
urines <- urines %>% AmpC_variable()
urines_assess <- urines_assess %>% AmpC_variable()

###Assign and save reference datasets
ur_util <- urines_assess
urines_ref <- urines
urines <- tibble(urines %>% ungroup() %>% select(AMP:VAN,pAMPr:`Hepatitis B Treatment- Nucleoside Analogs (Antiviral)`,AmpC))
write_csv(urines, "urines.csv")
write_csv(urines_ref,"urines_ref.csv")
write_csv(urines_assess,"urines_assess.csv")

###Dataset for enterococcal sensitivity analysis
urines_ref2 <- urines_ref %>% filter(!grepl("Enterococcus",org_fullname))
write_csv(urines_ref2,"urines_ref2.csv")
urines5_ent <- tibble(urines_ref2 %>% ungroup() %>% select(AMP:VAN,pAMPr:`Hepatitis B Treatment- Nucleoside Analogs (Antiviral)`,AmpC))

##Main binomial analysis model development dataset

###Split off dataset from stem
urines5 <- urines

###Select only variables to be used for model development
urines5 <- urines5 %>% select(1:pTPN,temperature:AmpC)
urines5_ent <- urines5_ent %>% select(1:pTPN,temperature:AmpC)

###Binarise AST results
urines5 <- urines5 %>% binarise_full_df("R","S")
urines5_ent <- urines5_ent %>% binarise_full_df("R","S")

###Make 2-agent combination susceptibilities from urines and ur_util dataframes
ab_columns <- c("AMP", "SAM", "TZP", "CZO", "CRO", "CAZ", "FEP", "MEM", "CIP", "GEN", "SXT", "NIT", "VAN")
ab_pairs <- combn(ab_columns, 2, simplify = FALSE)
ur_util <- ur_util %>% double_sens_columns()
urines5 <- urines5 %>% double_sens_columns()
urines5_ent <- urines5_ent %>% double_sens_columns()

###AmpC conversion
ur_util <- ur_util %>% AmpC_converter()
urines5 <- urines5 %>% AmpC_converter() %>% select(-AmpC)
urines5_ent <- urines5_ent %>% AmpC_converter() %>% select(-AmpC)

###Save to file
write_csv(urines5,"urines5.csv")
write_csv(urines5_ent,"urines5_ent.csv")
write_csv(ur_util,"ur_util.csv")

###Drugs data frames for training and testing weighting models
abx1 <- drugs %>% filter(grepl(
  "(Ampicillin|Piperacillin/tazobactam|Cefazolin|Ceftriaxone|^Ceftazidime$|Cefepime|Meropenem|Ciprofloxacin|Gentamicin|Trimethoprim/sulfamethoxazole|Nitrofurantoin)",
  ab_name) | (grepl("Vancomycin",ab_name) & route=="IV")
) %>% filter(!grepl("avibactam",ab_name)) %>%  
  filter(grepl("(PO|NG|IV)",route))
abx1 <- abx1 %>% mutate(charttime = starttime)
abx1 <- abx1 %>% filter(!is.na(starttime) & !is.na(stoptime))
abx_ref <- abx1 %>% semi_join(ur_util,by="subject_id")
abx <- abx1 %>% anti_join(ur_util,by="subject_id")

##CDI-related labels
micro <- read_csv("micro_clean2.csv")
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi <- cdi_ref %>% group_by(subject_id) %>% arrange(chartdate) %>% summarise_all(last)
cdi_ref <- cdi_ref %>% mutate(admittime = charttime-(60*60*24*28*3))

###Attach CDI in following 28d labels to abx dfs
abx <- abx %>% CDI_label(ab_name)
ur_util <- ur_util %>% CDI_label(AMP)

###Add previous hospital admission to abx dataframe
hadm <- read_csv("admissions.csv")
abx <- abx %>% hadm_label(ab_name)

###Check for previous CDI
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi_ref <- cdi_ref %>% mutate(admittime=charttime)
abx <- abx %>% pCDI_label(ab_name)
ur_util <- ur_util %>% pCDI_label(AMP)

##Check for age > 65
pats <- read_csv("patients.csv")
pats <- pats %>% mutate(age65 = case_when(
  anchor_age >=65 ~ TRUE, TRUE~FALSE
))
patskey <- pats %>% select(subject_id,age65)
abx <- abx %>% left_join(patskey,by="subject_id")
ur_util <- ur_util %>% left_join(patskey,by="subject_id")

##Nephrotoxicity-related labels

###Pre-adjusted AKI label
print(d_labitems %>% filter(grepl("creat",label,ignore.case=T)),n=25)
creats <- labevents %>% filter(itemid==50912) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

creats <- creats %>% group_by(subject_id) %>% mutate(
  baseline = min(valuenum),
  match_48h = case_when(
    charttime - lag(charttime) <= 2 ~ TRUE,
    TRUE~ FALSE),
  AKI1 = case_when(
    match_48h &
      valuenum >= 3.54 ~ TRUE,
    TRUE ~ FALSE),
  AKI2 = case_when(
    valuenum >= (3*baseline) ~ TRUE,
    TRUE ~ FALSE),
  AKI3 = case_when(AKI1|AKI2 ~ TRUE, TRUE~FALSE)) %>% 
  ungroup()
write_csv(creats,"creatinines.csv")
creats <- creats %>% mutate(admittime = #adjust to search after rather than before
                              charttime - (60*60*24*7))
abx <- abx %>% AKI_label(ab_name)
ur_util <- ur_util %>% AKI_label(AMP)

###Check for other nephrotoxins (as per BNF list)
nephrotoxins <- c("aceclofenac",	"aciclovir",	"adefovir",
                  "aspirin","captopril",	"carboplatin",	"cefaclor",
                  "cefadroxil",	"cefalexin","cefixime",	
                  "cefoxitin",	"cefradine",	"ceftaroline",
                  "ceftobiprole",	"ceftolozane","cefuroxime",	"celecoxib",
                  "ciclosporin",	"cidofovir",	"cisplatin",
                  "colistimethate",	"deferasirox",	"dexketoprofen",
                  "diclofenac",	"enalapril",	"etodolac",
                  "etoricoxib",	"flurbiprofen",	"foscarnet",
                  "fosinopril",	"ganciclovir","ibuprofen",	"ifosfamide",	"imidapril",
                  "indometacin",	"inotersen",	"ketoprofen",
                  "ketorolac",	"lisinopril",	"lithium",
                  "mefenamic acid",	"meloxicam",	"mesalazine",
                  "methotrexate",	"nabumetone",	"naproxen",
                  "neomycin",	"netilmicin",	"oxaliplatin",
                  "pamidronate",	"parecoxib",	"pemetrexed",
                  "pentamidine",	"perindopril",	"phenazone",
                  "piroxicam",	"quinapril",	"ramipril",
                  "rifampicin",	"streptomycin",	"streptozocin",
                  "sulfasalazine",	"sulindac",	"tacrolimus",
                  "tenofovir disoproxil",	"tenoxicam",	"tiaprofenic acid",
                  "tobramycin",	"tolfenamic acid",	"trandolapril",
                  "trimethoprim",	"valaciclovir",	"valganciclovir",
                  "voclosporin",	"zoledronate")

nephrotoxins <- str_to_title(str_to_lower(nephrotoxins))
nephrotoxics_key <- drugs %>%
  filter(drug %in% nephrotoxins) %>% distinct(hadm_id) %>% 
  mutate(Nephrotoxic_agent = TRUE)
abx <- abx %>% nephrotoxic_join()
ur_util <- ur_util %>% nephrotoxic_join()

###Check for recent IV contrast administration
d_icd_procedures <- read_csv("d_icd_procedures.csv")
procedures <- read_csv("procedures_icd.csv")
contrast_key <- d_icd_procedures %>% filter(grepl("contrast",long_title))
contrast <- procedures %>% semi_join(contrast_key) %>% left_join(hadm_key) %>% 
  distinct(hadm_id) %>% mutate(Contrast = TRUE)
abx <- abx %>% contrast_join()
ur_util <- ur_util %>% contrast_join()

###Adjust AKI label for other nephrotoxins and contrast
abx <- abx %>% AKI_adjusted_check()
ur_util <- ur_util %>% AKI_adjusted_check()

###Check for previous AKI
creats <- creats %>% mutate(admittime=charttime)
creats <- creats %>% mutate(AKI = AKI3)
abx <- abx %>% prAKI_label(ab_name)
ur_util <- ur_util %>% prAKI_label(AMP)

###Check for CKD
drgcodes <- read_csv("drgcodes.csv")
ckdkey <- drgcodes %>% filter(grepl(
  "(CHRONIC KIDNEY|DIAB|LIVER|CARD|HEART|MALIG|STROKE|SEPSIS)",
  description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(ckdkey,by="hadm_id")
abx <- abx %>% diag_label(ab_name,"CHRONIC KIDNEY",1e4,pCKD)
ur_util <- ur_util %>% diag_label(AMP,"CHRONIC KIDNEY",1e4,pCKD)

###Check for diabetes mellitus
abx <- abx %>% diag_label(ab_name,"DIAB",1e4,pDIAB)
ur_util <- ur_util %>% diag_label(AMP,"DIAB",1e4,pDIAB)

###Check for chronic liver disease
abx <- abx %>% diag_label(ab_name,"LIVER",1e4,pLIVER)
ur_util <- ur_util %>% diag_label(AMP,"LIVER",1e4,pLIVER)

###Check for cardiovascular disease
abx <- abx %>% diag_label(ab_name,"(HEART|CARDI)",1e4,pCARD)
ur_util <- ur_util %>% diag_label(AMP,"(HEART|CARDI))",1e4,pCARD)

###Check for malignancy in the last year
abx <- abx %>% diag_label(ab_name,"MALIG",365,pCA)
ur_util <- ur_util %>% diag_label(AMP,"MALIG",365,pCA)

###Check for previous stroke
abx <- abx %>% diag_label(ab_name,"STROKE",1e4,pCVA)
ur_util <- ur_util %>% diag_label(AMP,"STROKE",1e4,pCVA)

###Add male label to abx dataframe
patskey2 <- pats %>% mutate(MALE = case_when(gender=="M" ~ TRUE, TRUE~FALSE)) %>% 
  select(subject_id,MALE)
abx <- abx %>% left_join(patskey2,by="subject_id")

##Sepsis-related labels

###Check for current service provider
services <- read_csv("services.csv")
serv_key <- services %>% select(hadm_id,curr_service) %>% distinct(hadm_id,.keep_all = T)
abx <- abx %>% left_join(serv_key,by="hadm_id") %>% mutate(
  curr_service = case_when(is.na(curr_service) ~ "UNKNOWN",
                           TRUE ~ curr_service))
recipethis <- recipe(~curr_service,data=abx)
dummies <- recipethis %>% step_dummy(curr_service) %>% prep(training = abx)
dummy_data <- bake(dummies,new_data = NULL)
abx <- abx %>% cbind(dummy_data) %>% tibble()
recipethis <- recipe(~curr_service,data=ur_util)
dummies <- recipethis %>% step_dummy(curr_service) %>% prep(training = ur_util)
dummy_data <- bake(dummies,new_data = NULL)
ur_util <- ur_util %>% cbind(dummy_data) %>% tibble()

###Check for recent ICU admission
poe <- read_csv("poe_clean.csv")
icu <- poe %>% filter(field_value=="ICU") %>% mutate(
  field_value="ICU") %>% rename(admittime="ordertime")
abx <- abx %>% 
  prev_event_type_assign(pICU,icu,field_value,"ICU",28,1) %>%
  ungroup()

###Check for recent sepsis
abx <- abx %>% diag_label(ab_name,"SEPSIS",365,pSEPSIS)
ur_util <- ur_util %>% diag_label(AMP,"SEPSIS",365,pSEPSIS)

##Antimicrobial-related labels

###Check length of antimicrobial course
abx <- abx %>% mutate(course_length = as.numeric(stoptime-starttime)/24/60/60,
                      course_length=case_when(course_length<0 ~ 0,TRUE~course_length),
                      course_length=case_when(course_length>=365 ~ course_length-365,
                                              TRUE~course_length),
                      course_length = course_length/max(course_length)) %>% 
  group_by(ab_name) %>% mutate(median_course = median(course_length)) %>% ungroup()

med_course_key <- abx %>% select(ab_name,median_course) %>% 
  rename(Antimicrobial="ab_name") %>% distinct() %>% 
  mutate(Antimicrobial = str_replace(Antimicrobial,"/","-"))

###Engineer ddummy variables for antimicrobial agents
recipethis <- recipe(~ab_name,data=abx)
dummies <- recipethis %>% step_dummy(ab_name) %>% prep(training = abx)
dummy_data <- bake(dummies,new_data = NULL)
abx <- abx %>% cbind(dummy_data) %>% tibble()
abx <- abx %>% mutate(ab_name_Ampicillin = 
                        case_when(ab_name=="Ampicillin" ~
                                    1, TRUE ~ 0))

##Marrow suppression-related labels

###Check for leukopenia
print(d_labitems %>% filter(grepl("white",label,ignore.case=T)),n=25)
wbcs <- labevents %>% filter(itemid==51301) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

wbcs <- wbcs %>% new_lowvalue(new_leukopenia)
abx <- abx %>% abnormal_label(wbcs,leukopenia,new_leukopenia,ab_name)
ur_util <- ur_util %>% abnormal_label(wbcs,leukopenia,new_leukopenia,AMP)

###Check for anaemia
print(d_labitems %>% filter(grepl("Hemoglobin",label,ignore.case=T)),n=25)
hbs <- labevents %>% filter(itemid==51222) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()
hbs <- hbs %>% new_lowvalue(new_anaemia)
abx <- abx %>% abnormal_label(hbs,anaemia,new_anaemia,ab_name)
ur_util <- ur_util %>% abnormal_label(hbs,anaemia,new_anaemia,AMP)

###Check for thrombocytopenia
print(d_labitems %>% filter(grepl("Platelet",label,ignore.case=T)),n=25)
plts <- labevents %>% filter(itemid==51265) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()
plts <- plts %>% new_lowvalue(new_thrombocytopenia)
abx <- abx %>% abnormal_label(plts,thrombocytopenia,new_thrombocytopenia,ab_name)
ur_util <- ur_util %>% abnormal_label(plts,thrombocytopenia,new_thrombocytopenia,AMP)

###Check for previous bleeding
bleed_key <- d_icd_diagnoses %>% filter(grepl("bleed",long_title,ignore.case=T) &
                                          !grepl("without",long_title,ignore.case=T))
bleeding <- diagnoses %>% semi_join(bleed_key) %>% distinct(hadm_id,.keep_all = T) %>% 
  mutate(Bleeding_diagnosis = TRUE) %>% select(hadm_id,Bleeding_diagnosis)
abx <- abx %>% bleed_join()
ur_util <- ur_util %>% bleed_join()

###Check for other cytotoxic agents
cytotoxins <- c("adalimumab","aldesleukin","alemtuzumab",
                "amsacrine","arsenic trioxide","asparaginase",
                "axitinib","azacitidine","azathioprine",
                "belatacept","bendamustine","bevacizumab",
                "bexarotene","bleomycin","blinatumomab",
                "bortezomib","bosutinib","brentuximab vedotin",
                "busulfan","cabazitaxel","cabozantinib",
                "canakinumab","capecitabine","carboplatin",
                "carfilzomib","carmustine","ceritinib",
                "certolizumab pegol","chlorambucil","cisplatin",
                "cladribine","clofarabine","crisantaspase",
                "cyclophosphamide","cytarabine","dacarbazine",
                "dactinomycin","daratumumab","dasatinib",
                "daunorubicin","decitabine","dexrazoxane",
                "dinutuximab","docetaxel","doxorubicin",
                "epirubicin","eribulin","estramustine",
                "etoposide","fludarabine","fluorouracil",
                "ganciclovir","gemcitabine","gemtuzumab ozogamicin",
                "golimumab","hydroxycarbamide","ibrutinib",
                "idarubicin","ifosfamide","imatinib",
                "infliximab","inotuzumab ozogamicin","ipilimumab",
                "irinotecan","leflunomide","lenalidomide",
                "lomustine","melphalan","mercaptopurine",
                "methotrexate","mifamurtide","mitomycin",
                "mitotane","mitoxantrone","mogamulizumab",
                "nelarabine","nilotinib","niraparib",
                "nivolumab","obinutuzumab","olaparib",
                "oxaliplatin","paclitaxel","palbociclib",
                "panobinostat","pegaspargase","peginterferon alfa",
                "pembrolizumab","pemetrexed","pentostatin",
                "pixantrone","pomalidomide","procarbazine",
                "raltitrexed","ramucirumab","regorafenib",
                "ribociclib","rituximab","ropeginterferon alfa",
                "rucaparib","ruxolitinib","sorafenib",
                "streptozocin","sulfasalazine","sunitinib",
                "talazoparib","tegafur","temozolomide",
                "temsirolimus","thalidomide","thiotepa",
                "tioguanine","topotecan","trabectedin",
                "trastuzumab","trastuzumab","deruxtecan","trastuzumab","emtansine",
                "treosulfan","valganciclovir","vinblastine",
                "vincristine","vindesine","vinorelbine")
cytotoxins <- str_to_title(str_to_lower(cytotoxins))
cytotoxics_key <- drugs %>%
  filter(drug %in% cytotoxins) %>% distinct(hadm_id) %>% 
  mutate(Cytotoxic_agent = TRUE)
abx <- abx %>% cytotoxic_join()
ur_util <- ur_util %>% cytotoxic_join()

###Check for marrow suppression adjusted for bleeding & other cytotoxins
abx <- abx %>% marrow_check()
ur_util <- ur_util %>% marrow_check()

##Hepatotoxicity-related labels

###Check for elevated alkaline phosphatase
print(d_labitems %>% filter(grepl("Alkaline",label,ignore.case=T)),n=25)
alps <- labevents %>% filter(itemid==50863) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()
alps <- alps %>% new_highvalue(new_high_alp)
abx <- abx %>% abnormal_label(alps,high_alp,new_high_alp,ab_name)
ur_util <- ur_util %>% abnormal_label(alps,high_alp,new_high_alp,AMP)

###Check for elevated alanine aminotransferase
print(d_labitems %>% filter(grepl("transferase",label,ignore.case=T)),n=25)
alts <- labevents %>% filter(itemid==50861) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()
alts <- alts %>% new_highvalue(new_high_alt)
abx <- abx %>% abnormal_label(alts,high_alt,new_high_alt,ab_name)
ur_util <- ur_util %>% abnormal_label(alts,high_alt,new_high_alt,AMP)

###Check for elevated aspartate aminotransferase
print(d_labitems %>% filter(grepl("transferase",label,ignore.case=T)),n=25)
asts <- labevents %>% filter(itemid==50878) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()
asts <- asts %>% new_highvalue(new_high_ast)
abx <- abx %>% abnormal_label(asts,high_ast,new_high_ast,ab_name)
ur_util <- ur_util %>% abnormal_label(asts,high_ast,new_high_ast,AMP)

###Check for recent biliary procedures
biliary_proc_key <- d_icd_procedures %>% filter(grepl("Biliary",long_title,ignore.case=T) &
                                                  (grepl("Removal",long_title,ignore.case=T) |
                                                     grepl("Insertion",long_title,ignore.case=T) |
                                                     grepl("Inspection",long_title,ignore.case=T) |
                                                     grepl("Revision",long_title,ignore.case=T) |
                                                     grepl("Repair",long_title,ignore.case=T) |
                                                     grepl("Transplant",long_title,ignore.case=T) |
                                                     grepl("Introduction",long_title,ignore.case=T) |
                                                     grepl("Irrigation",long_title,ignore.case=T) |
                                                     grepl("Measurement",long_title,ignore.case=T) |
                                                     grepl("Endoscop",long_title,ignore.case=T) |
                                                     grepl("Fluoroscopy",long_title,ignore.case=T) |
                                                     grepl("Radiation",long_title,ignore.case=T)))
biliary <- procedures %>% semi_join(biliary_proc_key) %>% 
  distinct(hadm_id) %>% mutate(Biliary_procedure = TRUE)
abx <- abx %>% bil_proc_join()
ur_util <- ur_util %>% bil_proc_join()

###Check for hepatotoxicity adjusted for recent biliary procedures and chronic liver disease
abx <- abx %>% lft_check()
ur_util <- ur_util %>% lft_check()

###Check for unspecified adverse reactions
allerg_key <- d_icd_diagnoses %>% filter(grepl("antibiot",long_title,ignore.case=T)) %>%
  filter(grepl("adverse",long_title,ignore.case=T)) %>% select(icd_code,icd_version)
allerg_diags <- diagnoses %>% semi_join(allerg_key) %>% distinct(hadm_id) %>% 
  mutate(abx_ae = TRUE)
abx <- abx %>% left_join(allerg_diags) %>% mutate(abx_ae=case_when(abx_ae=TRUE~abx_ae,TRUE~FALSE))
ur_util <- ur_util %>% left_join(allerg_diags) %>% mutate(abx_ae=case_when(abx_ae=TRUE~abx_ae,TRUE~FALSE))

##Check for overall toxicity

abx <- abx %>% toxicity_check()
ur_util <- ur_util %>% toxicity_check()

##Sepsis adverse outcome-related labels

###Sepsis diagnosis on admission
sepsis_key <- d_icd_diagnoses %>% filter(grepl("SEPSIS",long_title,ignore.case=T))
sepsis_key <- diagnoses %>% semi_join(sepsis_key) %>% distinct(hadm_id) %>% 
  mutate(admission_sepsis=TRUE)
abx <- abx %>% left_join(sepsis_key)
ur_util <- ur_util %>% left_join(sepsis_key)

###Infection diagnosis on admission
infection_key <- d_icd_diagnoses %>% filter(grepl("infection",long_title,ignore.case=T))
infection_key <- diagnoses %>% semi_join(infection_key) %>% distinct(hadm_id) %>% 
  mutate(admission_infection=TRUE)
abx <- abx %>% left_join(infection_key)
ur_util <- ur_util %>% left_join(infection_key)

###High obs frequency in ED
obfreq_key <- urines_ref %>% select(hadm_id,ob_freq) %>% distinct(hadm_id,.keep_all = T)
abx <- abx %>% left_join(obfreq_key) 
abx <- abx %>% mutate(high_obfreq = case_when(ob_freq > median(abx$ob_freq,na.rm=T)~TRUE,
                                 TRUE~FALSE))
ur_util <- ur_util %>% mutate(high_obfreq = case_when(ob_freq > median(ur_util$ob_freq,na.rm=T)~TRUE,
                                                      TRUE~FALSE))

###High CRP on admission
highcrp_key <- urines_ref %>% select(hadm_id,highCRP) %>% distinct(hadm_id,.keep_all = T)
abx <- abx %>% left_join(highcrp_key) %>% mutate(highCRP=case_when(highCRP~highCRP,
                                                  TRUE~FALSE))

###Sepsis high-risk based on ED observations themselves
staykey <- edstays %>% select(hadm_id,stay_id) %>% 
  distinct(hadm_id,.keep_all = T) %>% filter(!is.na(hadm_id))
triage <- triage %>% left_join(staykey)
vitalsign <- vitalsign %>% left_join(staykey)
triagekey <- triage %>% keymaker()
vitalskey <- vitalsign %>% mutate(chiefcomplaint=NA) %>% keymaker()
ed_obskey <- tibble(rbind(triagekey,vitalskey))
ed_obskey <- ed_obskey %>% mutate(temperature = ((temperature -32)*5)/9 )
ed_obskey <- ed_obskey %>% pc_dummies()
ed_obskey <- ed_obskey %>% mutate(SIRS = case_when(
  pc_confusion==TRUE |
    resprate > 20 |
    sbp < 100 |
    heartrate > 90 |
    temperature < 36 ~ TRUE, TRUE~FALSE
))
ed_obskey <- ed_obskey %>% group_by(hadm_id) %>% 
  mutate(SIRS=max(SIRS),
         SIRS=case_when(SIRS==1~ TRUE,TRUE~FALSE)) %>% ungroup() %>% 
  distinct(hadm_id,.keep_all=T)

abx <- abx %>% left_join(ed_obskey)
ur_util <- ur_util %>% left_join(ed_obskey) %>% 
  mutate(SIRS=case_when(SIRS~SIRS,TRUE~FALSE))

###Update sepsis criterion with above variables
abx <- abx %>% update_sepsis()
ur_util <- ur_util %>% update_sepsis()

###Check for death in the next 28 days
ur_util$ob_freq %>% sort(T)
abx <- abx %>% death_check(hadm,death_28d,ab_name)
ur_util <- ur_util %>% death_check(hadm,death_28d,AMP)

###Check for ICU admission in the next 28 days
abx <- abx %>% icu_check(icu,icu_28d,ab_name)
ur_util <- ur_util %>% icu_check(icu,icu_28d,AMP)

###Check if still an inpatient 7 days later
abx <- abx %>% inpt7d_check()
ur_util <- ur_util %>% inpt7d_check()

###Check for overall sepsis adverse outcomes
abx <- abx %>% sepsis_ae_check()
ur_util <- ur_util %>% sepsis_ae_check()

##Combinations dataframe

###Remove drugs started twice in the same day
combos <- abx %>% mutate(startdate=as.Date(starttime)) %>% 
  arrange(subject_id,starttime) %>% 
  distinct(subject_id,ab_name,startdate,.keep_all = T)
combos_ref <- abx_ref %>% mutate(startdate=as.Date(starttime)) %>% 
  arrange(subject_id,starttime) %>% 
  distinct(subject_id,ab_name,startdate,.keep_all = T)

###Check for combinations
combos <- combocheck(combos)
combos_ref <- combocheck(combos_ref)
combos <- combodrugcheck(combos) %>% filter(combo_agent)
combos_ref <- combodrugcheck(combos_ref) %>% filter(combo_agent)
write_csv(combos,"combos.csv")

###Prior agent start and stop times
combos <- combos %>% prior_agent_stoptimes()
combos_ref <- combos_ref %>% prior_agent_stoptimes()
combos <- combos %>% starttime_check()
combos <- combos %>% mutate(combo_starttime=as_datetime(combo_starttime))
combos_ref <- combos_ref %>% starttime_check()
combos_ref <- combos_ref %>% mutate(combo_starttime=as_datetime(combo_starttime))

###Stoptime of current regime
combos <- combos %>% curr_regime_stoptimes()
combos_ref <- combos_ref %>% curr_regime_stoptimes()
combos <- combos %>% stoptime_check()
combos <- combos %>% mutate(combo_stoptime=as_datetime(combo_stoptime))
combos_ref <- combos_ref %>% stoptime_check()
combos_ref <- combos_ref %>% mutate(combo_stoptime=as_datetime(combo_stoptime))

###Check if second agent of combination
combos <- combos %>% combo_booleans()
combos_ref <- combos_ref %>% combo_booleans()

###Count number of drugs in combination
combos <- combos %>% count_drugs()
combos_ref <- combos_ref %>% count_drugs()
start_col <- "Ampicillin"
end_col <- "Vancomycin"

###Replace name, starttimes, and stoptimes with combination values
combos <- combos %>% combo_mutate()
combos_ref <- combos_ref %>% combo_mutate()
write_csv(combos,"combos.csv")
combos <- combos %>% times_amend()
combos_ref <- combos_ref %>% times_amend()

##Preprocessing for utility function

###Remove same-day antimicrobial startdate duplicates
abx <- abx %>% dup_remove()
abx_ref <- abx_ref %>% dup_remove()

###Amend abx dataframe to include combinations
abx <- abx %>% bind_combos(combos)
abx_ref <- abx_ref %>% bind_combos(combos_ref)

###Attach row ids
abx <- abx %>% rowids()
abx_ref <- abx_ref %>% rowids()

###Save interim with all antimicrobial combinations, then filter to 2 or fewer
write_csv(abx,"abx_all_combos.csv")
write_csv(abx_ref,"abx_ref_all_combos.csv")
abx <- abx %>% filter(str_count(ab_name,"_") <2)

###Filter out rare combinations (<500 entries)
common_combos <- abx %>% count(ab_name) %>% filter(n>=500) %>% pull(ab_name)
abx <- abx %>% filter(ab_name%in%common_combos)

###Filter out combinations that last less than 24 hours
abx <- abx %>% mutate(combo_24h = case_when(as.numeric(stoptime-starttime) > 1~TRUE,
                      TRUE~FALSE)) %>% filter(combo_24h)

###Re-engineer antimicrobial dummy variables#
abx <- abx %>% select(-(ab_name_Ampicillin.sulbactam:ab_name_Ampicillin))
recipethis <- recipe(~ab_name,data=abx)
dummies <- recipethis %>% step_dummy(ab_name) %>% prep(training = abx)
dummy_data <- bake(dummies,new_data = NULL)
abx <- abx %>% cbind(dummy_data) %>% tibble()
abx <- abx %>% mutate(ab_name_Ampicillin = 
                        case_when(ab_name=="Ampicillin" ~
                                    1, TRUE ~ 0))

##Antimicrobial-related variables in ur_util

###Check for currently prescribed antimicrobial agent
write_csv(ur_util,"extra_ur_util.csv")
ab_key <- abx_ref %>% select(subject_id,ab_name,starttime,stoptime)
ur_util <- ur_util %>% select(-ab_name)
urine_abx <- ur_util %>% left_join(ab_key,by=c("subject_id")) %>% 
  mutate(on_ab = case_when(
    storetime > starttime & storetime < stoptime ~ TRUE,
    TRUE ~ FALSE )) %>% filter(on_ab)
urine_abx <- urine_abx %>%
  group_by(micro_specimen_id) %>%
  summarise(ab_name = paste(ab_name, collapse = ", "), .groups = 'drop') %>% 
  ungroup() %>% mutate(on_ab=TRUE)
ur_util <- ur_util %>% left_join(urine_abx)

###Apply AST R result value based on currently-prescribed antimicrobial agent
ur_util <- ur_util %>% mutate(
  AMP_R_value = case_when(on_ab & grepl("(Amoxicillin|Benzylpenicillin|Ampicillin)",ab_name) 
                          & !grepl("(clavulanic|sulbactam)",ab_name) ~ TRUE, TRUE~FALSE),
  SAM_R_value = case_when(on_ab & grepl("(Amoxicillin|Benzylpenicillin|Ampicillin)",ab_name)
                          & !grepl("clavulanic",ab_name) ~ TRUE, TRUE~FALSE),
  TZP_R_value = case_when(on_ab & grepl("(Amoxicillin|Benzylpenicillin|Ampicillin|Piperacillin)",ab_name)
                          ~ TRUE, TRUE~FALSE),
  CZO_R_value = case_when(on_ab & grepl("Cefazolin",ab_name)
                          ~ TRUE, TRUE~FALSE),
  CRO_R_value = case_when(on_ab &((org_order=="Enterobacterales" &
                                     grepl("(Cefazolin|Cefalexin)",ab_name))|
                                    grepl("Ceftriaxone",ab_name) |
                                    (org_order=="Staphylococcus" &
                                       grepl("(Amoxicillin|Benzylpenicillin|Ampicillin|Piperacillin)",ab_name)))
                          ~ TRUE, TRUE~FALSE),
  CAZ_R_value = case_when(on_ab &((org_order=="Enterobacterales" &
                                     grepl("(Cefazolin|Cefalexin)",ab_name))|
                                    grepl("Ceftazidime",ab_name))
                          ~ TRUE, TRUE~FALSE),
  FEP_R_value = case_when(on_ab &((org_order=="Enterobacterales" &
                                     grepl("(Cefazolin|Cefalexin)",ab_name))|
                                    grepl("Cefepime",ab_name))
                          ~ TRUE, TRUE~FALSE),
  MEM_R_value = case_when(on_ab & grepl("(Meropenem|Ertapenem|Amoxicillin|Benzylpenicillin|Ampicillin|Piperacillin)",ab_name)
                          ~ TRUE, TRUE~FALSE),
  CIP_R_value = case_when(on_ab &((org_order=="Enterobacterales" &
                                     grepl("(Levofloxacin|Moxifloxacin)",ab_name))|
                                    grepl("Ciprofloxacin",ab_name))
                          ~ TRUE, TRUE~FALSE),
  GEN_R_value = case_when(on_ab & grepl("Gentamicin",ab_name)
                          ~ TRUE, TRUE~FALSE),
  SXT_R_value = case_when(on_ab & grepl("Trimethoprim",ab_name)
                          ~ TRUE, TRUE~FALSE),
  NIT_R_value = case_when(on_ab & grepl("Nitrofurantoin",ab_name)
                          ~ TRUE, TRUE~FALSE),
  VAN_R_value = case_when(on_ab & grepl("Vancomycin",ab_name)
                          ~ TRUE, TRUE~FALSE),
)

###Fixing continuous and numeric variables in abx dataframe
abx <- abx %>% mutate(resprate=case_when(resprate>30|resprate<12~
                        NA,TRUE~resprate))
abx <- abx %>% mutate(sbp=case_when(sbp>250|sbp<50~
                                           NA,TRUE~sbp))
abx <- abx %>% mutate(dbp=case_when(dbp>120|dbp<30~
                                           NA,TRUE~dbp))
abx <- abx %>% mutate(heartrate=case_when(heartrate>250|heartrate<30~
                                      NA,TRUE~heartrate))
abx <- abx %>% mutate(temperature=case_when(temperature>41|temperature<34~
                                            NA,TRUE~temperature))
abx <- abx %>% mutate(o2sat=case_when(o2sat>100|o2sat<70~
                                              NA,TRUE~o2sat))
abx <- abx %>% mutate(heartrate=case_when(is.na(heartrate)~mean(abx$heartrate,na.rm=T),
                                                        TRUE~heartrate),
                                    resprate=case_when(is.na(resprate)~mean(abx$resprate,na.rm=T),
                                                       TRUE~resprate),
                                    sbp=case_when(is.na(sbp)~mean(abx$sbp,na.rm=T),
                                                  TRUE~sbp),
                                    dbp=case_when(is.na(dbp)~mean(abx$dbp,na.rm=T),
                                                  TRUE~dbp),
                                    o2sat=case_when(is.na(o2sat)~mean(abx$o2sat,na.rm=T),
                                                    TRUE~o2sat),
                                    temperature=case_when(is.na(temperature)~mean(abx$temperature,na.rm=T),
                                                          TRUE~temperature),
)
abx <- abx %>% mutate(heartrate=standardize(heartrate),
                                    resprate=standardize(heartrate),
                                    sbp=standardize(heartrate),
                                    dbp=standardize(heartrate),
                                    acuity=standardize(heartrate),
                                    o2sat=standardize(heartrate),
                                    temperature=standardize(temperature))
abx <- abx %>% mutate(ob_freq = case_when(is.na(ob_freq) ~ mean(abx$ob_freq,na.rm=T),
                                          TRUE ~ ob_freq))

###Write interim CSVs
write_csv(abx,"interim_abx.csv")
write_csv(ur_util,"interim_ur_util.csv")
write_csv(combos,"combos.csv")

###Split abx dataframe into train and test dataframes
subjects <- abx %>% distinct(subject_id)
smp_size <- floor(0.8 * nrow(subjects))
set.seed(123)
train_ind <- sample(seq_len(nrow(subjects)), size = smp_size)
train_ids <- subjects[train_ind,]
test_ids <- subjects[-train_ind,]
abx <- abx %>% mutate(CDI = factor(CDI),
                      overall_tox = factor(overall_tox),
                      sepsis_ae=factor(sepsis_ae))
train_abx <- abx %>% semi_join(train_ids,by="subject_id")
test_abx <- abx %>% semi_join(test_ids,by="subject_id")

write_csv(train_abx,"train_abx.csv")
write_csv(test_abx,"test_abx.csv")


