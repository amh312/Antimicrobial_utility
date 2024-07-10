library("brms")
library("tidyverse")
library("tidymodels")
library("openxlsx")
library("mlogit")
library("AMR")
library("glue")


#####FUNCTIONS#################################################################

#Read-in and cleaning
read_in <- function(file_name) {
  
  file_path <- file.path(path_to_data, file_name)
  df <- fread(file_path)
  
}

micro_clean <- function(file_location,file_name) {
  
  path_to_data <- file_location
  
  read_in <- function(file_name) {
    
    file_path <- file.path(path_to_data, file_name)
    df <- fread(file_path)
    
  }
  
  df <- read_in(file_name)
  
  df <- df %>% mutate(org_name = str_remove(org_name,"POSITIVE FOR"),#---------------------------------------------------Cleaning
                      org_name = str_remove(org_name,"PRESUMPTIVELY"),
                      org_name = str_remove(org_name,"PRESUMPTIVE"),
                      org_name = str_remove(org_name,"PROBABLE"),
                      org_name = str_remove(org_name,"IDENTIFICATION"),
                      org_name = str_remove(org_name,"RESEMBLING"),
                      org_name = str_remove(org_name,"SEEN"),
                      org_name = str_remove(org_name,"MODERATE"),
                      org_name = str_remove(org_name,"FEW"),
                      org_name = str_remove(org_name,"BETA"),
                      org_name = str_remove(org_name,"METHICILLIN RESISTANT"),
                      org_name = str_remove(org_name,"NUTRITIONALLY VARIANT"),
                      org_name = str_remove(org_name,"NOT C. PERFRINGENS OR C. SEPTICUM"),
                      org_name = str_remove(org_name,"-LACTAMASE POSITIVE"),
                      org_name = str_remove(org_name,"-LACTAMASE NEGATIVE"),
                      org_name = str_remove(org_name,"VIRAL ANTIGEN"),
                      org_name = str_remove(org_name,"CANDIDA INCONSPICUA"),
                      org_name = str_remove(org_name,"/POSADASII"),
                      org_name = str_remove(org_name,"NOT FUMIGATUS, FLAVUS OR NIGER"),
                      org_name = str_remove(org_name,"MRSA NEGATIVE"),
                      org_name = str_remove(org_name,"HISTOLYTICA/DISPAR"),
                      org_name = case_when(grepl("NON-FERMENTER",org_name)~"PSEUDOMONADALES",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ABIOTROPHIA/GRANULICATELLA",org_name)~"STREPTOCOCCUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("S. AUREUS POSITIVE",org_name)~"STAPHYLOCOCCUS AUREUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("ASPERGILLUS FUMIGATUS COMPLEX",org_name)~"ASPERGILLUS FUMIGATUS",                     
                                           TRUE~org_name),
                      org_name = case_when(grepl("(CRYPTOSPORIDIUM PARVUM OOCYSTS|CUNNINGHAMELLA BERTHOLLETIAE|EPIDERMOPHYTON FLOCCOSUM|EXOPHIALA JEANSELMEI COMPLEX|SCEDOSPORIUM|NEOASCOCHYTA DESMAZIERI|NEOSCYTALIDIUM DIMIDIATUM|LOMENTOSPORA|NEUROSPORA|PERONEUTYPA SCOPARIA|SPOROTHRIX SCHENCKII COMPLEX|ZYGOSACCHAROMYCES FERMENTATI)",org_name)~"UNKNOWN FUNGUS",                     
                                           TRUE~org_name)
  ) %>%#-------------------------------------------------------------------------------------------------------------------Removal of AMR package non-interpretable rows
    filter(!grepl("(CANCELLED|VIRUS|SIMPLEX|PARAINFLUENZA|INFLUENZA A|INFLUENZA B|TICK|AFB GROWN|GRAM VARIABLE RODS|HYMENOLEPIS)",org_name)) %>% 
    mutate(ab_name=AMR::as.ab(ab_name)) %>%#-------------------------------------------------------------------------------AMR package parsing of antimicrobial and organism names
    mutate(org_name=AMR::as.mo(org_name)) %>% 
    MIMER::transpose_microbioevents(
      key_columns = c('subject_id','micro_specimen_id','isolate_num','org_name','ab_itemid','test_name','test_seq'),#------Transpose AST results to columns
      required_columns =c('subject_id','chartdate',"hadm_id","order_provider_id",
                          "charttime","micro_specimen_id","spec_itemid","spec_type_desc",
                          "storedate","storetime","test_itemid","test_name","org_itemid",
                          "isolate_num","org_name","comments",'test_seq'),
      transpose_key_column = 'ab_name',
      transpose_value_column = 'interpretation',
      fill = "N/A",
      non_empty_filter_column='subject_id') %>%
    add_column(AMX=NA, AMC=NA, TIC=NA,PME=NA, FOS=NA,TMP=NA,#---------------------------------------------------------------Add missing AMR package-recognised antimicrobial agent columns
               MFX=NA, NOR=NA,CPD = NA, FOX1=NA,TEC=NA,TLV=NA,ORI=NA,
               TGC=NA,AZM=NA,ATM=NA,CRB=NA,CTX=NA,CPT=NA,SPT=NA,TZD=NA,ERV=NA,OMC=NA,FDX=NA,
               CZT=NA,LEX=NA,CLR=NA,DAL=NA,CZA=NA,NOV=NA,ETP=NA,
               MTR=NA,QDA=NA,TEM=NA,COL=NA,CHL=NA,BPR=NA,CEC=NA) %>%
    mutate(org_fullname = AMR::mo_fullname(org_name),#----------------------------------------------------------------------Additional organism categorisation columns
           org_kingdom = AMR::mo_kingdom(org_name),
           org_phylum = AMR::mo_phylum(org_name),
           org_class = AMR::mo_class(org_name),
           org_order = AMR::mo_order(org_name),
           org_family = AMR::mo_family(org_name),
           org_genus = AMR::mo_genus(org_name),
           org_species = AMR::mo_species(org_name),
           org_gram = AMR::mo_gramstain(org_name),
           org_o2 = AMR::mo_oxygen_tolerance(org_name),
           org_path = AMR::mo_pathogenicity(org_name),
           UTI = case_when(grepl("URINE",spec_type_desc)~TRUE,
                           TRUE~FALSE)) %>%
    relocate(PEN,OXA,AMP,AMX,PIP,TIC,CRB,PME,SAM,AMC,TZP,TEM,#--------------------------------------------------------------AST column reorganisation
             ATM,
             LEX,CZO,CEC,CXM,FOX1,CTX,CRO,CAZ,CPD,FEP,CPT,BPR,CZA,CZT,
             ETP,MEM,IPM,
             LVX,MFX,CIP,NOR,
             GEN,TOB,SPT,
             TMP,SXT,
             COL,NIT,FOS,NOV,CHL,
             TGC,ERV,OMC,TCY,
             ERY,CLR,AZM,CLI,QDA,
             LNZ,TZD,TEC,VAN,DAL,TLV,ORI,DAP,RIF,
             FDX,MTR,
             .before = "AMK"
    ) %>% relocate(AMK,.after = "GEN") %>% 
    mutate_at(vars(PEN:MTR),as.sir)
  
  df %>% mutate(#----------------------------------------------------------------------------------------------------------Addition of breakpoint interpretation and UTI columns
    guideline=rep("Original CLSI",nrow(df)),
    urine_interp = case_when(spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               (org_path=="Potentially pathogenic" |
                                  grepl("(discontinued|MIXED)",comments))~ "Possible UTI",
                             spec_type_desc=="URINE" &
                               !is.na(org_name) &
                               org_path=="Pathogenic" &
                               comments=="" ~ "Probable UTI",
                             TRUE ~ "Unlikely UTI"),
    AMPC=case_when(grepl("Citrobacter braakii",org_fullname) |#-----------------------------------------------------------Addition of chromosomal AmpC column
                     grepl("Citrobacter freundii",org_fullname) |
                     grepl("Citrobacter gillenii",org_fullname) |
                     grepl("Citrobacter murliniae",org_fullname) |
                     grepl("Citrobacter rodenticum",org_fullname) |
                     grepl("Citrobacter sedlakii",org_fullname) |
                     grepl("Citrobacter werkmanii",org_fullname) |
                     grepl("Citrobacter youngae",org_fullname) |
                     grepl("Enterobacter",org_fullname) |
                     grepl("Hafnia alvei",org_fullname) |
                     grepl("Klebsiella aerogenes",org_fullname) |
                     grepl("Morganella morganii",org_fullname) |
                     grepl("Providencia",org_fullname) |
                     grepl("Serratia",org_fullname) |
                     org_order=="Enterobacterales"& (CAZ=="R"|CAZ=="I")&FEP=="S"~"R",
                   (org_order=="Enterobacterales"& (CAZ=="R" & is.na(FEP))) |
                     (org_order=="Enterobacterales"& (is.na(CAZ) & FEP=="S" )) ~ NA,
                   TRUE~"S")
  ) %>% relocate(comments,.before="guideline")
  
}

prescriptions_clean <- function(file_location,file_name) {
  
  path_to_data <- file_location
  
  read_in <- function(file_name) {
    
    file_path <- file.path(path_to_data, file_name)
    df <- fread(file_path)
    
  }
  
  df <- read_in(file_name)
  
  df <- df %>%
    MIMER::clean_antibiotics(df,drug_col=drug)
  
}

#Intrinsic resistance population
intr_mic <- function(df) {
  
  x <- custom_eucast_rules(genus=="Enterococcus"~cephalosporins=="R",#----------------------------------------------Add custom rules
                           genus=="Enterococcus"~aminoglycosides=="R",
                           genus=="Enterococcus"~macrolides=="R",
                           genus=="Enterococcus"~lincosamides=="R",
                           fullname=="Enterococcus faecium"~carbapenems=="R",
                           genus=="Enterococcus"&AMP=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~SAM=="S",
                           genus=="Staphylococcus"&OXA=="S"~TZP=="S",
                           genus=="Staphylococcus"&OXA=="S"~AMC=="S",
                           genus=="Staphylococcus"&OXA=="S"~cephalosporins=="S",
                           genus=="Staphylococcus"&OXA=="S"~carbapenems=="S",
                           genus=="Staphylococcus"~CAZ=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~SAM=="R",
                           genus=="Staphylococcus"&OXA=="R"~TZP=="R",
                           genus=="Staphylococcus"&OXA=="R"~AMC=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_1st=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_2nd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_3rd=="R",
                           genus=="Staphylococcus"&OXA=="R"~cephalosporins_4th=="R",
                           genus=="Staphylococcus"&OXA=="R"~carbapenems=="R",
                           genus=="Staphylococcus"&PEN=="S"~aminopenicillins=="S",
                           genus=="Staphylococcus"&PEN=="R"~aminopenicillins=="R",
                           genus=="Streptococcus"&PEN=="S"~aminopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~ureidopenicillins=="S",
                           genus=="Streptococcus"&PEN=="S"~cephalosporins_except_caz=="S",
                           kingdom=="Fungi"~aminoglycosides=="R",
                           kingdom=="Fungi"~antimycobacterials=="R",
                           kingdom=="Fungi"~betalactams=="R",
                           kingdom=="Fungi"~quinolones=="R",
                           kingdom=="Fungi"~lincosamides=="R",
                           kingdom=="Fungi"~macrolides=="R",
                           kingdom=="Fungi"~oxazolidinones=="R",
                           kingdom=="Fungi"~polymyxins=="R",
                           kingdom=="Fungi"~streptogramins=="R",
                           kingdom=="Fungi"~tetracyclines=="R",
                           kingdom=="Fungi"~trimethoprims=="R",
                           kingdom=="Fungi"~glycopeptides=="R",
                           kingdom=="Fungi"~MTR=="R",
                           kingdom=="Fungi"~FDX=="R",
                           kingdom=="Fungi"~NIT=="R",
                           kingdom=="Fungi"~FOS=="R",
                           kingdom=="Fungi"~NOV=="R",
                           kingdom=="Fungi"~CHL=="R",
                           kingdom=="Fungi"~DAP=="R",
                           genus=="Enterobacter"~AMP=="R",
                           genus=="Enterobacter"~AMX=="R",
                           genus=="Enterobacter"~SAM=="R",
                           genus=="Enterobacter"~AMC=="R",
                           genus=="Enterobacter"~LEX=="R",
                           genus=="Enterobacter"~CZO=="R",
                           PIP=="S"~TZP=="S",
                           TMP=="S"~SXT=="S",
                           PEN=="S"~AMP=="S",
                           AMP=="S"~SAM=="S",
                           SAM=="S"~TZP=="S",
                           TZP=="R"~SAM=="R",
                           SAM=="R"~AMP=="R",
                           AMP=="R"~PEN=="R",
                           phylum=="Pseudomonadota"~DAP=="R",
                           genus=="Pseudomonas"~NIT=="R",
                           genus=="Pseudomonas"~SXT=="R",
                           org_o2=="aerobe"~MTR=="R")
  
  df %>% #---------------------------------------------------------------------------------------------------------Fill intrinsic resistance
    eucast_rules(col_mo = "org_name",
                 ampc_cephalosporin_resistance = "R",
                 rules="all",
                 custom_rules = x) %>% 
    mutate(MTR=case_when(org_o2!="anaerobe"~"R",
                         TRUE~MTR),
           PEN=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~PEN),
           AMP=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         grepl("Staphylococcus",org_fullname) & PEN=="R"~"R",
                         TRUE~AMP),
           AMX=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~AMX),
           AMC=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~AMC),
           PIP=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~PIP),
           SAM=case_when(grepl("Staphylococcus",org_fullname) & OXA=="R"~"R",
                         TRUE~SAM),
           cleaning = rep("(w/intr. R)",nrow(df)))
  
}

#Imputing missing results
res_sim <- function(df,col,condition,col2,condition2,antibiotic,alpha_prior,beta_prior,antimicrobial_name,extra="") {
  
  antibiotic <- enquo(antibiotic)
  col <- enquo(col)
  col2 <- enquo(col2)
  
  df$isolate_id <- as.character(df$org_name)#------------------------------------------------------------------------Unique isolate id column
  df$isolate_id[!is.na(df$isolate_id)] <- 1:sum(!is.na(df$org_name))
  
  x <- nrow(df %>%#--------------------------------------------------------------------------------------------------Number of observed Rs
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !!antibiotic=="R"))
  
  N <- nrow(df %>%#--------------------------------------------------------------------------------------------------Total number of observed results
              dplyr::filter(grepl(condition,!!col) &
                              grepl(condition2,!!col2) &
                              !is.na(!!antibiotic)))
  
  if(N>1) {
    
    #Bayesian calculations
    p <- seq( from=0 , to=1 , length.out=1e4 )
    
    posterior_alpha <- alpha_prior + x
    
    posterior_beta <- beta_prior + N - x
    
    mean_prob <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_prob <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    
    prior <- (p ^ (alpha_prior - 1)) * ((1 - p) ^ (beta_prior - 1))
    
    likelihood <- (p ^ x) * ( (1 - p) ^ (N - x)  )
    
    posterior <-  p ^ (posterior_alpha - 1) * (1 - p) ^ (posterior_beta - 1)
    
    if (mean(posterior) != 0 & mean(posterior) != 1) {
      
      #Sampling posterior distribution
      
      prior_samples <- sample( p , prob=prior , size=1e4 , replace=TRUE )
      prior_samples <- tibble(Probability = prior_samples,Distribution=rep("Prior",length(prior_samples)))
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      likelihood_samples <- tibble(Probability = likelihood_samples,Distribution=rep("Likelihood",length(likelihood_samples)))
      
      post_samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      post_samples <- tibble(Probability = post_samples,Distribution=rep("Posterior",length(post_samples)))
      
      #Prior, likelihood and posterior density plot
      
      prior_2 <- prior/max(prior)
      prior_plot <- tibble(Density = prior_2,Distribution=rep("Prior",length(prior_2)),Probability=p)
      
      likelihood_2 <- likelihood/max(likelihood)
      likelihood_plot <- tibble(Density = likelihood_2,Distribution=rep("Likelihood",length(likelihood_2)),Probability=p)
      
      posterior_2 <- posterior/max(posterior)
      post_plot <- tibble(Density = posterior_2,Distribution=rep("Posterior",length(posterior_2)),Probability=p)
      
      post_df <- rbind(prior_plot,likelihood_plot,post_plot)
      post_df$Distribution <- factor(post_df$Distribution, levels=c("Prior","Likelihood","Posterior"))
      
      print(ggplot(post_df,aes(y=Density,x=Probability,group=Distribution,fill=Distribution,color=Distribution)) +
              geom_line() +
              theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) +
              geom_line(size=.5) +
              geom_ribbon(data=subset(post_df,Probability>0 & Probability<1),aes(x=Probability,ymax=Density),ymin=0,alpha=0.3) +
              scale_fill_manual(name='', values=c("Prior" = "red", "Likelihood" = "green4","Posterior"="blue")) +
              guides(color = FALSE) +
              labs(title=glue("Probability: {antimicrobial_name} resistance in {condition}{extra}")))
      
      N_star <- nrow(df %>%
                       dplyr::filter(grepl(condition,!!col) &
                                       grepl(condition2,!!col2) &
                                       is.na(!!antibiotic)))
      
      #Posterior predictive bar plot
      
      BetaBinom <- Vectorize(function(x_star){
        log.val <- lchoose(N_star, x_star) + lbeta(posterior_alpha+x_star,posterior_beta+N_star-x_star) - lbeta(posterior_alpha,posterior_beta)
        return(exp(log.val))
      })
      
      
      post_predictive <- BetaBinom(1:N_star)
      plot(1:N_star,BetaBinom(1:N_star),type="h",col="darkblue",xlab="Estimated prevalence of resistance",ylab="Probability density",
           main = glue("Estimated prevalence of {antimicrobial_name} resistance in {N_star} {condition}{extra} isolates"),cex.axis= 1.5,cex.lab=1.5,lwd=4)
      
      samples <- sample( p , prob=posterior , size=1e4 , replace=TRUE )
      
    } else {
      
      samples = mean(posterior)
      
    }
    
    #Summary statistics
    
    simu <- rbetabinom(nrow(df %>%
                              dplyr::filter(grepl(condition,!!col) &
                                              grepl(condition2,!!col2) &
                                              is.na(!!antibiotic))),
                       size = 1 ,
                       shape1 = posterior_alpha,
                       shape2 = posterior_beta)
    
    n_likelihood <- nrow(df %>%
                           dplyr::filter(grepl(condition,!!col) &
                                           grepl(condition2,!!col2) &
                                           !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    
    if(mean(likelihood)!=0) {
      
      likelihood_samples <- sample( p , prob=likelihood , size=1e4 , replace=TRUE )
      
    } else {
      
      likelihood_samples <- 0
      
    }
    
    
    
    prior_amr_rate <- alpha_prior/(alpha_prior+beta_prior)
    mean_likelihood <- mean(rbinom(1e4,
                                   size = 1 ,
                                   prob = likelihood_samples))
    mean_posterior <- posterior_alpha / (posterior_alpha + posterior_beta)
    mode_posterior <- (posterior_alpha - 1) / (posterior_alpha + posterior_beta - 2)
    HPDI_posterior <- ifelse(mean(samples)==0,NA,data.frame(t(HPDI(samples,prob=0.95))))
    HPDI_samples <- data.frame(cbind(n_likelihood,n_simulated,prior_amr_rate,mean_likelihood,mean_posterior,mode_posterior,HPDI_posterior))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    #Result simulation
    
    simu <- ifelse(simu==1,"R","S")
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic)) %>% data.frame()
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    sample_df <- tibble(cbind(tibble(samples=samples),tibble(org_name=rep(glue("{condition}"),length(samples)))))
    
    sample_df <-
      sample_df %>%
      group_by(org_name) %>%
      mutate(outlier = samples < quantile(samples, .25) - 1.5*IQR(samples) | samples > quantile(samples, .75) + 1.5*IQR(samples)) %>%
      ungroup
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_df"),sample_df,envir = .GlobalEnv)
    
    df %>% select(-isolate_id)
    
  } else {
    
    #Output for insufficient results to inform likelihood
    
    #Summary statistics
    
    simu <- rbetabinom(nrow(df %>%
                              dplyr::filter(grepl(condition,!!col) &
                                              grepl(condition2,!!col2) &
                                              is.na(!!antibiotic))),
                       size = 1 ,
                       shape1 = alpha_prior,
                       shape2 = beta_prior)
    
    n_likelihood <- nrow(df %>%
                           dplyr::filter(grepl(condition,!!col) &
                                           grepl(condition2,!!col2) &
                                           !is.na(!!antibiotic)))
    n_simulated <- nrow(df %>%
                          dplyr::filter(grepl(condition,!!col) &
                                          grepl(condition2,!!col2) &
                                          is.na(!!antibiotic)))
    
    
    #Result simulation
    
    simu <- ifelse(simu==1,"R","S")
    
    target <- df %>% dplyr::filter(grepl(condition,!!col) &
                                     grepl(condition2,!!col2)) %>% 
      select(isolate_id,!!antibiotic) %>% 
      mutate(!!antibiotic := as.character(!!antibiotic)) %>% data.frame()
    
    target[is.na(target)] <- simu
    
    target <- target %>% distinct(isolate_id,.keep_all = T)
    
    df <- df %>% mutate(!!antibiotic := as.character(!!antibiotic)) %>% 
      rows_update(target,by=c("isolate_id"))
    
    df %>% select(-isolate_id)
    
    print(glue("Insufficient results to calculate {antimicrobial_name} resistance likelihood for {condition}"))
    
    missing <- data.frame(Antimicrobial=glue("{antimicrobial_name}"),Organism=glue("{condition}"))
    
    assign(glue("missing"),missing,envir = .GlobalEnv)
    
    missings <- rbind(missings,missing)
    
    assign(glue("missings"),missings,envir = .GlobalEnv)
    
    samples <- rep(1,1e4)
    
    HPDI_samples <- data.frame(matrix(nrow=1,ncol=7))
    rownames(HPDI_samples) <- glue("{condition}_{antimicrobial_name}{extra}")
    colnames(HPDI_samples) <- c("n_likelihood","n_simulated","prior_amr_rate",
                                "mean_likelihood","mean_posterior","mode_posterior",
                                "HPDI_posterior")
    
    assign(glue("{condition}_{antimicrobial_name}{extra}"),HPDI_samples,envir = .GlobalEnv)
    
    sample_df <- tibble(cbind(tibble(samples=rep(5,length(samples)))),tibble(org_name=rep(glue("{condition}"),length(samples))),
                        tibble(outlier=rep(FALSE,length(samples))))
    
    assign(glue("{condition}_{antimicrobial_name}{extra}_df"),sample_df,envir = .GlobalEnv)
    
    
    df %>% select(-isolate_id)
    
  }
  
}

#Individual antimicrobial simulations
COL_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(colistin_bug = case_when(org_order=="Enterobacterales" & org_family!="Morganellaceae"
                                    & org_fullname!="Serratia marcescens" & org_fullname!="Hafnia"
                                    ~ "Enterobacterales",
                                    TRUE ~ "N")) %>% 
    res_sim(colistin_bug,"Enterobacterales",org_fullname,"",COL,1,10000,"Colistin",) %>%
    select(-colistin_bug) %>% 
    res_sim(org_fullname,"Pseudomonas aeruginosa",org_fullname,"",COL,1,10000,"Colistin",) %>% 
    res_sim(org_genus,"Acinetobacter",org_fullname,"",COL,1,10000,"Colistin",)
  
  
  COL_summary <- data.frame(rbind(
    `Enterobacterales_Colistin`,
    `Pseudomonas aeruginosa_Colistin`,
    `Acinetobacter_Colistin`
  ))
  
  COLdf <- data.frame(rbind(
    `Enterobacterales_Colistin_df`,
    `Pseudomonas aeruginosa_Colistin_df`,
    `Acinetobacter_Colistin_df`
  ))
  
  COLdf$org_name <- factor(COLdf$org_name, levels = COLdf %>% 
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  COL_plot <- ggplot(COLdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = COLdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Colistin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    scale_fill_brewer(palette = "Spectral") +
    theme_classic() +
    theme(legend.position = "none")
  
  print(COL_plot)
  
  assign("COL_summary",COL_summary,envir = .GlobalEnv)
  
  df
  
}
ETP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Salmonella Typhi",org_fullname,"",ETP,1,1e4,"Ertapenem") %>% 
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",ETP,1,1e4,"Ertapenem") %>%
    mutate(IPM = case_when(ETP=="S" & is.na(IPM) ~ "S",
                           TRUE ~ IPM),
           MEM = case_when(ETP=="S" & is.na(MEM) ~"S",
                           TRUE ~ MEM))
  
  
  ETP_summary <- data.frame(rbind(
    `Salmonella Typhi_Ertapenem`,
    `Haemophilus influenzae_Ertapenem`
  ))
  
  ETPdf <- data.frame(rbind(
    `Salmonella Typhi_Ertapenem_df`,
    `Haemophilus influenzae_Ertapenem_df`
  ))
  
  ETPdf$org_name <- factor(ETPdf$org_name, levels = ETPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  ETP_plot <- ggplot(ETPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = ETPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ertapenem"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(ETP_plot)
  
  assign("ETP_summary",ETP_summary,envir = .GlobalEnv)
  
  df
  
}
CRO_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    res_sim(org_fullname,"Moraxella catarrhalis",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    res_sim(org_fullname,"Neisseria meningitidis",org_fullname,"",CRO,1,1e4,"Ceftriaxone",) %>% 
    mutate(CTX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CTX) ~ "S",
                           TRUE ~ CTX),
           CAZ = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CAZ) ~ "S",
                           TRUE ~ CAZ),
           CZA = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CZA) ~ "S",
                           TRUE ~ CZA),
           CPD = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CRO=="S" & is.na(CPD) ~ "S",
                           TRUE ~ CPD))
  
  CRO_summary <- data.frame(rbind(
    `Haemophilus influenzae_Ceftriaxone`,
    `Moraxella catarrhalis_Ceftriaxone`,
    `Neisseria meningitidis_Ceftriaxone`
  ))
  
  CROdf <- data.frame(rbind(
    `Haemophilus influenzae_Ceftriaxone_df`,
    `Moraxella catarrhalis_Ceftriaxone_df`,
    `Neisseria meningitidis_Ceftriaxone_df`
  ))
  
  CROdf$org_name <- factor(CROdf$org_name, levels = CROdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CRO_plot <- ggplot(CROdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CROdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ceftriaxone"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CRO_plot)
  
  assign("CRO_summary",CRO_summary,envir = .GlobalEnv)
  
  df
  
}
CIP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  
  df <- df %>%
    res_sim(org_fullname,"Haemophilus influenzae",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    res_sim(org_fullname,"Moraxella catarrhalis",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    res_sim(org_fullname,"Neisseria meningitidis",org_fullname,"",CIP,1,1e4,"Ciprofloxacin",) %>% 
    mutate(LVX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(LVX) ~ "S",
                           TRUE ~ LVX),
           MFX = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(MFX) ~ "S",
                           TRUE ~ MFX),
           NOR = case_when(grepl("(Haemophilus influenzae|Moraxella catarrhalis|Neisseria meningitidis)",
                                 org_fullname) &
                             CIP=="S" & is.na(NOR) ~ "S",
                           TRUE ~ NOR))
  
  CIP_summary <- data.frame(rbind(
    `Haemophilus influenzae_Ciprofloxacin`,
    `Moraxella catarrhalis_Ciprofloxacin`,
    `Neisseria meningitidis_Ciprofloxacin`
  ))
  
  CIPdf <- data.frame(rbind(
    `Haemophilus influenzae_Ciprofloxacin_df`,
    `Moraxella catarrhalis_Ciprofloxacin_df`,
    `Neisseria meningitidis_Ciprofloxacin_df`
  ))
  
  CIPdf$org_name <- factor(CIPdf$org_name, levels = CIPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  CIP_plot <- ggplot(CIPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = CIPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ciprofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(CIP_plot)
  
  assign("CIP_summary",CIP_summary,envir = .GlobalEnv)
  
  df
  
}
SPT_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Neisseria gonorrhoeae",org_fullname,"",SPT,1,1e4,"Spectinomycin")
  
  SPT_summary <- data.frame(rbind(
    `Neisseria gonorrhoeae_Spectinomycin`
  ))
  
  SPTdf <- data.frame(rbind(
    `Neisseria gonorrhoeae_Spectinomycin_df`
  ))
  
  SPTdf$org_name <- factor(SPTdf$org_name, levels = SPTdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  SPT_plot <- ggplot(SPTdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = SPTdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ciprofloxacin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(SPT_plot)
  
  assign("SPT_summary",SPT_summary,envir = .GlobalEnv)
  
  df
}
VAN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    res_sim(org_fullname,"Staphylococcus aureus",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    mutate(CNS = case_when(org_genus=="Staphylococcus" & org_fullname!="Staphylococcus aureus"~ "CNS",
                           TRUE ~ "N")) %>% 
    res_sim(CNS,"CNS",org_fullname,"",VAN,1,1e4,"Vancomycin") %>%
    
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",VAN,1,1e4,"Vancomycin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",VAN,1,1e4,"Vancomycin") %>% 
    mutate(TEC = case_when((grepl("Staphylococcus aureus|Corynebacterium|Streptococcus pneumoniae)",org_fullname) |
                              grepl("BHS",BHS)) &
                             VAN=="S" & is.na(TEC) ~ "S",
                           TRUE ~ TEC),
           DAL = case_when(VAN=="S" & is.na(DAL) ~ "S",
                           TRUE ~ DAL),
           TLV = case_when(VAN=="S" & is.na(CAZ) ~ "S",
                           TRUE ~ TLV),
           ORI = case_when(VAN=="S" & is.na(CZA) ~ "S",
                           TRUE ~ ORI)) %>% 
    select(-CNS,-BHS)
  
  
  
  VAN_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Vancomycin`,
    `Staphylococcus aureus_Vancomycin`,
    `CNS_Vancomycin`,
    `BHS_Vancomycin`,
    `Corynebacterium_Vancomycin`
  ))
  
  VANdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Vancomycin_df`,
    `Staphylococcus aureus_Vancomycin_df`,
    `CNS_Vancomycin_df`,
    `BHS_Vancomycin_df`,
    `Corynebacterium_Vancomycin_df`
  ))
  
  VANdf$org_name <- factor(VANdf$org_name, levels = VANdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  VAN_plot <- ggplot(VANdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = VANdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Vancomycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(VAN_plot)
  
  assign("VAN_summary",VAN_summary,envir = .GlobalEnv)
  
  df
}
DAP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",DAP,1,1e4,"Daptomycin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",DAP,1,1e4,"Daptomycin") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",DAP,1,100,"Daptomycin")
  
  DAP_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Daptomycin`,
    `Staphylococcus_Daptomycin`,
    `BHS_Daptomycin`,
    `Corynebacterium_Daptomycin`,
    `Enterococcus_Daptomycin`
  ))
  
  DAPdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Daptomycin_df`,
    `Staphylococcus_Daptomycin_df`,
    `BHS_Daptomycin_df`,
    `Corynebacterium_Daptomycin_df`,
    `Enterococcus_Daptomycin_df`
  ))
  
  DAPdf$org_name <- factor(DAPdf$org_name, levels = DAPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  DAP_plot <- ggplot(DAPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = DAPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Daptomycin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(DAP_plot)
  
  assign("DAP_summary",DAP_summary,envir = .GlobalEnv)
  
  df
}
LNZ_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",LNZ,1,1e4,"Linezolid") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",LNZ,1,1e4,"Linezolid") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",LNZ,1,1e4,"Linezolid")
  
  LNZ_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Linezolid`,
    `Staphylococcus_Linezolid`,
    `BHS_Linezolid`,
    `Corynebacterium_Linezolid`,
    `Enterococcus_Linezolid`
  ))
  
  LNZdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Linezolid_df`,
    `Staphylococcus_Linezolid_df`,
    `BHS_Linezolid_df`,
    `Corynebacterium_Linezolid_df`,
    `Enterococcus_Linezolid_df`
  ))
  
  LNZdf$org_name <- factor(LNZdf$org_name, levels = LNZdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  LNZ_plot <- ggplot(LNZdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = LNZdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Linezolid"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(LNZ_plot)
  
  assign("LNZ_summary",LNZ_summary,envir = .GlobalEnv)
  
  df
}

PEN_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",PEN,1,1e4,"Benzylpenicillin") %>%
    mutate(OXA = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ OXA),
           AMP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMP),
           AMX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMX),
           PIP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ PIP),
           TIC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ TIC),
           CRB = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CRB),
           SAM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ SAM),
           AMC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ AMC),
           TZP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ TZP),
           LEX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ LEX),
           CZO = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CZO),
           CEC = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CEC),
           CXM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CXM),
           FOX1 = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                            TRUE ~ FOX1),
           CTX = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CTX),
           CRO = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CRO),
           CPD = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CPD),
           FEP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ FEP),
           CPT = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ CPT),
           BPR = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ BPR),
           ETP = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ ETP),
           MEM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ MEM),
           IPM = case_when(BHS=="BHS" & PEN=="S" ~ "S",
                           TRUE ~ IPM)
    ) %>% 
    select(-BHS)
  
  PEN_summary <- data.frame(rbind(
    `BHS_Benzylpenicillin`
  ))
  
  PENdf <- data.frame(rbind(
    `BHS_Benzylpenicillin_df`
  ))
  
  PENdf$org_name <- factor(PENdf$org_name, levels = PENdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  PEN_plot <- ggplot(PENdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = PENdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Benzylpenicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(PEN_plot)
  
  assign("PEN_summary",PEN_summary,envir = .GlobalEnv)
  
  df
}
AMP_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Enterococcus faecalis",org_fullname,"",AMP,1,1000,"Ampicillin")
  
  AMP_summary <- data.frame(rbind(
    `Enterococcus faecalis_Ampicillin`
  ))
  
  AMPdf <- data.frame(rbind(
    `Enterococcus faecalis_Ampicillin_df`
  ))
  
  AMPdf$org_name <- factor(AMPdf$org_name, levels = AMPdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  AMP_plot <- ggplot(AMPdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = AMPdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Ampicillin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(AMP_plot)
  
  assign("AMP_summary",AMP_summary,envir = .GlobalEnv)
  
  df
}
TEC_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    mutate(CNS = case_when(org_genus=="Staphylococcus" & org_fullname!="Staphylococcus aureus"~ "CNS",
                           TRUE ~ "N")) %>% 
    res_sim(CNS,"CNS",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>%
    
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",TEC,1,1e4,"Teicoplanin") %>% 
    select(-CNS,-BHS)
  
  
  
  TEC_summary <- data.frame(rbind(
    `CNS_Teicoplanin`,
    `BHS_Teicoplanin`,
    `Corynebacterium_Teicoplanin`
  ))
  
  TECdf <- data.frame(rbind(
    `CNS_Teicoplanin_df`,
    `BHS_Teicoplanin_df`,
    `Corynebacterium_Teicoplanin_df`
  ))
  
  TECdf$org_name <- factor(TECdf$org_name, levels = TECdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TEC_plot <- ggplot(TECdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TECdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Teicoplanin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TEC_plot)
  
  assign("TEC_summary",TEC_summary,envir = .GlobalEnv)
  
  df
}
RIF_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",RIF,1,1e4,"Rifampicin")
  
  RIF_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Rifampicin`
  ))
  
  RIFdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Rifampicin_df`
  ))
  
  RIFdf$org_name <- factor(RIFdf$org_name, levels = RIFdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  RIF_plot <- ggplot(RIFdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = RIFdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for RIFtamicin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(leRIFd.position = "none") 
  
  print(RIF_plot)
  
  assign("RIF_summary",RIF_summary,envir = .GlobalEnv)
  
  df
}
TGC_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",TGC,1,1e4,"Tigecycline") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",TGC,1,1e4,"Tigecycline") %>% 
    select(-BHS) %>% 
    res_sim(org_fullname,"Enterococcus",org_fullname,"",TGC,1,1000,"Tigecycline")
  
  TGC_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Tigecycline`,
    `Staphylococcus_Tigecycline`,
    `BHS_Tigecycline`,
    `Corynebacterium_Tigecycline`,
    `Enterococcus_Tigecycline`
  ))
  
  TGCdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Tigecycline_df`,
    `Staphylococcus_Tigecycline_df`,
    `BHS_Tigecycline_df`,
    `Corynebacterium_Tigecycline_df`,
    `Enterococcus_Tigecycline_df`
  ))
  
  TGCdf$org_name <- factor(TGCdf$org_name, levels = TGCdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  TGC_plot <- ggplot(TGCdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = TGCdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Tigecycline"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(TGC_plot)
  
  assign("TGC_summary",TGC_summary,envir = .GlobalEnv)
  
  df
}
QDA_simul <- function(df) {
  
  par(mfrow = c(1,1))
  
  df <- df %>%
    res_sim(org_fullname,"Streptococcus pneumoniae",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    res_sim(org_fullname,"Staphylococcus",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    mutate(BHS = case_when(grepl("Streptococcus Group",org_fullname) ~ "BHS",
                           TRUE ~ "N")) %>%
    res_sim(BHS,"BHS",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>%
    res_sim(org_fullname,"Corynebacterium",org_fullname,"",QDA,1,1e4,"Quinupristin-dalfopristin") %>% 
    select(-BHS)
  
  QDA_summary <- data.frame(rbind(
    `Streptococcus pneumoniae_Quinupristin-dalfopristin`,
    `Staphylococcus_Quinupristin-dalfopristin`,
    `BHS_Quinupristin-dalfopristin`,
    `Corynebacterium_Quinupristin-dalfopristin`
  ))
  
  QDAdf <- data.frame(rbind(
    `Streptococcus pneumoniae_Quinupristin-dalfopristin_df`,
    `Staphylococcus_Quinupristin-dalfopristin_df`,
    `BHS_Quinupristin-dalfopristin_df`,
    `Corynebacterium_Quinupristin-dalfopristin_df`
  ))
  
  QDAdf$org_name <- factor(QDAdf$org_name, levels = QDAdf %>%
                             distinct(org_name) %>% unlist()) %>% 
    fct_rev()
  
  QDA_plot <- ggplot(QDAdf, aes(x=org_name,y=samples,fill=org_name))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = QDAdf %>% dplyr::filter(outlier), position = 'jitter',alpha=0.05) +
    coord_flip() +
    ggtitle(glue("Posterior Resistance Probability for Quinupristin-dalfopristin"))+
    xlab("Organism Species/Group") +
    ylab("Probability") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.position = "none") 
  
  print(QDA_plot)
  
  assign("QDA_summary",QDA_summary,envir = .GlobalEnv)
  
  df
}

#I to R reassignment for sensitivity analysis

sensitivity_func <- function(df) {
  
  df2 <- df %>% select(PEN:MTR)
  df2[df2=="I"] <- "R"
  df[,17:81] <- df2
  
  df
  
}


################ Y Assign (R)

y_r_assign <- function(df,Y_var,abx,grouping_var,group) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>% mutate({{Y_var}} := case_when(!!abx == 'R' & 
                                         grepl(group,!!grouping_var,ignore.case=T) ~ TRUE,
                                       TRUE ~ FALSE))
  
}

################ Y Assign (S)

y_s_assign <- function(df,Y_var,abx,grouping_var,group) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>% mutate({{Y_var}} := case_when(!!abx == 'S' & 
                                         grepl(group,!!grouping_var,ignore.case=T) ~ 1,
                                       TRUE ~ 0))
  
}


############### Age

age_assign <- function(df,B_var,age_df,age_cutoff) {
  
  age_df %>%
    select('subject_id', 'anchor_age') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(anchor_age >= age_cutoff ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(anchor_age=NULL)
  
  
  
}

############## Gender

gender_assign <- function(df,B_var,gender_df) {
  
  gender_df %>%
    select('subject_id', 'gender') %>%
    right_join(df) %>%
    mutate({{B_var}} := case_when(gender=="M" ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(gender=NULL)
  
}

################ PREVIOUS ORGANISMS GROWN IN URINE

prev_org_assign <- function(df, B_var, org,no_days,no_events) {
  
  df %>%
    mutate(bug = case_when(
      grepl(org,org_fullname) ~ "Yes", 
      TRUE ~ "No"
    )) %>% 
    MIMER::check_previous_events(cols="bug", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_bug==TRUE ~ 2,
                                  TRUE ~ 1)) %>% 
    mutate(bug = NULL,pr_bug=NULL)
  
  
}




################ PREVIOUS EVENT

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
    mutate(event = NULL, pr_event=NULL)
  
}

################## PREVIOUS EVENT TYPE

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
    mutate(event = NULL, pr_event=NULL)
  
  
}

################ PREVIOUS ANTIMICROBIAL RESISTANCE (URINE)

prev_urine_r_assign <- function(df, B_var, abx , grouping_var, group, no_days, no_events) {
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  
  df %>%
    mutate(urine_resistance = case_when(!!abx == 'R' & 
                                          grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                        TRUE ~ "No")) %>%
    MIMER::check_previous_events(cols = "urine_resistance", sort_by_col = 'charttime',
                                 patient_id_col = 'subject_id', event_indi_value = 'Yes',
                                 new_col_prefix = "pr_R_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate({{B_var}} := case_when(pr_R_urine_resistance == TRUE ~ 2,
                                  TRUE ~ 1)) %>% 
    select(-urine_resistance, -pr_R_urine_resistance)
  
  
}

################ PREVIOUS ANTIMICROBIAL RESISTANCE (ALL)

prev_r_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_resistance = case_when(!!abx=='R' 
                                      & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                      TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_resistance", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_R_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_R_any_resistance==TRUE ~ TRUE,
                                TRUE ~ FALSE)) %>% 
    mutate(any_resistance=NULL,pr_R_any_resistance=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL SUSCEPTIBILITY (ALL)

prev_s_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='S' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_S_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_S_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_S_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL INTERMEDIATE(ALL)

prev_i_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='I' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_I_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_I_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_I_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL NOT TESTED (ALL)

prev_s_assign <- function(df, B_var, micro_df, abx , grouping_var, group, no_days, no_events) {
  
  ur_df <- df %>% mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S'))
  
  abx <- enquo(abx)
  grouping_var <- enquo(grouping_var)
  B_var <- enquo(B_var)
  
  micro_df %>%
    mutate(charttime=as.POSIXct(charttime,format='%Y-%m-%d %H:%M:%S')) %>% 
    bind_rows(ur_df) %>%  
    mutate(any_susceptibility = case_when(!!abx=='NT' 
                                          & grepl(group,!!grouping_var,ignore.case=T) ~ 'Yes',
                                          TRUE~"No")) %>%
    MIMER::check_previous_events(cols="any_susceptibility", sort_by_col='charttime',
                                 patient_id_col='subject_id', event_indi_value='Yes',
                                 new_col_prefix="pr_NT_",
                                 time_period_in_days = no_days, minimum_prev_events = no_events,
                                 default_na_date = '9999-12-31 00:00:00') %>% 
    mutate(!!B_var := case_when(pr_NT_any_susceptibility==TRUE ~ 2,
                                TRUE ~ 1)) %>% 
    mutate(any_susceptibility=NULL,pr_NT_any_susceptibility=NULL) %>% 
    filter(grepl('URINE', spec_type_desc))
  
}

################ PREVIOUS ANTIMICROBIAL TREATMENT

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

################ UNCERTAINTY PRIORITISER

R_unc_prioritiser = function(df,spec_id,panel_size) {
  
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(abs(0.5-R)) %>% select(Antimicrobial,R) %>% 
    mutate(R = round(R*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`% prob R` = "R")
  
}


################ AWARE CLASS PRIORITISER

aware_prioritiser = function(df,spec_id,panel_size,acs_weight=1,war_weight=1) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="SAM" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="TZP" ~ R * war_weight,
      as.ab(Antimicrobial)=="CZO" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="CRO" ~ R * war_weight,
      as.ab(Antimicrobial)=="CAZ" ~ R * war_weight,
      as.ab(Antimicrobial)=="FEP" ~ R * war_weight,
      as.ab(Antimicrobial)=="MEM" ~ R * war_weight,
      as.ab(Antimicrobial)=="CIP" ~ R * war_weight,
      as.ab(Antimicrobial)=="GEN" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="SXT" ~ (S+I) * acs_weight,
      as.ab(Antimicrobial)=="NIT" ~ (S+I) * acs_weight
    )) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}


################ AWARE mk2

aware_mk2 = function(df,spec_id,panel_size,acs_weight=1,war_cutoff=0.25) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="SAM" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="TZP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CZO" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="CRO" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CAZ" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="FEP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="MEM" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="CIP" & R > war_cutoff ~ 1,
      as.ab(Antimicrobial)=="GEN" ~ (S+I) * 1,
      as.ab(Antimicrobial)=="SXT" ~ S * 1,
      as.ab(Antimicrobial)=="NIT" ~ (S+I) * 1,
      TRUE ~ 0)) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}


################ AWARE mk3

aware_mk3 = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="SAM" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="TZP" ~ (S+I),
      as.ab(Antimicrobial)=="CZO" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="CRO" ~(S+I),
      as.ab(Antimicrobial)=="CAZ" ~(S+I),
      as.ab(Antimicrobial)=="FEP" ~(S+I),
      as.ab(Antimicrobial)=="MEM" ~(S+I),
      as.ab(Antimicrobial)=="CIP" ~(S+I),
      as.ab(Antimicrobial)=="GEN" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="SXT" & (S+I) > acs_cutoff ~ 1+S+I,
      as.ab(Antimicrobial)=="NIT" & (S+I) > acs_cutoff ~ 1+S+I,
      TRUE ~ 0)) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility*100,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

probs_df_overall <- read_csv("probs_df_overall.csv")

for (i in 1:100) {
  
  print(probs_df_overall %>% aware_mk3(probs_df_overall$micro_specimen_id[i],6))
  
}


################ AWARE mk3

aware_mkI = function(df,spec_id,panel_size,acs_cutoff=0.5) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    mutate(aware_utility = case_when(
      as.ab(Antimicrobial)=="AMP" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="SAM" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="TZP" ~ (S),
      as.ab(Antimicrobial)=="CZO" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="CRO" ~(S),
      as.ab(Antimicrobial)=="CAZ" ~(S),
      as.ab(Antimicrobial)=="FEP" ~(S),
      as.ab(Antimicrobial)=="MEM" ~(S),
      as.ab(Antimicrobial)=="CIP" ~(S),
      as.ab(Antimicrobial)=="GEN" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="SXT" & (S) > acs_cutoff ~ 1+S,
      as.ab(Antimicrobial)=="NIT" & (S) > acs_cutoff ~ 1+S,
      TRUE ~ S)) %>% 
    arrange(desc(aware_utility)) %>% select(Antimicrobial,aware_utility) %>% 
    mutate(aware_utility = round(aware_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Recommended tests` = "Antimicrobial",`AWaRe Utility` = "aware_utility")
  
}

####probs_df_final################ PROBABILITY PLOT

plot_probs <- function(df,chosen_test_df,spec_id) {
  
  plot_df <- chosen_test_df %>% rename(Antimicrobial="Recommended tests") %>%
    select(Antimicrobial) %>% mutate(Selected=TRUE) %>%
    right_join(df %>%
                 filter(micro_specimen_id==spec_id),by="Antimicrobial") %>% 
    mutate(Selected = case_when(is.na(Selected) ~ FALSE, TRUE~TRUE))
  
  plot_df$Antimicrobial <- factor(plot_df$Antimicrobial,
                                  levels=plot_df %>% filter(micro_specimen_id==spec_id) %>%
                                    arrange(R) %>% 
                                    select(Antimicrobial) %>% unlist())
  
  
  ggplot(plot_df %>% filter(micro_specimen_id==spec_id),
         aes(x=Antimicrobial, y=R,fill=Selected)) +
    geom_col() +
    coord_flip() +
    ylab("Resistance probability")+
    xlab("Antimicrobial agent")+
    ggtitle(glue("Resistance probability for specimen {spec_id}")) +
    geom_hline(yintercept=0.5,linetype="dashed",color="grey")+
    ylim(c(0,1))
  
}

#########################################################

#Read-in and clean up

setwd("/Users/alexhoward/Documents/Projects/UDAST_code/")
results <- read_csv("ADAPT-AST Factors influencing Antimicrobial Prescribing for Urinary Tract Infection.csv")
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
                             id = seq_rep1)




##Rank logit model

mlogit_data <- mlogit.data(scores, choice = "Rank", shape = "long", 
                           chid.var = "id", alt.var = "Antibiotic", 
                           ranked = TRUE)


formula <- Rank ~ CDI_highrisk + Toxicity_highrisk + UTI_specific + 
  Oral_option + IV_option + High_cost + Access + Reserve | 0

ranked_logit_model <- mlogit(formula, data = mlogit_data, 
                               rpar = c(CDI_highrisk = "n", Toxicity_highrisk = "n",
                                        UTI_specific = "n", Oral_option1 = "n",
                                        IV_option1 = "n", High_cost1 = "n",
                                        Access = "n", Reserve = "n"))


rank_coefs <- coef(ranked_logit_model) %>% data.frame()
rank_sds <- rank_coefs[grepl("sd.",rownames(rank_coefs)),] %>% abs()
rank_coefs$sd <- rank_sds
rank_coefs <- rank_coefs %>% slice(1:8) %>% rename(mean = ".") %>%
  mutate(mean=as.numeric(mean), sd=as.numeric(sd))
rownames(rank_coefs) <- c("High CDI risk","High toxicity risk","UTI-specific",
"Oral option","IV option","High cost", "Access category","Reserve category")
rank_coefs$Property <- rownames(rank_coefs)

rank_coefs$standardised_mean <- rank_coefs$mean / max(abs(rank_coefs$mean) + rank_coefs$sd/2)
rank_coefs$standardised_sd <- rank_coefs$sd / max(abs(rank_coefs$mean) + rank_coefs$sd/2)


rank_coefs$Property <- factor(rank_coefs$Property,
                               levels=rank_coefs %>%
                                 arrange(mean) %>% select(Property) %>%
                                           unlist())

rank_coefs <- rank_coefs %>% mutate(colour = case_when(
  mean > 0 ~ "B", TRUE ~ "A"
))



ggplot(rank_coefs,aes(x=Property,y=standardised_mean,fill=colour)) +
  geom_col() +
  geom_errorbar(aes(y=standardised_mean,ymin=standardised_mean-(standardised_sd/2),ymax=standardised_mean+(standardised_sd/2)),width=0.1) +
  coord_flip() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylim(min(rank_coefs$standardised_mean-(rank_coefs$standardised_sd/2)),-min(rank_coefs$standardised_mean-(rank_coefs$standardised_sd/2))) +
  ylab("Effect on drug preference") +
  ggtitle("The effect of different antimicrobial drug properties\non clinician prescribing preference in UTI scenario")


#Weighting base factors of respective drugs from data

ab_props <- read_csv("Ab_props.csv")
ab_props <- ab_props %>% mutate(Antimicrobial = ab_name(Antimicrobial),
                                Antimicrobial = str_replace(
                                  Antimicrobial,"/","-"
                                ))

ur_util <- read_csv("urines_assess.csv")
micro <- read_csv("micro_clean2.csv")
pos_urines <- read_csv("pos_urines.csv")
mic_ref <- micro %>% anti_join(ur_util,by="subject_id")


#C difficile base weight
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi <- cdi_ref %>% group_by(subject_id) %>% arrange(chartdate) %>% summarise_all(last) 

drugs <- read_csv("drugs_clean.csv")
abx <- drugs %>% filter(grepl(
  "(Ampicillin|Amoxicillin|Piperacillin/tazobactam|Cefazolin|Ceftriaxone|^Ceftazidime$|Cefepime|Meropenem|Ciprofloxacin|Gentamicin|Trimethoprim/sulfamethoxazole|Nitrofurantoin)",
  abx_name
)) %>% filter(!grepl("avibactam",abx_name)) %>% filter(!grepl("clavulan",abx_name)) %>% 
  mutate(abx_name = case_when(
    grepl("Amox",abx_name) ~ "Ampicillin",
          TRUE ~ abx_name)
  ) %>% 
  mutate(abx_name = case_when(
    grepl("^Ampicillin$",abx_name) ~ "Amp/Amoxicillin",
    TRUE ~ abx_name
  ))

abx <- abx %>% mutate(charttime = stoptime)
cdi_ref <- cdi_ref %>% mutate(admittime = charttime-(60*60*24*28))

#Attach CDI in following 28d labels to abx df
abx <- abx %>% 
  prev_event_type_assign(CDI,cdi_ref,org_fullname,"Clostridioides difficile",
                         28,1) %>% ungroup()  %>% 
  filter(!is.na(abx_name))
abx$CDI <- factor(abx$CDI)

#previous hospital admission
hadm <- read_csv("admissions.csv")
abx <- abx %>% 
  prev_event_assign(pHADM,hadm,hadm_id,28,1) %>% ungroup() %>% 
  filter(!is.na(abx_name))

#previous CDI
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi_ref <- cdi_ref %>% mutate(admittime=charttime)
abx <- abx %>% 
  prev_event_assign(pCDI,cdi_ref,org_fullname,1e4,1) %>% ungroup() %>% 
  filter(!is.na(abx_name))

#age > 65
pats <- read_csv("patients.csv")
pats <- pats %>% mutate(age65 = case_when(
  anchor_age >=65 ~ TRUE, TRUE~FALSE
))
patskey <- pats %>% select(subject_id,age65)
abx <- abx %>% left_join(patskey,by="subject_id")

#dummy variables for antimicrobials
recipethis <- recipe(~abx_name,data=abx)
dummies <- recipethis %>% step_dummy(abx_name) %>% prep(training = abx)
dummy_data <- bake(dummies,new_data = NULL)
abx <- abx %>% cbind(dummy_data) %>% tibble()
abx <- abx %>% mutate(abx_name_Amp.Amoxicillin = 
                            case_when(abx_name=="Amp/Amoxicillin" ~
                            1, TRUE ~ 0))



#CDI model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

log_reg_fit <- log_reg_spec %>%
  fit(CDI ~ pCDI+pHADM+age65+abx_name_Amp.Amoxicillin+
        abx_name_Ampicillin.sulbactam+abx_name_Cefazolin+
        abx_name_Cefepime+abx_name_Ceftazidime+
        abx_name_Ceftriaxone+abx_name_Ciprofloxacin+
        abx_name_Gentamicin+abx_name_Meropenem+
        abx_name_Nitrofurantoin+abx_name_Piperacillin.tazobactam+
        abx_name_Trimethoprim.sulfamethoxazole,
        data = abx)

log_reg_fit




















##Attach base utilities to probability dataframe

util_probs_df <- read_csv("probs_df_overall.csv")


for (i in 1:nrow(ab_props)) {
  
  ab_props[i,2:ncol(ab_props)] <- ab_props[i,2:ncol(ab_props)] *
    rank_coefs$standardised_mean
  
}

util_probs_df <- util_probs_df %>% left_join(ab_props,by="Antimicrobial")






############Identifying weighting factors




###REF: Antibiotics and hospital-acquired Clostridium difficile infection: update of systematic review and meta-analysis
#REF: Comparison of Different Antibiotics and the Risk for Community-Associated Clostridioides difficile Infection: A Case–Control Study

#Previous C difficile
micaborgs <- micro %>% filter(!is.na(org_name))
micabnas <- micro %>% filter(is.na(org_name))
micaborgab <- micaborgs %>% select(PEN:MTR)
micaborgab[is.na(micaborgab)] <- "NT"
micaborgs[,17:81] <- micaborgab
micro2 <- tibble(rbind(micaborgs,micabnas))
micro2 <- micro2 %>% rename(admittime = "charttime")

ur_util <- ur_util %>% 
  prev_event_type_assign(pCDI,micro2,org_fullname,"Clostridioides difficile",
                         1825,1) %>% ungroup()

#age > 65
pats <- read_csv("patients.csv")
pats <- pats %>% mutate(age65 = case_when(
  anchor_age >=65 ~ TRUE, TRUE~FALSE
))
patskey <- pats %>% select(subject_id,age65)
ur_util <- ur_util %>% left_join(patskey,by="subject_id")

#age > 80
pats <- read_csv("patients.csv")
pats <- pats %>% mutate(age80 = case_when(
  anchor_age >=80 ~ TRUE, TRUE~FALSE
))
patskey <- pats %>% select(subject_id,age80)
ur_util <- ur_util %>% left_join(patskey,by="subject_id")

#Any abx in the last 7d
ur_util <- ur_util %>% mutate(
  abx_7d = case_when(
    d7AMPrx|d7AMXrx|d7AMCrx|d7SAMrx|d7TZPrx|d7CZOrx|d7CROrx|d7CAZrx|
      d7FEPrx|d7MEMrx|d7ETPrx|d7ATMrx|d7CIPrx|d7GENrx|d7TOBrx|d7AMKrx|
      d7RIFrx|d7SXTrx|d7NITrx|d7ERYrx|d7CLRrx|d7AZMrx|d7CLIrx|d7VANrx|
      d7MTRrx|d7LNZrx|d7DAPrx|d7DOXrx ~ TRUE, TRUE~FALSE
  )
)

#diabetes
diagnoses <- read_csv("diagnoses_icd.csv")
drgcodes <- read_csv("drgcodes.csv")
diabkey <- drgcodes %>% filter(grepl("DIAB",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(diabkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pDIAB,hadm,description,
                         "DIAB",
                         1e4,1) %>% ungroup()

#heart failure
hfkey <- drgcodes %>% filter(grepl("HEART FAILURE",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(hfkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pHF,hadm,description,
                         "HEART FAILURE",
                         1e4,1) %>% ungroup()

#CKD
ckdkey <- drgcodes %>% filter(grepl("CHRONIC KIDNEY",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(ckdkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pCKD,hadm,description,
                         "CHRONIC KIDNEY",
                         1e4,1) %>% ungroup()

#LIVER DISEASE
liverkey <- drgcodes %>% filter(grepl("LIVER",description)&
                                !grepl("DELIVERY",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(liverkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pLD,hadm,description,
                         "LIVER",
                         1e4,1) %>% ungroup()

#CARDIOVASCULAR DISEASE
cardikey <- drgcodes %>% filter(grepl("HEART",description)|grepl("CARDI",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(cardikey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pCVD,hadm,description,
                         "(HEART|CARDI)",
                         1e4,1) %>% ungroup()

#HEART FAILURE
HFkey <- drgcodes %>% filter(grepl("HEART FAILURE",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(HFkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pHF,hadm,description,
                         "HEART FAILURE",
                         1e4,1) %>% ungroup()

#CANCER IN LAST YEAR
cancerkey <- drgcodes %>% filter(grepl("MALIG",description) & !grepl("EXCEPT MALI",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(cancerkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pCa,hadm,description,
                         "MALIG",
                         365,1) %>% ungroup()

#CIRRHOSIS
cirrkey <- drgcodes %>% filter(grepl("CIRRHO",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(cirrkey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pCirr,hadm,description,
                         "CIRRHO",
                         1e4,1) %>% ungroup()

#STROKE
strokekey <- drgcodes %>% filter(grepl("STROKE",description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(strokekey,by="hadm_id")
ur_util <- ur_util %>% 
  prev_event_type_assign(pStroke,hadm,description,
                         "STROKE",
                         1e4,1) %>% ungroup()



#Incorporating weighting factors

weight_key <- ur_util %>%
  select(micro_specimen_id,
         pCDI,
         age65,
         pHADM,
         abx_7d,
         pDIAB,
         pCKD,
         pLD,
         pCVD,
         pSurg,
         provider_id,
         pHF,
         age80,
         MALE
         ) #provider_id true if not admitted to hospital

util_probs_df <- util_probs_df %>%
  left_join(weight_key,by="micro_specimen_id")


#odds ratios (from literature)

#CDI ORs
util_probs_df <- util_probs_df %>% 
  mutate(
    pCDI_CDI = case_when(
      pCDI ~ 2.70, TRUE~1
    ),
    age65_CDI = case_when(
      age65 ~ 1.84, TRUE~1
    ),
    pHADM_CDI = case_when(
      pHADM ~ 1.98, TRUE~1
    ),
    abx_7d_CDI = case_when(
      abx_7d ~ 1.30, TRUE~1
    ),
    pSurg_CDI = case_when(
      pSurg ~ 1.78, TRUE~1
    ),
    CDI_util = CDI_highrisk*pCDI_CDI*age65_CDI*pHADM_CDI*pSurg_CDI*abx_7d_CDI
  )


#toxicity ORs
util_probs_df <- util_probs_df %>% 
  mutate(
    pDIAB_tox = case_when(
      pDIAB ~ 2.13, TRUE~1
    ),
    pCKD_tox = case_when(
      pCKD ~ 0.42, TRUE~1
    ),
    pLD_tox = case_when(
      pLD~ 1.83, TRUE~1
    ),
    pCVD_tox = case_when(
      pCVD ~ 1.47, TRUE~1
    ),
    pSurg_tox = case_when(
      pSurg ~ 1.31, TRUE~1
    ),
    Tox_util = Toxicity_highrisk*pDIAB_tox*pCKD_tox*pLD_tox*
      pCVD_tox*pSurg_tox
  )

#IV option
util_probs_df <- util_probs_df %>% 
  mutate(
    not_IP_IV = case_when(provider_id &
                            grepl("(Piperacillin-tazobactam|Cefazolin|Cefepime|Meropenem|Gentamicin)",Antimicrobial)
                          ~ 0, TRUE~1),
    IV_util = IV_option*not_IP_IV
  )

#Resistance utility

util_probs_df <- util_probs_df %>% 
  mutate(
    R_utility = 
      R + IV_option
  )


###Overall S utility score

util_probs_df <- util_probs_df %>% mutate(
  S_utility = 
    S * (CDI_util+Tox_util+UTI_specific+
           Oral_option+IV_util+High_cost+Access+Reserve)
)

#Overall AST utility
util_probs_df <- util_probs_df %>% mutate(
  AST_utility = 
    S_utility + R_utility 
)


#Formulary recommendations

util_probs_df %>% group_by(Antimicrobial) %>% 
  summarise(Mean_util=mean(S_utility)) %>% 
  arrange(desc(Mean_util))


#Individualised recommendations (S)

util_mk1 = function(df,spec_id,panel_size) {
  df %>% filter(micro_specimen_id==spec_id) %>%
    arrange(desc(S_utility)) %>% select(Antimicrobial,S_utility) %>% 
    mutate(S_utility = round(S_utility,1)) %>% slice(1:panel_size) %>% 
    rename(`Antimicrobial ranking` = "Antimicrobial",`Rx Utility` = "S_utility")
  
}


test_recs <-  data.frame(matrix(nrow=12,ncol=0))

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

for (i in 1:nrow(ur_util)) {
  
  rec <- util_probs_df %>% util_mk1(spec_id = ur_util$micro_specimen_id[i], panel_size = 12) %>% 
    select(1)
  
  test_recs <- cbind(test_recs,rec)
  
  print(glue("{round((i/nrow(ur_util)) * 100,0)}%"))
  
}

test_recs <- data.frame(t(test_recs))
test_recs <- data.frame(cbind(ur_util$micro_specimen_id,test_recs))
colnames(test_recs) <- c("micro_specimen_id","PDRx_1","PDRx_2","PDRx_3",
                         "PDRx_4","PDRx_5","PDRx_6","PDRx_7","PDRx_8",
                         "PDRx_9","PDRx_10","PDRx_11","PDRx_12")

ur_util <- ur_util %>% left_join(test_recs,by="micro_specimen_id")

ur_util <- ur_util %>% mutate(across(PDRx_1:PDRx_12,as.ab))

ur_util %>% count(PDRx_1) %>% arrange(desc(n))






