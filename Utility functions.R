library("brms")
library("tidyverse")
library("tidymodels")
library("openxlsx")
library("mlogit")
library("AMR")
library("glue")
library("reshape2")
library("car")
library("glmnet")
library("boot")
library("gtools")

options(error=NULL)


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
                             id = seq_rep1) %>% 
  mutate(across(Oral_option:High_cost,as.numeric))

mlogit_data <- mlogit.data(scores, choice = "Rank", shape = "long", 
                           chid.var = "id", alt.var = "Antibiotic", 
                           ranked = TRUE)






###ridge regression method

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

ggplot(scores,aes(x=OR_dif,y=Coefficient,fill=colour)) +
  geom_col() +
  theme(legend.position = "None") +
  geom_hline(aes(yintercept=0)) +
  ylab("Drug property") +
  xlab("Odds ratio for drug selection") +
  ggtitle("The effect of different antimicrobial drug properties\non clinician prescribing preference in UTI scenario")+
  scale_x_continuous(labels = function(x) x+1)+
  geom_vline(xintercept = 0)


#WEIGHTING FACTORS

ur_util <- read_csv("urines_assess.csv")
micro <- read_csv("micro_clean2.csv")
pos_urines <- read_csv("pos_urines.csv")
util_probs_df <- read_csv("probs_df_overall.csv")
mic_ref <- micro %>% anti_join(ur_util,by="subject_id")

##Drugs data frames for training and testing weighting models
drugs <- read_csv("drugs_clean.csv")
abx <- drugs %>% filter(grepl(
  "(Ampicillin|Piperacillin/tazobactam|Cefazolin|Ceftriaxone|^Ceftazidime$|Cefepime|Meropenem|Ciprofloxacin|Gentamicin|Trimethoprim/sulfamethoxazole|Nitrofurantoin)",
  abx_name
)) %>% filter(!grepl("avibactam",abx_name)) %>% anti_join(ur_util,by="subject_id") %>% 
  filter(grepl("(PO|NG|IV)",route))
abx <- abx %>% mutate(charttime = starttime,
                      ab_name = abx_name)
abx <- abx %>% filter(!is.na(starttime) & !is.na(stoptime))






###########C DIFFICILE COST
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi <- cdi_ref %>% group_by(subject_id) %>% arrange(chartdate) %>% summarise_all(last)
cdi_ref <- cdi_ref %>% mutate(admittime = charttime-(60*60*24*28))

#Attach CDI in following 28d labels to abx dfs

CDI_label <- function(df,filter_term) {

filter_term <- enquo(filter_term)
  
df <- df %>% 
  prev_event_type_assign(CDI,cdi_ref,org_fullname,"Clostridioides difficile",
                         28,1) %>% ungroup()  %>% 
  filter(!is.na(!!filter_term))
df$CDI <- factor(df$CDI)

df

}

abx <- abx %>% CDI_label(abx_name)
ur_util <- ur_util %>% CDI_label(AMP)

#previous hospital admission
hadm <- read_csv("admissions.csv")

hadm_label <- function(df,filter_term) {
  
  filter_term <- enquo(filter_term)
  
  df %>% 
    prev_event_assign(pHADM,hadm,hadm_id,28,1) %>% ungroup() %>% 
    filter(!is.na(!!filter_term))
  
}

abx <- abx %>% hadm_label(abx_name)
ur_util <- ur_util %>% hadm_label(AMP)

#previous CDI
cdi_ref <- micro %>% filter(org_fullname=="Clostridioides difficile")
cdi_ref <- cdi_ref %>% mutate(admittime=charttime)

pCDI_label <- function(df,filter_term) {
  
  filter_term <- enquo(filter_term)
  
  df %>% 
    prev_event_assign(pCDI,cdi_ref,org_fullname,1e4,1) %>% ungroup() %>% 
    filter(!is.na(!!filter_term))
  
}

abx <- abx %>% pCDI_label(abx_name)
ur_util <- ur_util %>% pCDI_label(AMP)

#age > 65
pats <- read_csv("patients.csv")
pats <- pats %>% mutate(age65 = case_when(
  anchor_age >=65 ~ TRUE, TRUE~FALSE
))
patskey <- pats %>% select(subject_id,age65)
abx <- abx %>% left_join(patskey,by="subject_id")
ur_util <- ur_util %>% left_join(patskey,by="subject_id")



########NEPHROTOXICITY COST
diagnoses <- read_csv("diagnoses_icd.csv")
drgcodes <- read_csv("drgcodes.csv")
d_icd_diagnoses <- read_csv("d_icd_diagnoses.csv")

labevents <- read_csv("labevents.csv")
d_labitems <- read_csv("d_labitems.csv")

#AKI
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
      valuenum >= (lag(valuenum)+0.3) ~ TRUE,
    TRUE ~ FALSE),
  AKI2 = case_when(
      valuenum >= (1.5*baseline) ~ TRUE,
    TRUE ~ FALSE),
  AKI3 = case_when(AKI1|AKI2 ~ TRUE, TRUE~FALSE)) %>% 
  ungroup()

write_csv(creats,"creatinines.csv")

creats <- creats %>% mutate(admittime = #adjust to search after rather than before
                              charttime - (60*60*24*7))

#check between 48h and 7d after prescription

AKI_label <- function(df,filter_term) {
  
  filter_term <- enquo(filter_term)
  
  df %>% 
    prev_event_type_assign(AKI,creats,AKI3,TRUE,
                           5,1) %>% ungroup() %>% 
    mutate(AKI = factor(AKI)) %>% 
    filter(!is.na(!!filter_term))
  
}

abx <- abx %>% AKI_label(abx_name)
ur_util <- ur_util %>% AKI_label(AMP)

#Check for other nephrotoxins & contrast
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
nephrotoxics_key <- drugs %>%
  filter(drug %in% nephrotoxins) %>% distinct(hadm_id) %>% 
  mutate(Nephrotoxic_agent = TRUE)
nephrotoxic_join <- function(df) {
  df %>% left_join(nephrotoxics_key) %>% mutate(
    Nephrotoxic_agent = case_when(is.na(Nephrotoxic_agent) ~ FALSE, TRUE ~ Nephrotoxic_agent)
  )
}

abx <- abx %>% nephrotoxic_join()
ur_util <- ur_util %>% nephrotoxic_join()

d_icd_procedures <- read_csv("d_icd_procedures.csv")
procedures <- read_csv("procedures_icd.csv")
contrast_key <- d_icd_procedures %>% filter(grepl("contrast",long_title))
contrast <- procedures %>% semi_join(contrast_key) %>% left_join(hadm_key) %>% 
  distinct(hadm_id) %>% mutate(Contrast = TRUE)

contrast_join <- function(df) {
  df %>% left_join(contrast) %>% mutate(
    Contrast = case_when(is.na(Contrast) ~ FALSE, TRUE ~ Contrast)
  )
}

abx <- abx %>% contrast_join()
ur_util <- ur_util %>% contrast_join()

AKI_adjusted_check <- function(df) {
  
  df %>% mutate(AKI = case_when(
    AKI==TRUE & Nephrotoxic_agent==FALSE & Contrast==FALSE ~ TRUE,
    TRUE ~ FALSE
  ))
  
}

abx <- abx %>% AKI_adjusted_check()
ur_util <- ur_util %>% AKI_adjusted_check()

#Previous AKI
creats <- creats %>% mutate(admittime=charttime)

prAKI_label <- function(df,filter_term) {
  
  filter_term <- enquo(filter_term)
  
  df %>% 
    prev_event_assign(prAKI,creats,AKI,1e4,1) %>% ungroup() %>% 
    filter(!is.na(!!filter_term))
  
}
creats <- creats %>% mutate(AKI = AKI3)
abx <- abx %>% prAKI_label(abx_name)
ur_util <- ur_util %>% prAKI_label(AMP)

#CKD
drgcodes <- read_csv("drgcodes.csv")
ckdkey <- drgcodes %>% filter(grepl(
  "(CHRONIC KIDNEY|DIAB|LIVER|CARD|HEART|MALIG|STROKE|SEPSIS)",
  description)) %>% 
  select(hadm_id,description)
hadm <- read_csv("admissions.csv")
hadm <- hadm %>% left_join(ckdkey,by="hadm_id")
diag_label <- function(df,filter_term,label,timeframe,newvar) {
  
  filter_term <- enquo(filter_term)
  newvar <- enquo(newvar)
  
  df <- df %>% 
    prev_event_type_assign(!!newvar,hadm,description,
                           label,
                           timeframe,1) %>% ungroup() %>% 
    filter(!is.na(!!filter_term))
  
}

abx <- abx %>% diag_label(abx_name,"CHRONIC KIDNEY",1e4,pCKD)
ur_util <- ur_util %>% diag_label(AMP,"CHRONIC KIDNEY",1e4,pCKD)

#diabetes
abx <- abx %>% diag_label(abx_name,"DIAB",1e4,pDIAB)
ur_util <- ur_util %>% diag_label(AMP,"DIAB",1e4,pDIAB)

#LIVER DISEASE
abx <- abx %>% diag_label(abx_name,"LIVER",1e4,pLIVER)
ur_util <- ur_util %>% diag_label(AMP,"LIVER",1e4,pLIVER)

#CARDIOVASCULAR DISEASE
abx <- abx %>% diag_label(abx_name,"(HEART|CARDI)",1e4,pCARD)
ur_util <- ur_util %>% diag_label(AMP,"(HEART|CARDI))",1e4,pCARD)


#CANCER IN LAST YEAR
abx <- abx %>% diag_label(abx_name,"MALIG",365,pCA)
ur_util <- ur_util %>% diag_label(AMP,"MALIG",365,pCA)

#STROKE
abx <- abx %>% diag_label(abx_name,"STROKE",1e4,pCVA)
ur_util <- ur_util %>% diag_label(AMP,"STROKE",1e4,pCVA)

#MALE
patskey2 <- pats %>% mutate(MALE = case_when(gender=="M" ~ TRUE, TRUE~FALSE)) %>% 
                            select(subject_id,MALE)
abx <- abx %>% left_join(patskey2,by="subject_id")

#CURRENT SERVICE
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

#RECENT ICU ADMISSION (last 28 days)
poe <- read_csv("poe_clean.csv")
icu <- poe %>% filter(field_value=="ICU") %>% mutate(
  field_value="ICU") %>% rename(admittime="ordertime")
abx <- abx %>% 
  prev_event_type_assign(pICU,icu,field_value,"ICU",28,1) %>%
  ungroup()

#RECENT SEPSIS
abx <- abx %>% diag_label(abx_name,"SEPSIS",365,pSEPSIS)
ur_util <- ur_util %>% diag_label(AMP,"SEPSIS",365,pSEPSIS)

#LENGTH OF ANTIBIOTIC COURSE
abx <- abx %>% mutate(course_length = as.numeric(stoptime-starttime)/24/60/60,
                      course_length=case_when(course_length<0 ~ 0,TRUE~course_length),
                      course_length=case_when(course_length>=365 ~ course_length-365,
                                              TRUE~course_length),
                      course_length = course_length/max(course_length)) %>% 
  group_by(abx_name) %>% mutate(median_course = median(course_length)) %>% ungroup()

med_course_key <- abx %>% select(abx_name,median_course) %>% 
  rename(Antimicrobial="abx_name") %>% distinct() %>% 
  mutate(Antimicrobial = str_replace(Antimicrobial,"/","-"))

#dummy variables for antimicrobials

recipethis <- recipe(~abx_name,data=abx)
dummies <- recipethis %>% step_dummy(abx_name) %>% prep(training = abx)
dummy_data <- bake(dummies,new_data = NULL)
abx <- abx %>% cbind(dummy_data) %>% tibble()
abx <- abx %>% mutate(abx_name_Ampicillin = 
                            case_when(abx_name=="Ampicillin" ~
                            1, TRUE ~ 0))


########MARROW SUPPRESSION COST

new_lowvalue <- function(df,new_colname) {
  
  new_colname <- enquo(new_colname)
  
  df %>% group_by(subject_id) %>% mutate(
    !!new_colname := case_when(
      (valuenum < as.numeric(ref_range_lower)) &
      !(lag(valuenum) < lag(as.numeric(ref_range_lower))) ~ TRUE,
      TRUE~FALSE)
  ) %>% ungroup()
  
}

new_highvalue <- function(df,new_colname) {
  
  new_colname <- enquo(new_colname)
  
  df %>% group_by(subject_id) %>% mutate(
    !!new_colname := case_when(
      (valuenum > as.numeric(ref_range_upper)) &
        !(lag(valuenum) > lag(as.numeric(ref_range_upper))) ~ TRUE,
      TRUE~FALSE)
  ) %>% ungroup()
  
  abnormal_label <- function(df,df2,new_column,search_term,filter_term) {
    
    filter_term <- enquo(filter_term)
    new_column <- enquo(new_column)
    search_term <- enquo(search_term)
    
    df2 <- df2 %>% mutate(admittime = #adjust to search after rather than before
                                  charttime - (60*60*24*7))
    
    df %>% 
      prev_event_type_assign(!!new_column,df2,!!search_term,TRUE,
                             5,1) %>% ungroup() %>% 
      mutate(!!new_column := factor(!!new_column)) %>% 
      filter(!is.na(!!filter_term))
    
  }
  
}

#WBCs
print(d_labitems %>% filter(grepl("white",label,ignore.case=T)),n=25)
wbcs <- labevents %>% filter(itemid==51301) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

wbcs <- wbcs %>% new_lowvalue(new_leukopenia)
abx <- abx %>% abnormal_label(wbcs,leukopenia,new_leukopenia,abx_name)
ur_util <- ur_util %>% abnormal_label(wbcs,leukopenia,new_leukopenia,AMP)

#Hb
print(d_labitems %>% filter(grepl("Hemoglobin",label,ignore.case=T)),n=25)
hbs <- labevents %>% filter(itemid==51222) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

hbs <- hbs %>% new_lowvalue(new_anaemia)
abx <- abx %>% abnormal_label(hbs,anaemia,new_anaemia,abx_name)
ur_util <- ur_util %>% abnormal_label(hbs,anaemia,new_anaemia,AMP)

#Plts
print(d_labitems %>% filter(grepl("Platelet",label,ignore.case=T)),n=25)
plts <- labevents %>% filter(itemid==51265) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

plts <- plts %>% new_lowvalue(new_thrombocytopenia)
abx <- abx %>% abnormal_label(plts,thrombocytopenia,new_thrombocytopenia,abx_name)
ur_util <- ur_util %>% abnormal_label(plts,thrombocytopenia,new_thrombocytopenia,AMP)

#Check for previous bleeding and other cytotoxic agents
bleed_key <- d_icd_diagnoses %>% filter(grepl("bleed",long_title,ignore.case=T) &
                             !grepl("without",long_title,ignore.case=T))
bleeding <- diagnoses %>% semi_join(bleed_key) %>% distinct(hadm_id,.keep_all = T) %>% 
  mutate(Bleeding_diagnosis = TRUE) %>% select(hadm_id,Bleeding_diagnosis)


bleed_join <- function(df) {
  df %>% left_join(bleeding) %>% mutate(
    Bleeding_diagnosis = case_when(is.na(Bleeding_diagnosis) ~ FALSE, TRUE ~ Bleeding_diagnosis)
  )
}

abx <- abx %>% bleed_join()
ur_util <- ur_util %>% bleed_join()

cytotoxins <- c("Allopurinol","Aprepitant",
"Azathioprine",
"Carmustine",
"Cisplatin",
"Cyclophosphamide",
"Dacarbazine",
"Daunorubicin",
"Dexamethasone",
"Doxorubicin hydrochloride",
"Epirubicin hydrochloride",
"Estramustine phosphate",
"Etoposide",
"Febuxostat",
"Fluorouracil",
"Idarubicin hydrochloride",
"Ifosfamide",
"Lomustine",
"Lorazepam",
"Melphalan",
"Mercaptopurine",
"Methotrexate",
"Metoclopramide hydrochloride",
"Mitomycin",
"Mitoxantrone",
"Pixantrone",
"Rasburicase",
"Vinblastine sulfate")

cytotoxics_key <- drugs %>%
  filter(drug %in% cytotoxins) %>% distinct(hadm_id) %>% 
  mutate(Cytotoxic_agent = TRUE)

cytotoxic_join <- function(df) {
  df %>% left_join(cytotoxics_key) %>% mutate(
    Cytotoxic_agent = case_when(is.na(Cytotoxic_agent) ~ FALSE, TRUE ~ Cytotoxic_agent)
  )
}

abx <- abx %>% cytotoxic_join()
ur_util <- ur_util %>% cytotoxic_join()

#Marrow suppression overall
marrow_check <- function(df) {
  
  df %>% mutate(marrow_suppress = case_when(
    (leukopenia==TRUE|anaemia==TRUE|thrombocytopenia==TRUE)
    & Bleeding_diagnosis==FALSE & Cytotoxic_agent==FALSE~TRUE, TRUE~FALSE
  ))
  
}

abx <- abx %>% marrow_check()
ur_util <- ur_util %>% marrow_check()


########LFT DERANGEMENT COST
#ALP
print(d_labitems %>% filter(grepl("Alkaline",label,ignore.case=T)),n=25)
alps <- labevents %>% filter(itemid==50863) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

alps <- alps %>% new_highvalue(new_high_alp)
abx <- abx %>% abnormal_label(alps,high_alp,new_high_alp,abx_name)
ur_util <- ur_util %>% abnormal_label(alps,high_alp,new_high_alp,AMP)

#ALT
print(d_labitems %>% filter(grepl("transferase",label,ignore.case=T)),n=25)
alts <- labevents %>% filter(itemid==50861) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

alts <- alts %>% new_highvalue(new_high_alt)
abx <- abx %>% abnormal_label(alts,high_alt,new_high_alt,abx_name)
ur_util <- ur_util %>% abnormal_label(alts,high_alt,new_high_alt,AMP)

#AST
print(d_labitems %>% filter(grepl("transferase",label,ignore.case=T)),n=25)
asts <- labevents %>% filter(itemid==50878) %>% group_by(subject_id) %>% 
  distinct(charttime,.keep_all = T) %>% ungroup()

asts <- asts %>% new_highvalue(new_high_ast)
abx <- abx %>% abnormal_label(asts,high_ast,new_high_ast,abx_name)
ur_util <- ur_util %>% abnormal_label(asts,high_ast,new_high_ast,AMP)

#Check for previous liver disease and recent biliary instrumentation
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

bil_proc_join <- function(df) {
  df %>% left_join(biliary) %>% mutate(
    Biliary_procedure = case_when(is.na(Biliary_procedure) ~ FALSE, TRUE ~ Biliary_procedure)
  )
}

abx <- abx %>% bil_proc_join()
ur_util <- ur_util %>% bil_proc_join()


#Overall LFT derangement
lft_check <- function(df) {
  
  df %>% mutate(deranged_lfts = case_when(
    (high_alp==TRUE|high_alt==TRUE|high_ast==TRUE) &
      pLIVER==FALSE&Biliary_procedure==FALSE~TRUE, TRUE~FALSE
  ))
  
}
abx <- abx %>% lft_check()
ur_util <- ur_util %>% lft_check()

#Overall toxicity
toxicity_check <- function(df) {
  
  df %>% mutate(overall_tox = case_when(
    AKI==TRUE|marrow_suppress==TRUE|deranged_lfts==TRUE ~TRUE, TRUE~FALSE
  ))
  
}
abx <- abx %>% toxicity_check()
ur_util <- ur_util %>% toxicity_check()

#row ids
abx <- abx %>% mutate(row_id = seq(1,nrow(abx))) %>% 
  relocate(row_id,.before = "subject_id")

write_csv(abx,"interim_abx.csv")
write_csv(ur_util,"interim_ur_util.csv")

abx <- read_csv("interim_abx.csv")
util_probs_df <- read_csv("probs_df_overall.csv")
ur_util <- read_csv("interim_ur_util.csv")

#split abx into train_test
subjects <- abx %>% distinct(subject_id)
smp_size <- floor(0.8 * nrow(subjects))
set.seed(123)
train_ind <- sample(seq_len(nrow(subjects)), size = smp_size)
train_ids <- subjects[train_ind,]
test_ids <- subjects[-train_ind,]

abx <- abx %>% mutate(CDI = factor(CDI),
                      overall_tox = factor(overall_tox))

train_abx <- abx %>% semi_join(train_ids,by="subject_id")
test_abx <- abx %>% semi_join(test_ids,by="subject_id")







##############CDI model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

#fit model (ampicillin excluded)
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

#Evaluate on test data
cdi_test_probs <- predict(underlying_cdi, test_abx,type="response")

#Attach to probs_df
cdi_util_key <- ur_util %>% select(micro_specimen_id,pHADM,MALE,pICU,
                                   CDI:pSEPSIS) %>% 
  select(-AKI)
cdi_util_key
util_probs_df <- util_probs_df %>% 
  left_join(cdi_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

############Toxicity model
log_reg_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glm") %>%
  set_mode("classification")

#fit model (ampicillin excluded)
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

#evaluate on test data
tox_test_probs <- predict(underlying_tox, test_abx,type="response")

#Attach to probs_df
tox_util_key <- ur_util %>% select(micro_specimen_id,overall_tox)
util_probs_df <- util_probs_df %>% 
  left_join(tox_util_key,by="micro_specimen_id",
            relationship = "many-to-one")

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

####CURRENT POSITION############

cdi_util_probs <- predict(underlying_cdi, util_probs_df,type="response")
cdi_value <- scores[rownames(scores)=="CDI_highrisk",] %>% 
  select(stan_OR) %>% unlist()

tox_util_probs <- predict(underlying_tox, util_probs_df,type="response")
tox_value <- scores[rownames(scores)=="Toxicity_highrisk",] %>% 
  select(stan_OR) %>% unlist()

uti_specifics <- c("Nitrofurantoin")
uti_value <- scores[rownames(scores)=="UTI_specific",] %>% 
  select(stan_OR) %>% unlist()

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
access_value <- scores[rownames(scores)=="Access",] %>% 
  select(stan_OR) %>% unlist()

oral_abs <- c("AMP","SAM","CIP",
                "GEN","SXT","NIT") %>% ab_name() %>% 
  str_replace("/","-")
oral_value <- scores[rownames(scores)=="Oral_option",] %>% 
  select(stan_OR) %>% unlist()

iv_abs <- c("AMP","SAM","TZP","CIP","FEP","CAZ","CRO","CZO","MEM",
              "GEN","SXT","VAN") %>% ab_name() %>% 
  str_replace("/","-")
iv_value <- scores[rownames(scores)=="IV_option",] %>% 
  select(stan_OR) %>% unlist()

reserve_abs <- c()
reserve_value <- scores[rownames(scores)=="Reserve",] %>% 
  select(stan_OR) %>% unlist()

highcost_abs <- c()
cost_value <- scores[rownames(scores)=="High_cost",] %>% 
  select(stan_OR) %>% unlist()



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
         Highcost_agent = case_when(Antimicrobial %in% highcost_abs ~ 1, TRUE~0),
         value_highcost = cost_value,
         util_highcost = Highcost_agent * value_highcost,
         )





#OVERALL UTILITY SCORE

util_probs_df <- util_probs_df %>% mutate(overall_util = util_CDI + util_tox +
                                            util_uti + util_access +
                                            util_oral + util_iv +
                                            util_reserve + util_highcost,
                      S_utility = S*overall_util)






#RECOMMENDATIONS

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


test_recs <-  data.frame(matrix(nrow=13,ncol=0))

access_abs <- c("AMP","SAM","CZO",
                "GEN","SXT","NIT")

for (i in 1:nrow(ur_util)) {
  
  rec <- util_probs_df %>% util_mk1(spec_id = ur_util$micro_specimen_id[i], panel_size = 13) %>% 
    select(1)
  
  test_recs <- cbind(test_recs,rec)
  
  print(glue("{round((i/nrow(ur_util)) * 100,0)}%"))
  
}

test_recs <- data.frame(t(test_recs))
test_recs <- data.frame(cbind(ur_util$micro_specimen_id,test_recs))
colnames(test_recs) <- c("micro_specimen_id","PDRx_1","PDRx_2","PDRx_3",
                         "PDRx_4","PDRx_5","PDRx_6","PDRx_7","PDRx_8",
                         "PDRx_9","PDRx_10","PDRx_11","PDRx_12","PDRx_13")

ur_util <- ur_util %>% left_join(test_recs,by="micro_specimen_id")

ur_util <- ur_util %>% mutate(across(PDRx_1:PDRx_13,as.ab))

ur_util %>% count(PDRx_1) %>% arrange(desc(n))



