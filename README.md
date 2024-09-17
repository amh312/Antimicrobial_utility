This repository houses the code for the academic paper **"A universal antimicrobial utility function for personalised infection care"**, for the purpose of peer review and subsequent open-sourcing.

If you use this code please cite this repository.

***Instructions for use:***

The source data can be obtained from PhysioNet at https://physionet.org/content/mimiciv/2.2/ once the terms of access are met. The csv filenames used in this code match the following default filenames that can be downloaded from the *hosp* folder at the bottom of the page: *"prescriptions.csv", "d_icd_diagnoses.csv", "diagnoses_icd.csv", "d_icd_procedures.csv", "procedures_icd.csv", "labevents.csv", "d_labitems.csv", "microbiologyevents.csv", "poe_detail.csv", "poe.csv", "omr.csv", "admissions.csv", "patients.csv"*, and *"services.csv"*.

PhysioNet MIMIC-IV citations:

*Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark, R. (2023). MIMIC-IV (version 2.2). PhysioNet. https://doi.org/10.13026/6mm1-ek67.*

*Johnson, A.E.W., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data 10, 1 (2023). https://doi.org/10.1038/s41597-022-01899-x*

*Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.*

This code was written and run using *R* version 4.3.2 and *Python* version 3.12.0, on a laptop computer running macOS Sonoma version 14.5 with an Apple M1 Pro processor, 16GB random-access memory and 10 cores. The code may need to be run in chunks, depending on application memory. The typical run time of all code was approximately 3-4 hours.

Before running the code, the data and *Python* files should be saved into a secure local directory - the filepath of this directory should then be substituted into the file in place of **#FILEPATH#** in all R scripts before running the code. The required package versions are included in the *packages.csv* file within this directory. A conda environment was used to run the *Reticulate* interface package - the local environment used should be substituted for **#CONDAENV_FILEPATH#** in the **PDAST_microsimulation.R**, **PDAST_microsim_I_sensitivity.R**, **PDAST_time_sensitivity.R**, and **app.R** scripts.

***Reproducing the study***

1. To reproduce the study, all scripts above must be stored in the working directory. The following scripts must then be run:  
   a. **UF_packages&setup.R**  
   b. **UF_cleaning.R**  
   c. **UF_prediction&preprocessing.R**  
   d. **UF_utility_calculations.R**  
   e. **UF_recommendations.R**  
   f. **UF_microsimulation.R**
