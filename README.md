This repository houses the code for the academic paper **"A utility function for antibiotic decision support systems: a simulation study"**, for the purpose of peer review and subsequent open-sourcing.

If you use this code please cite this repository.

***Instructions for use:***

The electronic healthcare record source data can be obtained from PhysioNet at https://physionet.org/content/mimiciv/2.2/ and https://physionet.org/content/mimic-iv-ed/2.2/ once the terms of access are met. The csv filenames used in this code match the following default filenames that can be downloaded from the site: *"prescriptions.csv", "d_icd_diagnoses.csv", "diagnoses_icd.csv", "d_icd_procedures.csv", "procedures_icd.csv", "labevents.csv", "d_labitems.csv", "microbiologyevents.csv", "poe_detail.csv", "poe.csv", "omr.csv", "admissions.csv", "patients.csv"*, *"services.csv"*, *"edstays.csv"*, *"triage.csv"*, *"pyxis.csv"*, *"medrecon.csv"*, and *"vitalsign.csv"*.

PhysioNet MIMIC-IV citations:

*Johnson, A., Bulgarelli, L., Pollard, T., Horng, S., Celi, L. A., & Mark, R. (2023). MIMIC-IV (version 2.2). PhysioNet. https://doi.org/10.13026/6mm1-ek67.*

*Johnson, A.E.W., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data 10, 1 (2023). https://doi.org/10.1038/s41597-022-01899-x*

*Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.*

This code was written and run using *R* version 4.3.2 on a laptop computer running macOS Sonoma version 14.5 with an Apple M1 Pro processor, 16GB random-access memory and 10 cores. The code may need to be run in chunks, depending on application memory. The typical run time of all code was approximately 5 days.

Before running the code, the data should be saved into a secure local directory - the filepath of this directory should then be substituted into the file in place of **#FILEPATH#** in all R scripts before running the code. The required package versions are included in the *packages.csv* file within this directory.

***Reproducing the study***

Please note that this code will exactly reproduce the clinical prediction model results of the study and descriptive data (scripts 1-8 and 12), but scripts 9 and 10 run probability-weighted randomly-generated discrete choice experiment results, so the results from scripts 9-11 will randomly vary each time the code is run.

To run the code, the *"DCE_characteristics.csv"* file must be stored in the working directory, and list that directory in place of #FILEPATH# in the script(s) marked with an asterisk below. The scripts must then be run in this order:  

   1. **UF_packages&setup.R***
   2. **UF_cleaning.R**
   3. **UF_preprocessing.R**
   4. **UF_urine_model_tuning.R**
   5. **UF_urine_model_validate.R**
   6. **UF_prescription_model_tuning.R**
   7. **UF_prescription_model_validate.R**
   8. **UF_model_testing.R**  
   9. **UF_DCE.git.R**
   10. **UF_spec_git.R**
   11. **UF_microsimulation.R**
   12. **UF_descriptive.R**

The code in its entirety typically takes approximately 72-96 hours to run. Please note that the code also includes additional analyses that are not included in the paper for the purpose of brevity, that pertain to combination treatment.

Scripts 4 and 6 (hyperparameter tuning) are the most computationally time-consuming and may need to be run in chunks, depending on RAM. To skip these tuning steps, the csv files in the **hyperparameters** folder can be downloaded to provide the hyperparameters, and all scripts except 4 and 6 run.
