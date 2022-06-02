#### Tucson Data Cleaning
#### heidiesteiner@arizona.edu
#### Jan 30 2022



#### call packages
library(tidyverse)
library(REDCapR)

#### load data 
#### helpful resource for REDCapR 
#### https://github.com/nutterb/redcapAPI/wiki 
tucson_clinical <- redcap_read(redcap_uri="https://redcap.uahs.arizona.edu/api/", 
                               token="023F8FB7D25A42F1C4E8996484A0CC6F",
                               verbose = F,
                               raw_or_label_headers = "label")$data %>% 
  select(-`Date of Birth`, -`Initials...3`, -`Complete?...87`, -`Complete?...99`) %>% # remove PHI before it loads
  discard(~all(is.na(.) | . == 0)) # remove everything entirely NA or without variability 


#### clean data
tucson_clinical2 <- tucson_clinical %>% 
  select(id = `Study ID`, # select columns desired
         sex = Sex...4, 
         age = Age, 
         race_AI = `Race/Ethnicity (check all that apply) (choice=American Indian or Alaska Native)`,
         race_white = `Race/Ethnicity (check all that apply) (choice=White)`,
         height = `Height (inches)`, weight = `Weight (Pounds)`,
         smoke = `Do you smoke cigarettes?`,
         drink = `Do you drink alcohol?`,
         dose = `Weekly warfarin dose (mg/week)`,
         target2_3 = `Target INR (choice=2.0-3.0)`,
         target25_35 = `Target INR (choice=2.5-3.5)`,
         targetother = `Target INR (choice=other)`,
         indication_dvt = `Indication(s) for warfarin (check all that apply) (choice=Deep venous thrombosis (DVT))`,
         indication_pe = `Indication(s) for warfarin (check all that apply) (choice=Pulmonary embolism (PE))`,
         indication_tia = `Indication(s) for warfarin (check all that apply) (choice=Stroke/transient ischemic attack (TIA))`,
         indicataion_af = `Indication(s) for warfarin (check all that apply) (choice=Atrial fibrillation (A.Fib) or flultter)`,
         indication_valve = `Indication(s) for warfarin (check all that apply) (choice=Prosthetic Valve)`,
         indication_other = `Indication(s) for warfarin (check all that apply) (choice=Other)`,
         diabetes = `Past Medical History (check all that apply) (choice=Diabetes)`,
         htn = `Past Medical History (check all that apply) (choice=HTN)`,
         medlist = `Current Medication List`) %>% 
  mutate(met = if_else(grepl("metf", medlist, ignore.case = T),1,0), # create a metformin variable
         amio = if_else(grepl("amiod", medlist, ignore.case = T), 1,0), # create an amiodarone variable
         ei = if_else(grepl("carba|pheny|rif", medlist, ignore.case = T),1,0), # create an enzyme inducer varible 
         dose = parse_number(dose), # extract only digits from dose column
         IID = toupper(id), # make ids uniform 
         ethnicity_latino = 1, # create ethnicity variable
         bsa =  sqrt((height*weight)/ 3600),
         height_cm = 2.54*height,
         weight_kg = 0.453592*weight) %>% 
  mutate_at(vars(matches(c("race", "indic", "target", "smoke", "drink", "diabetes", "htn", "met", "amio", "ei", "dose", "age", "height", "weight"))), as.numeric) %>% 
  select(-id, -height, -weight, -medlist)

#### write out data
write_delim(tucson_clinical2, "results/datasets/tucson_redcap.tsv", delim = "\t")


