rm(list = ls())

source("Scripts/Helpers/init.R")
source("Scripts/Helpers/differentialAnalysis.R")
source("Scripts/Helpers/hierClust.R")
source("Scripts/Helpers/demoHist.R")
source("Scripts/Helpers/enrichAnalysis.R")
source("Scripts/Helpers/wgcna.R")

dir <- "~/Developer/Data/proteomics_project_data"
get_biodata <- function(path) {
  read.csv(paste(dir, path, sep = "/"))
}

# same thing with stanford data? ----
stanford_barcodes <- get_biodata("ADRC/barcodes.csv")
stanford_CSF_patients <- get_biodata("ADRC/CSF_patients.csv")
stanford_CSF_prots <- get_biodata("ADRC/CSF_prots.csv")
stanford_plasma_patients <- get_biodata("ADRC/plasma_patients.csv")
stanford_plasma_prots <- get_biodata("ADRC/plasma_prots.csv")

stanford_CSF_patients$Age <- gsub('^"&"$', '', stanford_CSF_patients$Age)
stanford_CSF_patients$Age <- as.numeric(stanford_CSF_patients$Age)

stanford_plasma_prots <- stanford_plasma_prots[stanford_plasma_prots$Barcode %in% stanford_plasma_patients$Barcode, ]
stanford_CSF_prots <- data.frame(Barcode = stanford_CSF_patients$Barcode, stanford_CSF_prots)

stanford_common_prots <- intersect(names(stanford_CSF_prots), names(stanford_plasma_prots))

stanford_plasma_prots <- stanford_plasma_prots[stanford_common_prots]
stanford_CSF_prots <- stanford_CSF_prots[stanford_common_prots]

stanford_plasma_prots <- stanford_plasma_prots[stanford_plasma_prots$Barcode %in% stanford_barcodes$Plasma_Barcode, ]
stanford_plasma_patients <- stanford_plasma_patients[stanford_plasma_patients$Barcode %in% stanford_barcodes$Plasma_Barcode, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_CSF_prots$Barcode %in% stanford_barcodes$CSF_Barcode, ]
stanford_CSF_patients <- stanford_CSF_patients[stanford_CSF_patients$Barcode %in% stanford_barcodes$CSF_Barcode, ]

stanford_plasma_prots <- stanford_plasma_prots %>%
  filter(Barcode %in% stanford_barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, stanford_barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, stanford_barcodes$Plasma_Barcode)))

stanford_plasma_patients <- stanford_plasma_patients %>%
  filter(Barcode %in% stanford_barcodes$Plasma_Barcode) %>% 
  arrange(match(Barcode, stanford_barcodes$Plasma_Barcode)) %>% 
  slice(order(match(Barcode, stanford_barcodes$Plasma_Barcode)))

stanford_gendermatch_index <- stanford_CSF_patients$Gender == stanford_plasma_patients$Gender
stanford_CSF_patients <- stanford_CSF_patients[stanford_gendermatch_index, ]
stanford_plasma_patients <- stanford_plasma_patients[stanford_gendermatch_index, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_gendermatch_index, ]
stanford_plasma_prots <- stanford_plasma_prots[stanford_gendermatch_index, ]
stanford_barcodes <- stanford_barcodes[stanford_gendermatch_index, ]

stanford_patients <- stanford_plasma_patients %>% 
  mutate(Plasma_Barcode = Barcode, CSF_Barcode = stanford_CSF_patients$Barcode, 
         Plasma_drawage = Age, CSF_drawage = stanford_CSF_patients$Age, 
         Plasma_storagedays = Storage_days, 
         CSF_storagedays = stanford_CSF_patients$Storage_days,
         CSF_days = stanford_CSF_patients$Storage_days, 
         Plasma_drawdate = strptime(Date.of.draw, format = "%m/%d/%y"), 
         CSF_drawdate = strptime(stanford_CSF_patients$Date_of_CSF_sample, format = "%m/%d/%y"), 
         Plasma_drawstatus = Diagnosis_group, 
         CSF_drawstatus = stanford_CSF_patients$Diagnosis_group) %>% 
  group_by(Plasma_Barcode) %>%
  mutate(avg_drawage = (Plasma_drawage + as.numeric(CSF_drawage)) / 2, 
         drawdate_diff = difftime(Plasma_drawdate, CSF_drawdate, units = "days"), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>%
  mutate(Plasma_drawdate2 = as.numeric(difftime(today(), Plasma_drawdate)), batch_effect = 1) %>%
  dplyr::select(Plasma_Barcode, CSF_Barcode, Gender, avg_drawage, Plasma_drawdate, 
                Plasma_drawage, Plasma_drawstatus, Plasma_storagedays, CSF_drawdate, 
                CSF_drawage, CSF_drawstatus, CSF_storagedays, drawdate_diff,
                Plasma_drawdate2, avg_drawdate, batch_effect) %>%
  as.data.frame()


# actually drawdate diff is all below 90!!! we can proceed without setting a time window. 
stanford_patients <- stanford_patients %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "ADMCI", "AD", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "MCI-AD", "AD", CSF_drawstatus)) %>% 
  mutate(Plasma_drawstatus = 
           if_else(Plasma_drawstatus == "HC", "CO", Plasma_drawstatus)) %>% 
  mutate(CSF_drawstatus = 
           if_else(CSF_drawstatus == "HC", "CO", CSF_drawstatus)) %>% 
  mutate(final_status = 
           if_else(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus))

stanford_AD_index <- stanford_patients$final_status == "AD" 
stanford_CO_index <- stanford_patients$final_status == "CO"
stanford_AD_index[is.na(stanford_AD_index)] <- FALSE
stanford_CO_index[is.na(stanford_CO_index)] <- FALSE
stanford_index <- stanford_AD_index | stanford_CO_index

# clean WashU data 
knight_plasma_meta <- get_biodata("Knight/plasma/WashU_Plasma_SomaScan7K_sample_metadata.csv")
knight_plasma_prots <- get_biodata("Knight/plasma/WashU_Plasma_SomaScan7K_sample_protein_expression.csv")
knight_plasma_patients <- read.xlsx(paste(dir, "Knight/plasma/Plasma_ BasicDemo.xlsx", sep = "/"))
knight_CSF_meta <- get_biodata("Knight/CSF/WashU_CSF_SomaScan7K_sample_metadata.csv")
knight_CSF_prots <- get_biodata("Knight/CSF/WashU_CSF_SomaScan7K_sample_protein_expression.csv")
knight_CSF_patients <- read.xlsx(paste(dir, "Knight/CSF/CSF_BasicDemo.xlsx", sep = "/"))

knight_patients_meta <- get_biodata("Knight/WashU_ADNI_demographic.csv")

knight_plasma_prots[4:ncol(knight_plasma_prots)] <- log10(knight_plasma_prots[4:ncol(knight_plasma_prots)])
knight_CSF_prots[4:ncol(knight_CSF_prots)] <- log10(knight_CSF_prots[4:ncol(knight_CSF_prots)])

protMeta <- read_csv("~/Developer/Data/proteomics_project_data/Knight/ProteinMetadata_with_LOD.csv", col_types = )
protMeta$SeqId = str_replace(protMeta$SeqId, "-", ".")
protMeta <- filter(protMeta, (Organism == "Human" | Organism == "HIV-1" | Organism == "HIV-2") & Type == "Protein")
all_prots_ID <- names(knight_plasma_prots)[4:ncol(knight_plasma_prots)] 
all_prots <- str_replace(all_prots_ID, ".*?(\\d+)\\.(\\d+)\\.\\d+$", "\\1.\\2")
my_order <- match(protMeta$SeqId, all_prots)
all_prots <- protMeta$Key_2
nprots <- length(all_prots)
my_order <- my_order + 3

knight_plasma_prots <- knight_plasma_prots[c(1:3, my_order)]
colnames(knight_plasma_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)

knight_CSF_prots <- knight_CSF_prots[c(1:3, my_order)]
colnames(knight_CSF_prots) <- c("PA_DB_UID", "Barcode2d", "ExtIdentifier", all_prots)


all_prots <- intersect(names(knight_plasma_prots)[4:ncol(knight_plasma_prots)], 
                       names(knight_CSF_prots)[4:ncol(knight_CSF_prots)])
nprots <- length(all_prots)

knight_plasma_patients <- knight_plasma_patients %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed AD", 
                   "AD", final_cc_status.updated)) %>% 
  mutate(final_cc_status.updated = 
           if_else(final_cc_status.updated == "Neuropath Confirmed Control", 
                   "CO", final_cc_status.updated)) %>% 
  mutate(DateOfBirth = if_else(grepl("/", DateOfBirth), strptime(DateOfBirth, format = "%m/%d/%Y"), 
                               convertToDateTime(as.numeric(DateOfBirth)))) %>% 
  mutate(drawdate = convertToDateTime(as.numeric(drawdate)))

knight_CSF_patients <- knight_CSF_patients %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed AD", 
                   "AD", Last.status)) %>% 
  mutate(Last.status = 
           if_else(Last.status == "Neuropath Confirmed Control", 
                   "CO", Last.status)) %>% 
  mutate(DrawDate = convertToDateTime(as.numeric(DrawDate)))

knight_plasma_meta <- arrange(knight_plasma_meta, PA_DB_UID)
knight_plasma_prots <- arrange(knight_plasma_prots, PA_DB_UID)
knight_plasma_patients <- arrange(knight_plasma_patients, PA_DB_UID)
knight_CSF_meta <- arrange(knight_CSF_meta, PA_DB_UID)
knight_CSF_prots <- arrange(knight_CSF_prots, PA_DB_UID)
knight_CSF_patients <- arrange(knight_CSF_patients, PA_DB_UID)

knight_common_ID <- intersect(knight_CSF_patients$PA_DB_UID, knight_plasma_patients$PA_DB_UID) %>% 
  data.frame(PA_DB_UID = .) %>% arrange(PA_DB_UID) 
knight_common_ID <- knight_common_ID$PA_DB_UID

knight_plasma_meta <- knight_plasma_meta[knight_plasma_meta$PA_DB_UID %in% knight_common_ID, ]
knight_plasma_prots <- knight_plasma_prots[knight_plasma_prots$PA_DB_UID %in% knight_common_ID, ]
knight_plasma_patients <- knight_plasma_patients[knight_plasma_patients$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_meta <- knight_CSF_meta[knight_CSF_meta$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_prots <- knight_CSF_prots[knight_CSF_prots$PA_DB_UID %in% knight_common_ID, ]
knight_CSF_patients <- knight_CSF_patients[knight_CSF_patients$PA_DB_UID %in% knight_common_ID, ]

knight_patients <- merge(knight_CSF_patients, knight_plasma_patients, by = "PA_DB_UID") %>% 
  select(PA_DB_UID = PA_DB_UID, Sex = gender, DOB = DateOfBirth, 
         Plasma_drawdate = drawdate, Plasma_drawage = Age_at_blood_draw.updated, 
         Plasma_drawstatus = final_cc_status.updated, 
         CSF_drawdate = DrawDate, CSF_drawage = age_at_csf_draw, 
         CSF_drawstatus = Last.status) %>%
  group_by(PA_DB_UID) %>% 
  mutate(drawdate_diff = abs(difftime(Plasma_drawdate, CSF_drawdate, units = "days")), 
         avg_drawdate = as.Date(mean(c(Plasma_drawdate, CSF_drawdate)))) %>% 
  mutate(storage_days = as.numeric(difftime(today(), Plasma_drawdate, units = "days")), batch_effect = 0) %>% 
  as.data.frame() 

date_window <- 120
knight_date_index <- knight_patients$drawdate_diff <= date_window
knight_CSF_prots <- knight_CSF_prots[knight_date_index, ]
knight_plasma_prots <- knight_plasma_prots[knight_date_index, ]
knight_patients <- knight_patients[knight_date_index, ]

# get rid of non AD or COs
knight_index <- knight_patients$final_status == "AD" | knight_patients$final_status == "CO" 
stanford_index <- stanford_patients$final_status == "AD" | stanford_patients$final_status == "CO"
knight_patients <- knight_patients[knight_index, ]
stanford_patients <- stanford_patients[stanford_index, ]
knight_plasma_prots <- knight_plasma_prots[knight_index, ]
stanford_plasma_prots <- stanford_plasma_prots[stanford_index, ]
knight_CSF_prots <- knight_CSF_prots[knight_index, ]
stanford_CSF_prots <- stanford_CSF_prots[stanford_index, ]

knight_patients <- knight_patients %>% 
  mutate(final_status = ifelse(Plasma_drawdate < CSF_drawdate, CSF_drawstatus, Plasma_drawstatus), 
         avg_drawage = abs(Plasma_drawage + CSF_drawage) / 2, storage_days2 = storage_days) %>% 
  select(PA_DB_UID, Sex, avg_drawage, drawdate_diff, Plasma_drawdate, CSF_drawdate, 
         final_status, storage_days, storage_days2, batch_effect) 








