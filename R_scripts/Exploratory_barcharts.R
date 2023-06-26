# Exploratory bar charts script
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load the required data #####

# Breast Cancer ###
breast_cancer_study_chars = read.xlsx("data/Study_characteristics.xlsx", sheet = 1)
breast_cancer_study_chars$Dataset = c("C1", "C2", "C3", "E1_1", "E1_2", "E1_3",
                                      "E2", "E3", "E4_1", "E4_2", "XVC1", "XVC2",
                                      "XVE", "XVC3")
breast_cancer_Pheno = read.xlsx("data/Output sets/Pheno.xlsx", sheet = 1)
breast_cancer_extval = read.xlsx("data/Output sets/Ext_val.xlsx")
breast_cancer_full_pheno = rbind(breast_cancer_Pheno, breast_cancer_extval)
breast_cancer_full_pheno = breast_cancer_full_pheno %>% 
  inner_join(breast_cancer_study_chars %>%
               dplyr::select(Dataset, Publication, Year, Samples, Patients,
                             Location, ER.status, Menopause.status,
                             `Age.(years)`, Design = `Timepoints`),
             by = "Dataset")

# Fix Mammaprint column
breast_cancer_full_pheno$new_Mamma = NA
breast_cancer_full_pheno$new_Mamma[breast_cancer_full_pheno$Mammaprint_risk == 1] = "Risk"
breast_cancer_full_pheno$new_Mamma[breast_cancer_full_pheno$Mammaprint_risk == 0] = "No risk"
breast_cancer_full_pheno = breast_cancer_full_pheno %>%
  dplyr::select(-Mammaprint_risk) %>%
  dplyr::rename(Mammaprint_risk = new_Mamma)

rm(breast_cancer_Pheno, breast_cancer_extval, breast_cancer_study_chars); gc()

# Plotting #####

# Samples and Patients
n_patients = distinct(breast_cancer_full_pheno[, c("Dataset", "Patient.ID")]) %>%
  count(Dataset)
labels = breast_cancer_full_pheno %>%
  count(Dataset) %>%
  mutate(n_patients = n_patients$n); rm(n_patients); gc()

sample_patient_barchart = ggplot(data = labels, aes(x = Dataset)) +
  geom_bar(aes(y = n, fill = "dodgerblue2"), color = "turquoise", 
           stat = "identity", position = "identity")+
  geom_text(aes(y = n, label = n), vjust = 1.4, 
            color = "white", size = 2*1.5, fontface = "bold")+
  geom_bar(aes(y = n_patients, fill = "orange2"), color = "orange",
           stat = "identity", position = "identity")+
  geom_text(aes(y = n_patients, label = n_patients), vjust = 1.4, 
            color="white", size = 2*1.5, fontface = "bold")+
  scale_fill_manual(name = "", values = c("dodgerblue2" = "dodgerblue2",
                                              "orange2" = "orange2"),
                    labels = c("Samples", "Patients"))+
  ylab("Counts")+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 2*5, hjust = 0.5),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 2*4.5),
        axis.text.x = element_text(size = 2*3.5, margin = ggplot2::margin(t = -7)),
        axis.text.y = element_text(size = 2*3.5),
        legend.position = "none",
        panel.grid.major.y = element_line(linewidth = 0.1, color = "grey80"),
        panel.grid.minor.y = element_line(linewidth = 0.1, color = "grey95"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(color = "white")) +
  labs(title = "Samples and Patients bar plot across studies")

sample_patient_barchart
ggsave(filename = "sample_patient_barchart.tiff",
       path = "data/Exploratory barcharts",
       width = 2*1920, height = 2*1920, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()
rm(labels); gc()

# Response
response_barchart = ggplot(data = breast_cancer_full_pheno, aes(x = Response,
                                                                fill = Response))+
  scale_x_discrete(labels = c("Non-responder", "Responder"))+
  geom_bar(linewidth = 0.25, width = 0.35, color = "grey20")+
  geom_text(aes(label = after_stat(count)), vjust = 1.6, 
            color = "white", size = 2*1.5, stat = "count", fontface = "bold") +
  scale_fill_manual(name = "Group", values = c(Responder = "dodgerblue4",
                                               Non_responder = "deeppink4"),
                    labels = c("Responders", "Non-responders"))+
  ylab("Number of samples")+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 2*5, hjust = 0.5),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 2*4.5),
        axis.text.x = element_text(size = 2*3.5, margin = ggplot2::margin(t = -3)),
        axis.text.y = element_text(size = 2*3.5),
        legend.position = "none",
        panel.grid.major.y = element_line(linewidth = 0.1, color = "grey80"),
        panel.grid.minor.y = element_line(linewidth = 0.1, color = "grey95"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(color = "white")) +
  labs(title = "Response bar plot")

response_barchart
ggsave(filename = "Response_barchart.tiff",
       path = "data/Exploratory barcharts",
       width = 2*960, height = 2*960, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Treatment
treatment_barchart = ggplot(data = breast_cancer_full_pheno, aes(x = Treatment,
                                                                fill = Treatment))+
  scale_x_discrete(labels = c("Chemotherapy", "Endocrine"))+
  geom_bar(width = 0.35, color = "grey20", linewidth = 0.25)+
  geom_text(aes(label = after_stat(count)), vjust = 1.6, 
            color = "white", size = 2*1.5, stat = "count", fontface = "bold") +
  scale_fill_manual(name = "Group", values = c(Chemotherapy = "paleturquoise3",
                                               Endocrine_treatment = "sienna2"),
                    labels = c("Chemotherapy", "Endocrine"))+
  ylab("Number of samples")+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 2*5, hjust = 0.5),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 2*4.5),
        axis.text.x = element_text(size = 2*3.5, margin = ggplot2::margin(t = -3)),
        axis.text.y = element_text(size = 2*3.5),
        legend.position = "none",
        panel.grid.major.y = element_line(linewidth = 0.1, color = "grey80"),
        panel.grid.minor.y = element_line(linewidth = 0.1, color = "grey95"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(color = "white")) +
  labs(title = "Treatment bar plot")

treatment_barchart
ggsave(filename = "Treatment_barchart.tiff",
       path = "data/Exploratory barcharts",
       width = 2*960, height = 2*960, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# pam50
pam50_barchart = ggplot(data = breast_cancer_full_pheno, aes(x = pam50,
                                                            fill = pam50))+
  scale_x_discrete(labels = c("Basal", "HER2+", "Luminal A", "Luminal B", "Normal"))+
  geom_bar(width = 0.7, color = "grey20", linewidth = 0.25)+
  geom_text(aes(label = after_stat(count)), vjust = 1.6, 
            color = "white", size = 2*1.5, stat = "count", fontface = "bold") +
  scale_fill_manual(name = "Subtype", values = c(Basal = "red4", LumB = "skyblue",
                                                 LumA = "darkblue", Her2 = "purple4", 
                                                 Normal = "lightgreen"),
                    labels = c("Basal", "Luminal B", "Luminal A", "HER2+", "Normal"))+
  ylab("Number of samples")+
  xlab("pam50 subtype")+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 2*5, hjust = 0.5),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 2*4),
        axis.text.x = element_text(size = 2*3.5, margin = ggplot2::margin(t = -3)),
        axis.text.y = element_text(size = 2*3.5),
        legend.position = "none",
        panel.grid.major.y = element_line(linewidth = 0.1, color = "grey80"),
        panel.grid.minor.y = element_line(linewidth = 0.1, color = "grey95"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(color = "white")) +
  labs(title = "pam50 subtypes bar plot")

pam50_barchart
ggsave(filename = "pam50_barchart.tiff",
       path = "data/Exploratory barcharts",
       width = 2*960, height = 2*960, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Mammaprint risk
Mammaprint_risk_barchart = ggplot(data = breast_cancer_full_pheno, aes(x = Mammaprint_risk,
                                                                      fill = Mammaprint_risk))+
  scale_x_discrete(labels = c("No risk", "Risk"))+
  geom_bar(width = 0.35, color = "grey20", linewidth = 0.25)+
  geom_text(aes(label = after_stat(count)), vjust = 1.6, 
            color = "white", size = 2*1.5, stat = "count", fontface = "bold") +
  scale_fill_manual(name = "Subtype", values = c(`No risk` = "darkgreen",
                                                 Risk = "red3"),
                    labels = c("No risk", "Risk"))+
  ylab("Number of samples")+
  xlab("Mammaprint risk group")+
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 2*5, hjust = 0.5),
        plot.background = element_rect(fill = "white"),
        axis.title = element_text(face = "bold", size = 2*4),
        axis.text.x = element_text(size = 2*3.5, margin = ggplot2::margin(t = -3)),
        axis.text.y = element_text(size = 2*3.5),
        legend.position = "none",
        panel.grid.major.y = element_line(linewidth = 0.1, color = "grey80"),
        panel.grid.minor.y = element_line(linewidth = 0.1, color = "grey95"),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_line(color = "white")) +
  labs(title = "Mammaprint risk groups bar plot")

Mammaprint_risk_barchart
ggsave(filename = "Mammaprint_risk_barchart.tiff",
       path = "data/Exploratory barcharts",
       width = 2*960, height = 2*960, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Multiplots #####
# Vertical multiplot
ggarrange1 = ggarrange(sample_patient_barchart, labels = "A", ncol = 1,
                       font.label = list(size = 2*6, face = "bold"))
ggarrange2 = ggarrange(treatment_barchart, response_barchart, nrow = 1, ncol = 2,
                       labels = c("B", "C"), font.label = list(size = 2*6, face = "bold"))
ggarrange3 = ggarrange(pam50_barchart, Mammaprint_risk_barchart, nrow = 1, ncol = 2,
                       labels = c("D", "E"), font.label = list(size = 2*6, face = "bold"))
ggarrange(ggarrange1, ggarrange2, ggarrange3,
          ncol = 1, nrow = 3, widths = c(1, 1, 1), heights = c(2, 1, 1))
ggsave(filename = "Multiplot_barcharts_vertical.tiff",
       path = "data/Exploratory barcharts",
       width = 2120*2, height = 1920*4, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()

# Horizontal multiplot
hggarrange1 = ggarrange(sample_patient_barchart, labels = "A", ncol = 1,
                       font.label = list(size = 2*6, face = "bold"))
hggarrange2 = ggarrange(treatment_barchart, response_barchart, nrow = 1, ncol = 2,
                       labels = c("B", "C"), font.label = list(size = 2*6, face = "bold"))
hggarrange3 = ggarrange(pam50_barchart, Mammaprint_risk_barchart, nrow = 1, ncol = 2,
                       labels = c("D", "E"), font.label = list(size = 2*6, face = "bold"))
hggarrange_interm = ggarrange(hggarrange2, hggarrange3, ncol = 1, nrow = 2)
ggarrange(hggarrange1, hggarrange_interm,
          ncol = 2, nrow = 1)
ggsave(filename = "Multiplot_barcharts_horizontal.tiff",
       path = "data/Exploratory barcharts",
       width = 2120*4, height = 1920*2, device = 'tiff', units = "px",
       dpi = 700, compression = "lzw")
dev.off()