
#Rithik Castelino (Walters Lab, NCI)

#Loading in Packages
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(magrittr)) { install.packages("magrittr"); library(magrittr)}

#Loading in Data Rep #1
C1_rawdata_rep1 <- read.csv("./Input/C1-26S_P1535.csv")
deltaUIM_rawdata_rep1 <- read.csv("./Input/UIM-26S_P1535.csv")
Combined_withGoTerms_rawdata_rep1 <- read.csv("./Input/P1535-Appendix3-C126S-vs-UIM26S-comparion-summary.csv")

#Fill in NA Values which by definition in this context are 0
Combined_withGoTerms_NAstoZero_rep1 <- Combined_withGoTerms_rawdata_rep1 %>% 
  mutate(., C126S...PSMs = ifelse(is.na(C126S...PSMs), 0, C126S...PSMs)) %>% 
  mutate(., UIM26S...PSMs = ifelse(is.na(UIM26S...PSMs), 0, UIM26S...PSMs))

#Filtering each group by certainty
Combined_withGoTerms_filtered_rep1 <- Combined_withGoTerms_NAstoZero_rep1 %>% 
  filter(!((.$C126S...PSMs > 0 | .$UIM26S...PSMs > 0) & (.$Accession %in% deltaUIM_rawdata_rep1[which(deltaUIM_rawdata_rep1$Exp..q.value..Combined > 0.05),]$Accession |
                                                           .$Accession %in% C1_rawdata_rep1[which(C1_rawdata_rep1$Exp..q.value..Combined > 0.05),]$Accession)))

#Loading & Filtering Data Rep #2
C1_rawdata_rep2 <- read.csv("./Input/C1-26S_P1499.csv")
deltaUIM_rawdata_rep2 <- read.csv("./Input/UIM-26S_P1511.csv")

C1_filtered_rep2 <- C1_rawdata_rep2 %>% 
  filter(.$Exp..q.value..Combined < 0.05) %>% 
  dplyr::select(Accession, contains("PSMs"), contains("kDa")) %>% 
  rename(C1_PSMs_rep2 = contains("PSMs"), mw_C1 = contains("kDa")) %>% 
  mutate(mw_C1 = mw_C1 * 1000)
  
deltaUIM_filtered_rep2 <- deltaUIM_rawdata_rep2  %>% 
  filter(.$Exp..q.value..Combined < 0.05) %>% 
  select(Accession, contains("PSMs"), contains("kDa")) %>% 
  rename(deltaUIM_PSMs_rep2 = contains("PSMs"), mw_deltaUIM = contains("kDa")) %>% 
  mutate(mw_deltaUIM = mw_deltaUIM * 1000)

#Filtering Combined for Rep #1 by common proteins in Rep #2
Combined_withGoTerms_filtered_commonalities <- Combined_withGoTerms_filtered_rep1 %>% 
  filter(.$Accession %in% C1_filtered_rep2$Accession | .$Accession %in% deltaUIM_filtered_rep2$Accession)

#Renaming Columns to Distinguish Rep #1 and Bringing Values from Rep #2
Combined_withGoTerms_filtered_commonalities_merged <- Combined_withGoTerms_filtered_commonalities %>% 
  rename(C1_PSMs_rep1 = "C126S...PSMs", deltaUIM_PSMs_rep1 = "UIM26S...PSMs") %>% 
  left_join(C1_filtered_rep2, by="Accession") %>% 
  left_join(deltaUIM_filtered_rep2, by="Accession") %>% 
  replace(., is.na(.), 0) %>% 
  mutate(mw = pmax(mw_C1, mw_deltaUIM))


#DEFINING FUNCTIONS 
#Creating template with conversions between GenSymbol and Gene Names
if (file.exists('./Input/Accession_GeneName_Mapping.csv')) {
  Accession_GeneName_Mapping <- read.csv('./Input/Accession_GeneName_Mapping.csv')
  Accession_GeneName_Mapping <- Accession_GeneName_Mapping %>% 
    mutate(., Gene_Name = ifelse((Gene_Name == ""), NA, Gene_Name))
} else {
  Accession_GeneName_Mapping = data.frame(Accession = c(), Gene_Name = c()) 
  for (i in 1:(length(CombinedComparison_rawdata$Accession))) {
    print(i)
    Accession_i <- CombinedComparison_rawdata$Accession[i]
    Accession_GeneName_Mapping_GeneName_i <- ConvertID(Accession_i, ID_from = "UniProtKB_AC-ID", ID_to = "Gene_Name")
    
    if ((Accession_GeneName_Mapping_GeneName_i == "Resource not found" | Accession_GeneName_Mapping_GeneName_i == "Internal server error")) {
      Accession_GeneName_Mapping_GeneName_i <- data.frame(Gene_Name = c(NA))
    } else {
      Accession_GeneName_Mapping_i = data.frame(Accession = c(Accession_i), Gene_Name = c(Accession_GeneName_Mapping_GeneName_i)) 
    }
    Accession_GeneName_Mapping <- rbind(Accession_GeneName_Mapping, Accession_GeneName_Mapping_i)
  }
  write.csv(Accession_GeneName_Mapping, "./Input/Accession_GeneName_Mapping.csv")
}

#FUNCTION: Conversion from GenSymbol to Gene Names
Conversion_to_GeneName <- function(Data) {
  Data <- Data %>% 
    mutate(., GeneName = "")
    
  for (i in 1:length(Data$Accession)) {
    if (Data$Accession[i] %in% Accession_GeneName_Mapping$Accession) {
      Data$GeneName[i] <- Accession_GeneName_Mapping[which(Accession_GeneName_Mapping$Accession == Data$Accession[i]),]$Gene_Name
      print(Accession_GeneName_Mapping[which(Accession_GeneName_Mapping$Accession == Data$Accession[i]),]$Gene_Name)
    } else {
      Data$GeneName[i] <- ""
    }
  }
  
  Data <- Data %>% 
    mutate(., WaltersName = .$GeneName)
  
  Output <- Data %>% 
    relocate(GeneName, .before = 2) %>% 
    relocate(WaltersName, .before = 1) %>% 
    mutate(., WaltersName = ifelse(.$Accession == "P62979", "UB ", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$Accession == "P61221", "ABCE1", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD3", "Rpn3", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD1", "Rpn2", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD12", "Rpn5", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD11", "Rpn6", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD6", "Rpn7", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD7", "Rpn8", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD13", "Rpn9", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD14", "Rpn11", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD8", "Rpn12", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD9", "Rpn15", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD2", "Rpn1", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD4", "Rpn10", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "ADRM1", "Rpn13", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC2", "Rpt1", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC1", "Rpt2", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC4", "Rpt3", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC6", "Rpt4", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC3", "Rpt5", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMC5", "Rpt6", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA6", "α1", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA2", "α2", .$WaltersName)) %>%
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA4", "α3", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA7", "α4", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA5", "α5", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA1", "α6", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMA3", "α7", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB6", "β1", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB9", "β1i", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB7", "β2", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB10", "β2i", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB3", "β3", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB2", "β4", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB5", "β5", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB8", "β5i", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB1", "β6", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB4", "β7", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD5", "S5b", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD9", "p27", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMD10", "p28", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "PSMB1", "β6", .$WaltersName)) %>% 
    mutate(., WaltersName = ifelse(.$GeneName == "UBE3A", "E6AP", .$WaltersName))
  return(Output)
  
}

#FUNCTION: Including type/function
IncludingType <- function(Data) {
  Output <- Data %>% 
    mutate(., Type = ifelse(((grepl("proteasome binding", .$Gene.ontology..molecular.function., ignore.case=TRUE))), "Proteasome Accessory Proteins", "Other")) %>% 
    mutate(., Type = ifelse((grepl("chaperone", .$Function..CC., ignore.case=TRUE)), "Chaperone", Type)) %>% 
    mutate(., Type = ifelse(((grepl("proteasome", .$Protein.names, ignore.case=TRUE) &
                                grepl("subunit", .$Protein.names, ignore.case=TRUE)) |
                               grepl("Proteasomal ubiquitin receptor", .$Protein.names, ignore.case=TRUE)), "20S Proteasome Subunit", Type)) %>% 
    mutate(., Type = ifelse((grepl("20S Proteasome Subunit", .$Type, ignore.case=TRUE) & grepl("Rp", .$WaltersName, ignore.case=TRUE)), "19S Proteasome Subunit", Type)) %>% 
    mutate(., Type = ifelse((grepl("Acts as a chaperone during the assembly of the 26S proteasome", .$Function..CC., ignore.case=TRUE) |
                               grepl("Inhibits proteasome 26S assembly", .$Function..CC., ignore.case=TRUE)), "Proteasome Chaperone", Type)) %>% 
    mutate(., Type = ifelse((grepl("40S ribosomal protein", .$Description, ignore.case=TRUE) &
                               !grepl("UB", .$WaltersName, ignore.case=TRUE)), "40S Ribosome Subunit", Type)) %>% 
    mutate(., Type = ifelse((grepl("60S ribosomal protein", .$Description, ignore.case=TRUE) &
                               !grepl("UB", .$WaltersName, ignore.case=TRUE)), "60S Ribosome Subunit", Type)) %>% 
    mutate(., Type = ifelse((grepl("UB ", .$WaltersName, ignore.case=TRUE)), "Ubiquitin", Type)) %>%
    mutate(., WaltersName = ifelse((grepl("UB ", .$WaltersName, ignore.case=TRUE)), "UB", WaltersName)) %>%
    mutate(., Type = ifelse((grepl("TXNL1", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("UBE3C", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("E6AP", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("UBQLN", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("RAD23A", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("UBR4", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("DDI1", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("DDI2", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    mutate(., Type = ifelse((grepl("ZFAND5", .$WaltersName, ignore.case=TRUE)), "Proteasome Accessory Proteins", Type)) %>% 
    relocate(Type, .before = 5) 
}


#APPLYING FUNCTIONS & Further Processing
#REP #1
#Converting GenSymbols to Gene Names
Combined_GeneNames <- Conversion_to_GeneName(Combined_withGoTerms_filtered_commonalities_merged)

#Scaling by MW
Combined_GeneNames_scaled <- Combined_GeneNames %>% 
  mutate(C1_scaled_rep1 = C1_PSMs_rep1/(mw)) %>% 
  mutate(C1_scaled_rep2 = C1_PSMs_rep2/(mw)) %>% 
  mutate(deltaUIM_scaled_rep1 = deltaUIM_PSMs_rep1/(mw)) %>% 
  mutate(deltaUIM_scaled_rep2 = deltaUIM_PSMs_rep2/(mw))

Combined_GeneNames_scaled <- Combined_GeneNames_scaled %>% 
  mutate(deltaUIM_scaled_rep1 = ifelse(WaltersName == "Rpn10", deltaUIM_scaled_rep1*(40737.66/21151.30), deltaUIM_scaled_rep1)) %>% 
  mutate(deltaUIM_scaled_rep2 = ifelse(WaltersName == "Rpn10", deltaUIM_scaled_rep2*(40737.66/21151.30), deltaUIM_scaled_rep2))

Combined_GeneNames_normalized <- Combined_GeneNames_scaled %>% 
  mutate(C1_norm_rep1 = C1_scaled_rep1/.[which(.$WaltersName == "Rpn1"),]$C1_scaled_rep1) %>% 
  mutate(C1_norm_rep2 = C1_scaled_rep2/.[which(.$WaltersName == "Rpn1"),]$C1_scaled_rep2) %>% 
  mutate(deltaUIM_norm_rep1 = deltaUIM_scaled_rep1/.[which(.$WaltersName == "Rpn1"),]$deltaUIM_scaled_rep1) %>% 
  mutate(deltaUIM_norm_rep2 = deltaUIM_scaled_rep2/.[which(.$WaltersName == "Rpn1"),]$deltaUIM_scaled_rep2)

#Averaging Values Together
Combined_averages <- Combined_GeneNames_normalized %>% 
  mutate(C1_mean = rowMeans(select(., C1_norm_rep1, C1_norm_rep2))) %>% 
  mutate(deltaUIM_mean = rowMeans(select(., deltaUIM_norm_rep1, deltaUIM_norm_rep2)))

#Adding tag based on cell line variation
Combined_averaged_celllinevariation <- Combined_averages %>% 
  mutate(Tag = ifelse(.$deltaUIM_mean > 0 & .$C1_mean == 0, "deltaUIM_unique", "tbd")) %>% 
  mutate(Tag = ifelse(.$deltaUIM_mean == 0 & .$C1_mean > 0, "C1_unique", .$Tag)) %>% 
  mutate(Tag = ifelse(.$deltaUIM_mean > 0 & .$C1_mean > 0, "PresentinBoth", .$Tag))

#--------------------------------------------------------------------------------------------------------------------------------------------------------
#GRAPHING
#Loading in additional packages for graph aesthetics
if(!require(ggrepel)) { install.packages("ggrepel"); library(ggrepel)}
if(!require(scales)) { install.packages("scales"); library(scales)}
if(!require(ggh4x)) { install.packages("ggh4x"); library(ggh4x)}
if(!require(plotly)) { install.packages("plotly"); library(plotly)}
if(!require(cowplot)) { install.packages("cowplot"); library(cowplot)}

#Comparing WT & ∆UIM
Combined_averaged_celllinevariation_type <- IncludingType(Combined_averaged_celllinevariation)
Combined_averaged_celllinevariation_type$Type <- factor(Combined_averaged_celllinevariation_type$Type, levels = c("20S Proteasome Subunit", "19S Proteasome Subunit", "Proteasome Accessory Proteins", "Proteasome Chaperone", "Ubiquitin", "60S Ribosome Subunit", "40S Ribosome Subunit", "Chaperone", "Other"))

data_toplot <- Combined_averaged_celllinevariation_type[which(Combined_averaged_celllinevariation_type$Tag == "PresentinBoth"),]

graph <- ggplot(data = data_toplot, aes(x=C1_mean, y = deltaUIM_mean, label=WaltersName, colour=Type)) +
  geom_point(data = data_toplot[which(data_toplot$Type == "Other"),], size = 0.6, key_glyph = "point") +
  geom_point(data = data_toplot[which(data_toplot$Type != "Other"),], size = 0.6, key_glyph = "point") +
  geom_text_repel(data = data_toplot[which(data_toplot$Type == "Other" & (data_toplot$deltaUIM_mean > 0 & data_toplot$C1_mean > 0)),],size = 2, max.overlaps = 10, min.segment.length = 1e-10, force_pull = 1, aes(family = "Arial Narrow"), nudge_y = 0, nudge_x = 0, show.legend = FALSE, segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20) +
  theme_classic() +
  guides(colour = guide_legend(ncol = 2)) +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(0, 2)) +
  labs(x = "C1 (PSM/MW)", y = "ΔUIM (PSM/MW)") +
  theme(axis.text = element_text(size=5, colour = "black"), 
        axis.title = element_text(size=5, colour = "black"), 
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.length = unit(0.087, "cm"),
        axis.ticks = element_line(size = 0.087),
        legend.position = "none", legend.text = element_blank(), legend.title = element_blank(), legend.key.spacing.y = unit(1, 'cm'), legend.key.spacing.x = unit(1, 'cm'),
        text = element_text(family="Arial Narrow"),
        plot.background = element_blank(), panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill='transparent'))  +
  scale_colour_discrete(breaks=c("20S Proteasome Subunit", "19S Proteasome Subunit", "Proteasome Accessory Proteins", "Proteasome Chaperone", "Ubiquitin", "60S Ribosome Subunit", "40S Ribosome Subunit", "Chaperone", "Other"), type = c("grey", "purple", "black", "deeppink3", "brown", "red", "darkgreen", "lightgreen", "darkblue"))
graph

ggsave("./Output/PresentinBoth_ScaledbyMW_rec.svg", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")
ggsave("./Output/PresentinBoth_ScaledbyMW_rec.tiff", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")

#Plotting proteins only found in the C1/WT Sample
C1_unique_type <- IncludingType(Combined_averaged_celllinevariation[which(Combined_averaged_celllinevariation_type$Tag == "C1_unique"),])
C1_unique_type$Type <- factor(C1_unique_type$Type, levels = c("20S Proteasome Subunit", "19S Proteasome Subunit", "Proteasome Accessory Proteins", "Proteasome Chaperone", "Ubiquitin", "60S Ribosome Subunit", "40S Ribosome Subunit", "Chaperone", "Other"))
C1_unique_type_toorder <- C1_unique_type$C1_mean
C1_unique_type_toorder_orders <- rank(-C1_unique_type_toorder, ties.method="random")
C1_unique_type$rank <- C1_unique_type_toorder_orders

data_toplot <- C1_unique_type
graph <- ggplot(data = data_toplot, aes(x= rank, y = C1_mean, label=WaltersName, colour=Type)) +
  geom_point(data = data_toplot[which(data_toplot$Type == "Other"),], size = 0.6, key_glyph = "point") +
  geom_point(data = data_toplot[which(data_toplot$Type != "Other"),], size = 0.6, key_glyph = "point") +
  geom_text_repel(data = data_toplot[which(data_toplot$Type != "Other"),],size = 2, max.overlaps = Inf, min.segment.length = 1e-10, force_pull = 1, aes(family = "Arial Narrow"), nudge_y = 0.01, nudge_x = 0.01, show.legend = FALSE, segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20) +
  theme_classic() +
  scale_x_continuous(limits = c(1, 100), breaks=seq(1, 100, 20), minor_breaks=seq(1, 100, 10), guide="axis_minor", expand = c(0,1)) +
  scale_y_continuous(limits = c(0, 0.15), breaks=seq(0, 0.15, 0.03), minor_breaks=seq(0, 0.15, 0.015), guide="axis_minor", expand = c(0,0)) +
  labs(x = "Protein Rank", y = "C1 PSMs") +
  theme(axis.text = element_text(size=5, colour = "black"), 
        axis.title = element_text(size=5, colour = "black"), 
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.length = unit(0.087, "cm"),
        axis.ticks = element_line(size = 0.087),
        legend.position = c(0.75, 0.9), legend.text = element_text(size = 5, margin=margin(l=-1, r=-1)), legend.title = element_blank(), legend.key.spacing.y = unit(-0.4, 'cm'), legend.key.spacing.x = unit(-0.65, 'cm'),
        text = element_text(family="Arial Narrow"),
        plot.background = element_blank(), panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill='transparent'))  +
  scale_colour_discrete(breaks=c("Proteasome Accessory Proteins", "40S Ribosome Subunit", "Chaperone", "Other"), type = c("grey", "deeppink3", "lightgreen", "darkblue"))

graph
ggsave("./Output/C1_uniques.svg", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")
ggsave("./Output/C1_uniques.tiff", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")

#Plotting proteins only found in the ∆UIM Sample
deltaUIM_unique_type <- IncludingType(Combined_averaged_celllinevariation[which(Combined_averaged_celllinevariation_type$Tag == "deltaUIM_unique"),])
deltaUIM_unique_type$Type <- factor(deltaUIM_unique_type$Type, levels = c("20S Proteasome Subunit", "19S Proteasome Subunit", "Proteasome Accessory Proteins", "Proteasome Chaperone", "Ubiquitin", "60S Ribosome Subunit", "40S Ribosome Subunit", "Chaperone", "Other"))
deltaUIM_unique_type_toorder <- deltaUIM_unique_type$deltaUIM_mean
deltaUIM_unique_type_toorder_orders <- rank(-deltaUIM_unique_type_toorder, ties.method="random")
deltaUIM_unique_type$rank <- deltaUIM_unique_type_toorder_orders

data_toplot <- deltaUIM_unique_type

graph <- ggplot(data = data_toplot, aes(x= rank, y = deltaUIM_mean, label=WaltersName, colour=Type)) +
  geom_point(data = data_toplot[which(data_toplot$Type == "Other"),], size = 0.6, key_glyph = "point") +
  geom_point(data = data_toplot[which(data_toplot$Type != "Other"),], size = 0.6, key_glyph = "point") +
  geom_text_repel(data = data_toplot[which(data_toplot$Type != "Other"),],size = 2, max.overlaps = Inf, min.segment.length = 1e-10, force_pull = 1, aes(family = "Arial Narrow"), nudge_y = 0.01, nudge_x = 0.01, show.legend = FALSE, segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20) +
  theme_classic() +
  scale_x_continuous(limits = c(1, 10), breaks=seq(1, 10, 2), minor_breaks=seq(1, 10, 1), guide="axis_minor", expand = c(0,1)) +
  scale_y_continuous(limits = c(0, 0.011), breaks=seq(0, 0.011, 0.01), minor_breaks=seq(0, 0.011, 0.005), guide="axis_minor", expand = c(0,0)) +
  labs(x = "Protein Rank", y = "ΔUIM PSMs") +
  theme(axis.text = element_text(size=5, colour = "black"), 
        axis.title = element_text(size=5, colour = "black"), 
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.length = unit(0.087, "cm"),
        axis.ticks = element_line(size = 0.087),
        legend.position = c(0.75, 0.9), legend.text = element_text(size = 5, margin=margin(l=-1, r=-1)), legend.title = element_blank(), legend.key.spacing.y = unit(-0.4, 'cm'), legend.key.spacing.x = unit(-0.65, 'cm'),
        text = element_text(family="Arial Narrow"),
        plot.background = element_blank(), panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill='transparent'))  +
  scale_colour_discrete(breaks=c("Chaperone", "Other"), type = c("grey", "darkblue"))

graph
ggsave("./Output/deltaUIMs_uniques.svg", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")
ggsave("./Output/deltaUIMs_uniques.tiff", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")

#Saving output as a .csv for future viewing
write.csv(C1_unique_type, './Output/C1_unique_filtered_GeneNames.csv')
write.csv(deltaUIM_unique_type, './Output/deltaUIMs_unique_filtered_GeneNames.csv')
write.csv(Combined_averaged_celllinevariation_type, './Output/PresentinBoth_unique_filtered_GeneNames.csv')



Combined_averaged_celllinevariation_type <- IncludingType(Combined_averaged_celllinevariation)
Combined_averaged_celllinevariation_type$Type <- factor(Combined_averaged_celllinevariation_type$Type, levels = c("20S Proteasome Subunit", "19S Proteasome Subunit", "Proteasome Accessory Proteins", "Proteasome Chaperone", "Ubiquitin", "60S Ribosome Subunit", "40S Ribosome Subunit", "Chaperone", "Other"))
data_toplot <- Combined_averaged_celllinevariation_type %>% 
  select(WaltersName, Type, C1_norm_rep1, C1_norm_rep2, deltaUIM_norm_rep1, deltaUIM_norm_rep2) %>% 
  pivot_longer(!c(WaltersName, Type), names_to = "CellLine", values_to = "Ratio") %>% 
  #filter(.$Type =="Chaperone") %>% 
  filter(.$Type =="Proteasome Accessory Proteins" | .$Type =="Proteasome Chaperone") %>% 
  mutate(., Ratio = ifelse(.$Ratio == 0, NA, .$Ratio))

graph <- ggplot(data = data_toplot, aes(y=(reorder(WaltersName, Ratio)), x=CellLine, fill=Ratio)) +
  theme_classic() +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_viridis_c(option = "magma", direction = 1, na.value = "grey") +
  #scale_fill_gradient2(low = 'white', mid = "darkgrey", high='brown', na.value = "black") +
  #scale_fill_gradient2(low = 'white', mid = "darkgrey", high='brown', midpoint= max(data_toplot$Ratio)/2) +
  labs(x="Fraction") +
  theme(axis.text = element_text(size=5, colour = "black"), 
        axis.title = element_text(size=5, colour = "black"), 
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        axis.line = element_line(size = 0.25, colour = "black"),
        axis.ticks.length = unit(0.087, "cm"),
        axis.ticks = element_line(size = 0.087),
        text = element_text(family="Arial Narrow"),
        plot.background = element_blank(), panel.background = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill='transparent'))
graph

ggsave("./Output/ScaledbyMW_HEATMAP_rec.svg", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")
ggsave("./Output/ScaledbyMW_HEATMAP_rec.tiff", graph, width = 3456/1000, height = 2234/(1000/1.547), dpi = 1600, bg = "transparent")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
