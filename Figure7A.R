
#Rithik Castelino (Walters Lab, NCI)

if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(ggrepel)) { install.packages("ggrepel"); library(ggrepel)}
if(!require(scales)) { install.packages("scales"); library(scales)}
if(!require(ggh4x)) { install.packages("ggh4x"); library(ggh4x)}
if(!require(plotly)) { install.packages("plotly"); library(plotly)}
if(!require(svglite)) { install.packages("svglite"); library(svglite)}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(version = "3.19")

#BiocManager::install("TissueEnrich")
library(TissueEnrich)

#See extended data for TMT-MS data used for this analysis and extract and save the sheet as a .csv. 
TMTMS_rawdata <- read.csv("Abundance_vs_pvalues.csv", skip = 1)
TMTMS_data <- TMTMS_rawdata %>% 
  dplyr::select(., 1, 4, 5) %>% 
  dplyr::rename(., Log2Abundance = Log2_.Ab) %>% 
  dplyr::rename(., Neglog10pvalue = Neg_log10..pvalue..1)

GeneNames <- unique(tail(TMTMS_data$Gene, -1))

gs <- GeneSet(geneIds=GeneNames, organism="Homo Sapiens", geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes=gs, rnaSeqDataset = 1, tissueSpecificGeneType = 1)

TissueSpecificGenes <- list()

for (i in 1:length(output[[3]])) {
  enrichmentOutput<-data.frame(assay(output[[3]][[i]]))$Gene
  TissueSpecificGenes <- c(TissueSpecificGenes, enrichmentOutput)
}

Uncategorized <- geneIds(output[[4]])

TMTMS_data_labelled <- TMTMS_data %>% 
  mutate(., TissueSpecificity = ifelse(Gene %in% TissueSpecificGenes, "Yes", "No")) %>% 
  mutate(., TissueSpecificity = ifelse(Gene %in% Uncategorized, "NA", TissueSpecificity))

graph <- ggplot(data=TMTMS_data_labelled, aes(x=Log2Abundance, y=Neglog10pvalue, colour = TissueSpecificity, label=Gene)) +
  theme_classic() +
  geom_point(data = TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$Neglog10pvalue >= 1.301),], size = 2.5) +
  geom_point(data = TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) < 1 & TMTMS_data_labelled$Neglog10pvalue >= 1.301),], size = 1) +
  geom_point(data = TMTMS_data_labelled[which(TMTMS_data_labelled$Neglog10pvalue < 1.301),], size = 1, colour = "grey") +
  geom_text_repel(data = TMTMS_data_labelled[which((TMTMS_data_labelled$Log2Abundance) > 1 & TMTMS_data_labelled$Neglog10pvalue > 1.301),], size = 3.5, max.overlaps = 3, min.segment.length = 1e-6, force_pull = 0.3, aes(family = "Arial Narrow"), segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20, nudge_y = 0.2, nudge_x = 0.2) +
  geom_text_repel(data = TMTMS_data_labelled[which((TMTMS_data_labelled$Log2Abundance) < -1 & TMTMS_data_labelled$Neglog10pvalue > 1.301),], size = 3.5, max.overlaps = 3, min.segment.length = 1e-6, force_pull = 0.3, aes(family = "Arial Narrow"), segment.size = unit(0.08, "cm"), box.padding = 0.2, segment.curvature = -0.3, segment.ncp = 10, segment.angle = 20, nudge_y = 0.2, nudge_x = -0.2) +
  geom_vline(xintercept = c(-1, 1), linetype="dotted", colour = "black") +
  geom_hline(yintercept = c(1.301), linetype="dotted", colour = "black") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, 1), guide="axis_minor", expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), breaks=seq(0, 7, 1), guide="axis_minor", expand = c(0,0)) +
  labs(x= "Log2(∆UIM/WT)", y="-Log10(p-value)") +
  theme(legend.position = "None") + 
  scale_colour_discrete(breaks=c("NA", "No", "Yes"), type = c("black", "cornflowerblue", "maroon"))
graph

ggsave("./Output/tissuespecificity_∆UIM.tiff", graph, width = 3456/250, height = 2234/(250/1.547), dpi = 300, bg = "transparent")
#ggsave("./Output/tissuespecificity.svg", graph, width = 3456/250, height = 2234/(250/1.547), dpi = 300, bg = "transparent")

ggplotly(graph)

Num_Specific_less1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) <= 1 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_less1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) <= 1 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_less1 = Num_Specific_less1/Num_Total_less1

Num_Specific_more1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_more1 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 1 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_more1 = Num_Specific_more1/Num_Total_more1

Num_Specific_more2 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 2 & TMTMS_data_labelled$TissueSpecificity == "Yes" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Num_Total_more2 <- nrow(TMTMS_data_labelled[which(abs(TMTMS_data_labelled$Log2Abundance) >= 2 & TMTMS_data_labelled$TissueSpecificity != "NA" & TMTMS_data_labelled$Neglog10pvalue > 1.301),])
Percent_Specific_Total_more2 = Num_Specific_more2/Num_Total_more2

#_____________________________________________
