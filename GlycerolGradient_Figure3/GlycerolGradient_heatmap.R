
#GlycerolGradient_Heatmap Script
#Rithik Castelino (Walters Lab, NCI)
#Designed for ∆UIM Manuscript, taking in pixel intensity of western blots for varying proteasomal subunits after HCT116 Cell Lysates (WT, ∆UIM, and ∆RAZUL)
#were separated by a glycerol gradient.

#Loading in packages
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(ggpubr)) { install.packages("ggpubr"); library(ggpubr)}

#Loading in Data
raw_data <- read.csv("Input/20241112_GlycerolGradient_Figure3_Quantified.csv")

#Subtracting Background & Finding % of Lane
Processed_data_ECPAS_WT <- raw_data %>% 
  filter(., (CellLine == "WT" & Protein == "beta5" & Type == "Background")) %>% 
  inner_join(., filter(raw_data, (CellLine == "WT" & Protein == "beta5" & Type == "Sample")), by="Lane") %>% 
  mutate(., Difference = abs(Value.x - Value.y)) %>% 
  mutate(., Population = Difference/sum(Difference))

Processed_data_ECPAS_deltaUIM <- raw_data %>% 
  filter(., (CellLine == "deltaUIM" & Protein == "beta5" & Type == "Background")) %>% 
  inner_join(., filter(raw_data, (CellLine == "deltaUIM" & Protein == "beta5" & Type == "Sample")), by="Lane") %>% 
  mutate(., Difference = abs(Value.x - Value.y)) %>% 
  mutate(., Population = Difference/sum(Difference))

Processed_data_ECPAS <- rbind(Processed_data_ECPAS_WT, Processed_data_ECPAS_deltaUIM)

Data_toplot <- Processed_data_ECPAS %>% 
  mutate(., CellLine = factor(CellLine.x, levels = c("deltaUIM", "WT")))

ECPAS <- ggplot(data = Data_toplot, aes(y=CellLine.x, x=Lane, fill=Population)) +
  theme_classic() +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = '#fdfefe', mid = "#fdca69", high="#ff0000", midpoint=0.13) +
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
ECPAS
ggsave("./Output/beta5.svg", ECPAS)
ggsave("./Output/beta5.tiff", ECPAS)
