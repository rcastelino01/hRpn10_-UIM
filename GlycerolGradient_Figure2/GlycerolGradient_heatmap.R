
#GlycerolGradient_Heatmap Script
#Rithik Castelino (Walters Lab, NCI)
#Designed for ∆UIM Manuscript, taking in pixel intensity of western blots for varying proteasomal subunits after HCT116 Cell Lysates (WT, ∆UIM, and ∆RAZUL)
#were separated by a glycerol gradient.

#Loading in packages
if(!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse)}
if(!require(ggpubr)) { install.packages("ggpubr"); library(ggpubr)}

#Loading in Data
raw_data <- read.csv("Input/20241020_GlycerolGradient_Figure2.csv")

#Finding % of Lane (Background already subtracted in 'raw_data')
Processed_data_Rpn10_WT <- raw_data %>% 
  filter(., (CellLine == "WT" & Protein == "Rpn10")) %>% 
  mutate(., Population = Value/sum(Value))

Processed_data_Rpn10_deltaRAZUL <- raw_data %>% 
  filter(., (CellLine == "Razul" & Protein == "Rpn10")) %>% 
  mutate(., Population = Value/sum(Value))

Processed_data_Rpn10_deltaUIM <- raw_data %>% 
  filter(., (CellLine == "UIM" & Protein == "Rpn10")) %>% 
  mutate(., Population = Value/sum(Value))

Processed_data_ECPAS <- rbind(Processed_data_Rpn10_WT, Processed_data_Rpn10_deltaUIM, Processed_data_Rpn10_deltaRAZUL)

Data_toplot <- Processed_data_ECPAS %>% 
  mutate(., CellLine = factor(CellLine, levels = c("UIM", "Razul", "WT")))

Rpn10 <- ggplot(data = Data_toplot, aes(y=CellLine, x=Lane, fill=Population)) +
  theme_classic() +
  geom_tile(color = "black", linewidth = 0.5) +
  #scale_fill_gradient2(low = '#fdfefe', mid = "#fdca69", high="#ff0000", midpoint=0.17) +
  scale_fill_viridis_c(option = "magma", direction = -1, na.value = "white") +
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
Rpn10
ggsave("./Output/Rpn8.svg", Rpn10)
ggsave("./Output/Rpn8.tiff", Rpn10)
