################################################################################
#MANUAL VS AUTOMATIC - Density Validation
################################################################################
#load libraries
library(extrafont)
font_import() #need to do once at beginning before the rest of the script can run
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(stats)
library(RColorBrewer)
library(stringr)

#inputs
Cell_Type <- "mESCs"
New_Line_Name <- "WT_auto"
File_Path <- "C:\\Users\\hmroa\\Documents\\RNA-SPLIT_Results\\Manual_vs_Auto\\Pulse_Chase_Analysis"

################################################################################
#STEP 1: load manual and auto WT centroid data which has had specific data removed

setwd(paste(File_Path, Cell_Type, "Density", "All_Cell_Lines_Merged", sep="\\"))

#loads original WT files
WT_Original <- read_csv("Cleaned_Other_Cell_Lines_Density_Compile.csv") %>%
  filter(Cell_Line == "WT") %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_original"))

#loads manual compile files
WT_Manual_Line <- read_csv("Manual_Compile_Density_Data.csv")

#loads automated data
setwd(paste(File_Path, Cell_Type, "Density", New_Line_Name,sep="\\"))

WT_Auto_Line <- read_csv(paste(New_Line_Name, "Density_Compile.csv", sep="_")) %>%
  filter(Data_Set == "Dynamic")
  

################################################################################
#STEP 2:DENSITY ANALSYIS - plot the data

#combine auto and manual data to orignal
Present_Density <- WT_Auto_Line %>%
  select(-Cloud_Name, -Time) %>%
  bind_rows(WT_Manual_Line) %>% 
  select(Cell_Line, Distance) %>%
  #bind_rows(WT_Original) %>%
  mutate(X = "")

#order data to be presented in box plot 
Present_Density$Cell_Line <- factor(Present_Density$Cell_Line, 
                                 levels = c(#"WT_original",
                                            "WT_manual", "WT_auto"), ordered = TRUE) #use names defined in line 51

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri", size = (11)), 
               legend.text = element_text(family = "Calibri", size = (11)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (9)))


#make box plot
density_plot <- ggplot(Present_Density, aes(x = X, y = Distance, fill = Cell_Line)) + #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                        #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1) +                                           #shows box plots for selected cell lines
  scale_fill_manual(values = c(#"#F768A1",                                      #WT_original exp and SS colours 
                               "#FB6A4A",                                      #WT_manual exp and SS colours 
                               "#969696")) +                                   #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((250),(1300) )) +                                   #adjusts scale so whiskers don't touch the end graph 
  labs(title = "B) Validate Density Analysis",                                             #sets name of the axes
       x = "",
       y = "Distance between Xist RNPs [nm]") +                                          
  theme                                                                        #adds the theme which defines the font, text size etc


#view the box plot
density_plot

################################################################################
#STEP 2d: TOTAL MOLECULE COUNT ANALSYIS 
#preform t Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff <- t.test(Present_Density[which(Present_Density$Cell_Line==Refernece_Line),]$Distance,
                         Present_Density[which(Present_Density$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Density[which(Present_Density$Cell_Line==Refernece_Line),]$Distance,
                         Present_Density[which(Present_Density$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Density[which(Present_Density$Cell_Line==Refernece_Line),]$Distance,
                         Present_Density[which(Present_Density$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Both_Phases_WT_vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- "Overall, the automated compile has statistically greater distances between Xist RNPs (lower density) than the manual compile"
  } else {
    outcome_1 <- "Overall, the new automated compile has statistically shorter distances between Xist RNPs (higher density) than the manual compile"
  }} else {
    outcome_1 <- "There is no statistical difference in Xist RNP distances (density) between the automated compile and the manual compile"
  }


#save p values

setwd(paste(File_Path, Cell_Type, "Density", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Density_t-test_p_values.csv", sep="_"))

#return outcome of the stats test
print(outcome_1)

