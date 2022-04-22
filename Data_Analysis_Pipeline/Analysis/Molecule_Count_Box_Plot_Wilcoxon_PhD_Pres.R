#load libraries 

#run these two lines on their own first
library(extrafont)
#font_import() 

#finish running rest of the script from this line
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(RColorBrewer)

################################################################################
#load western blot data to be analysed

file_path <- choose.files(default = "", caption = "Select WB_DATA_R_Analysis_CSV file",
                          multi = TRUE, filters = Filters,
                          index = nrow(Filters))

save_path <- choose.dir(default = "", caption = "Select p-values folder")

All_Molecule_Count_Data <- read_csv(file_path)


#######################
#SELECT DATA TO PRESENT

#set name of new cell line, and up to 3 others to present in the plot
#all names should be spelled the same as file names within the directory
#can leave names empty if don't want to present them all
New_Line_Name <- "Mettl3_dTAG"       

Refernece_Line <- "WT"        #reference_line will be used to perform stats test with

Cell_Line_3 <- "Ciz1_KO"

Cell_Line_4 <- "SPEN_RRM_del"  

#define the threshold used to determine a significant result during wilcoxon test
p <- 0.05


#######################
#DEFINE ASTHETICS OF PLOT

#may need to adjust scale for y-axis on line 

#define colours for each cell line - leave blank if not presenting all 4 cell lines
Colour_1_exp <- "#FA9FB5" #sets colour for expansion phase new reference_line
Colour_1_ss <- "#DD3497"  #sets colour for steady state phase new reference_line

Colour_2_exp <- "#A6BDDB" #sets colour for expansion phase new cell line
Colour_2_ss <- "#3690C0"  #sets colour for steady state phase new cell line

Colour_3_exp <- "#FEB24C" #sets colour for expansion phase new cell_line_3
Colour_3_ss <- "#FC4E2A"  #sets colour for steady state phase new cell_line_3

Colour_4_exp <- "#ADDD8E" #sets colour for expansion phase new cell_line_3
Colour_4_ss <- "#41AB5D"  #sets colour for steady state phase new cell_line_3



#set name of how cell lines should be presented in the plot - ensures name consistency with other papers
#if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
Name_1 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line
Name_2 <- bquote("WT")                   #sets name for reference_line
Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
Name_4 <- bquote("SPEN"^~Delta*"RRM")    #sets name for cell_line_4

#set name of plot titles
Title <- "Number of Xist RNPs"
x_axis <- "Phase"
y_axis <- "Total Xist RNA [no. centroids]"

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of box plots

#chose which data to present in box plot
Present_Data <- All_Molecule_Count_Data %>%
  filter(Cell_Line %in%
           c(Refernece_Line, New_Line_Name #, Cell_Line_3, Cell_Line_4
             ) &  #select which cell lines to present
           Phase %in%
           c(#"Initiation", 
             "Maintenance")) %>%                            #select which Phases to present
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%                 
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                         #creates key variable to base colour coding off


#order data to be presented in box plot 
Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Refernece_Line, New_Line_Name, Cell_Line_3, Cell_Line_4), ordered = TRUE,
                                 labels = c(Name_2, Name_1, Name_3, Name_4)) 

#put in order for the graphs to appear on the x-axis
Present_Data$Key <- factor(Present_Data$Key, 
                           levels = c(paste("Expansion", Refernece_Line, sep="_"), paste("Steady_State", Refernece_Line, sep="_"),
                                      paste("Expansion", New_Line_Name, sep="_"), paste("Steady_State", New_Line_Name, sep="_"),
                                      paste("Expansion", Cell_Line_3, sep="_"), paste("Steady_State", Cell_Line_3, sep="_"),
                                      paste("Expansion", Cell_Line_4, sep="_"), paste("Steady_State", Cell_Line_4, sep="_"), ordered = TRUE))


#################################################################################
#STEP 3: PLOT THE DATA

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (20)),
               legend.title = element_text(family = "Calibri", size = (16)), 
               legend.text = element_text(family = "Calibri", size = (14)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (16)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (12)),
               strip.text.x = element_text(family = "Calibri", size = (16)))

#calculate the y scale for the plot
Min_y = min(Present_Data$Total_Centroids, na.rm = TRUE) 
Max_y = max(Present_Data$Total_Centroids, na.rm = TRUE)  

box_plot <- ggplot(Present_Data, aes(x= Cell_Line, y = Total_Centroids, fill = Key)) + #defines data to present
  geom_boxplot(outlier.shape = NA) +                                               #defines type of plot
  #facet_wrap(~Cell_Line, nrow = 1, labeller = label_parsed) +                      #shows box plots for selected cell lines
  scale_fill_manual(values = c(Colour_1_exp, #Colour_1_ss,                          #defines the colour of each boxplot 
                               Colour_2_exp, Colour_2_ss,                                               
                               Colour_3_exp, Colour_3_ss,                                                  
                               Colour_4_exp, Colour_4_ss)) +                                                     
  coord_cartesian(ylim = c((Min_y),(200) )) +                                      #adjusts scale so whiskers don't touch the end graph 
  labs(title = Title,                                                              #sets name of the axes
       x = x_axis,
       y = y_axis) +                                          
  theme   

#view the box plot
box_plot
