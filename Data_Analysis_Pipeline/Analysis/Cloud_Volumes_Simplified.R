#################################################################################
#SRCIPT CREATES BOX PLOT FOR CLOUD VOLUME DATA AND PREFORMS WILCOX STAT TEST
#################################################################################
#load libraries 

#run these two lines on their own first
library(extrafont)
#font_import() 

#finish running rest of the script from this line
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(RColorBrewer)
library(stringr)


################################################################################
#USER INPUT REQUIRED - check lines 21-63

#define type of cells used in experiment
Cell_Type <- "mESCs"                     #either "mESCs" or "NPCs"

#######################
#SELECT DATA TO PRESENT

#set name of new cell line, and up to 3 others to present in the plot
#all names should be spelled the same as file names within the directory
#can leave names empty if don't want to present them all
New_Line_Name <- "Mettl3_dTAG"       

Refernece_Line <- "WT"        #reference_line will be used to perform stats test with

#Cell_Line_3 <- "Ciz1_KO"

Cell_Line_4 <- "SPEN_RRM_del"  

#define the threshold used to determine a significant result during wilcoxon test
p <- 0.05


#######################
#DEFINE ASTHETICS OF PLOT

#may need to adjust scale for y-axis on line 

#define colours for each cell line - leave blank if not presenting all 4 cell lines
Colour_1 <- "#F768A1" #sets colour for new cell line
Colour_2 <- "#67A9CF" #sets colour for reference_line
Colour_3 <- "#FD8D3C" #sets colour for cell_line_3
Colour_4 <- "#78C679" #sets colour for cell_line_4


#set name of how cell lines should be presented in the plot - ensures name consistency with other papers
#if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
Name_1 <- bquote("Methylation_Removed") #sets name for new cell line
Name_2 <- bquote("Normal_Cell")                   #sets name for reference_line
#Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
Name_4 <- bquote("Delocalised_Mutant")    #sets name for cell_line_4

#set name of plot titles
Title <- "Volume of Nucleus Occupied byXist RNA"
#x_axis <- "Cell Type"
y_axis <- bquote("Volume"~"of"~"Xist"~"RNA"~"["*mu*"m"^"3"*"]")


#################################################################################
#STEP 1: load Cloud_Volume data for plotting

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data is stored")

#open Cloud_Volume data for all existing cell lines
if (file.exists(paste(File_Path, Cell_Type, "Cloud_Volume", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Cloud_Volume_Compile.csv", sep="/"))) { 
  All_Cloud_Volume_Data <- read_csv(paste(File_Path, Cell_Type, "Cloud_Volume", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Cloud_Volume_Compile.csv", sep="/"))
} else{
  stop("New_All_Cell_Lines_Cloud_Volume_Compile.csv file does not exist in directory
       - check File_Path input")
}


#checks if new cell line data is stored in the dataframe containing all Cloud_Volume data
if (any(All_Cloud_Volume_Data$Cell_Line == New_Line_Name)) {
  print("New cell line data is found within Main Dataframe")
} else {
  stop("This new cell line does not exist in Main Dataframe (All_Cloud_Volume_Data)
       - make sure to run Cloud_Volume_Compilation_Manipulation.R ")
}

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of box plots

#chose which data to present in box plot
Present_Data <- All_Cloud_Volume_Data %>%
  filter(Cell_Line %in%
           c(Refernece_Line, New_Line_Name, Cell_Line_4)) 

#order data to be presented in box plot 
Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Refernece_Line, New_Line_Name, Cell_Line_4), ordered = TRUE,
                                 labels = c(Name_2, Name_1, Name_4))  


#################################################################################
#STEP 3: PLOT THE DATA

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (20)),
               legend.title = element_text(family = "Calibri", size = (16)), 
               legend.text = element_text(family = "Calibri", size = (14)), 
               axis.title = element_text(family = "Calibri", size = (16)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (12)),
               strip.text.x = element_text(family = "Calibri", size = (16)))

#sets name of labels for each cell line
Labels = c(Cell_Line = Name_1,
           Cell_Line = Name_2,
           #Cell_Line = Name_3,
           Cell_Line = Name_4)

#calculate the y scale for the plot
Min_y = min(Present_Data$Cloud_Volume, na.rm = TRUE) 
Max_y = max(Present_Data$Cloud_Volume, na.rm = TRUE)  

#make box plot
box_plot <- ggplot(Present_Data, aes(x = Cell_Line, y = Cloud_Volume, fill = Cell_Line)) + #defines data to present
  geom_boxplot(outlier.shape = NA) +                                            #defines type of plot
  #facet_wrap(~Cell_Line, nrow = 1, labeller = label_parsed) +                   #separates cell lines, and sets labels
  scale_fill_manual(values = c(Colour_1, Colour_2, Colour_4)) +                                                     
  coord_cartesian(ylim = c((Min_y-10),(1200) )) +                               #adjusts scale so whiskers don't touch the end graph 
  labs(title = Title,                                                           #sets name of the axes
      y = y_axis) +                                          
  theme                                                                         #adds the theme which defines the font, text size etc

#view the box plot
box_plot


################################################################################
#STEP 4: need to manually export the plots

################################################################################
#STEP 5: preform Wilcox Test on Data to asses if there is a statistical difference

#is there a general difference between the data sets?
Stat_Diff <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Both_Phases", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < p) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- paste("Overall, the new cell line has statistically bigger clouds than the", Refernece_Line, sep=" ")
  } else {
    outcome_1 <- paste("Overall, the new cell line has statistically smaller clouds than the", Refernece_Line, sep=" ")
  }} else {
    outcome_1 <- paste("Overall, there is no statistical difference in cloud size between the new cell line and the", Refernece_Line, sep=" ")
  }

#is there a general difference during the expansion phase between the data sets?
Expansion_Data <- filter(All_Cloud_Volume_Data, Phase == "Initiation")

Exp_Stat_Diff <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                             Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
Exp_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Exp_Diff_Less <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                             Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
Exp_Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Exp_Diff_More <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                             Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
Exp_Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
e <- Exp_Stat_Diff$p.value
f <- Exp_Diff_Less$p.value
g <- Exp_Diff_More$p.value
Exp_p_values <- tibble(e, f, g) %>%
  transmute(Comparision = paste("Compare_Exp_Phase", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = e,
            Diff_Statistically_Less = f,
            Diff_Statistically_More = g) 

if (Exp_Stat_Diff$p.value < p) {
  if (Exp_Diff_Less$p.value < Exp_Diff_More$p.value){
    outcome_2 <- paste("The new cell line has statistically bigger clouds, during the expansion phase, than the", Refernece_Line, sep=" ")
  } else {
    outcome_2 <- paste("The new cell line has statistically smaller clouds, during the expansion phase, than the", Refernece_Line, sep=" ")
  }} else {
    outcome_2 <- paste("There is no statistical difference in cloud size, during the expansion phase, between the new cell line and the", Refernece_Line, sep=" ")
  }

#is there a general difference during the Steady_State phase between the data sets?
Steady_State_Data <- filter(All_Cloud_Volume_Data, Phase == "Maintenance")

SS_Stat_Diff <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                            Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
SS_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
SS_Diff_Less <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                            Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
SS_Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
SS_Diff_More <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Refernece_Line),]$Cloud_Volume,
                            Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
SS_Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
h <- SS_Stat_Diff$p.value
i <- SS_Diff_Less$p.value
j <- SS_Diff_More$p.value
SS_p_values <- tibble(h, i, j) %>%
  transmute(Comparision = paste("Compare_SS_Phase", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = h,
            Diff_Statistically_Less = i,
            Diff_Statistically_More = j) 

if (SS_Stat_Diff$p.value < p) {
  if (SS_Diff_Less$p.value < SS_Diff_More$p.value){
    outcome_3 <- paste("The new cell line has statistically bigger clouds, during the Steady_State phase, than the", Refernece_Line, sep=" ")
  } else {
    outcome_3 <- paste("The new cell line has statistically smaller clouds, during the Steady_State phase, than the", Refernece_Line, sep=" ")
  }} else {
    outcome_3 <- paste("There is no statistical difference in cloud size, during the Steady_State phase, between the new cell line and the", Refernece_Line, sep=" ")
  }

#combine all p values
p_values <- p_values %>%
  bind_rows(Exp_p_values,
            SS_p_values)

#save p values
Save_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", New_Line_Name, sep="/")

write_csv(p_values, paste(Save_Path, paste(New_Line_Name, "_Cloud_Volume_Wilcox_p_values.csv", sep=""), sep="/"))


#provide outcomes of Wilcoxon tests
if (Exp_Stat_Diff$p.value < p) {
  print(outcome_1)
  print(outcome_2)
  print(outcome_3)
} else {
  print(outcome_1)
}

