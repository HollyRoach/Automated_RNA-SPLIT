################################################################################
#MANUAL VS AUTOMATIC - Molecule Count and Transcription Dynamics Validation
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

setwd(paste(File_Path, Cell_Type, "Turnover", "All_Cell_Lines_Merged", sep="\\"))

#loads manual data
WT_Manual <- read_csv("Original_WT_Centroids_Compile_Includes_Turnover_Pulse_2.csv")

#loads automated data
WT_Auto <- read_csv("Updated_WT_Auto_Centroid_Compile.csv") 
  

################################################################################
#TOTAL MOLECULE cOUNT ANALSYIS 
#calculate total number of centroids per cloud - WT manual compile
#separate into pulse 1 and pulse 2 data

#calculate total molecules for dynamic manual WT
Manual_D_Pulse_1 <- WT_Manual %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Manual_D_Pulse_2 <- WT_Manual %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_2") %>%
  rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Manual_D_Total_Count <- Manual_D_Pulse_1 %>%
  mutate(P2_No_Centroids = Manual_D_Pulse_2$P2_No_Centroids,
         Total_Centroids = Manual_D_Pulse_1$P1_No_Centroids + Manual_D_Pulse_2$P2_No_Centroids) %>%
  select(-Normalised_No_Centroids)

#calculate total molecules for Turnover manual WT
Manual_T_Pulse_1 <- WT_Manual %>%
  filter(Data_Set == "Turnover" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Manual_T_Pulse_2 <- WT_Manual %>%
   filter(Data_Set == "Turnover" & Pulse == "Pulse_2") %>%
   rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Manual_T_Total_Count <- Manual_T_Pulse_1 %>%
   mutate(P2_No_Centroids = Manual_T_Pulse_2$P2_No_Centroids,
          Total_Centroids = Manual_T_Pulse_1$P1_No_Centroids + Manual_T_Pulse_2$P2_No_Centroids) %>%
   select(-Normalised_No_Centroids)

#combine dynamic and turnover data for manual compile
Manual_Total_Count <- Manual_D_Total_Count %>%
  bind_rows(Manual_T_Total_Count) %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_manual"))

################################################################################
#STEP b: TOTAL MOLECULE COUNT ANALSYIS 
##calculate total number of centroids per cloud - WT auto compile
#separate into pulse 1 and pulse 2 data

#calculate total molecules for dynamic Auto WT
Auto_D_Pulse_1 <- WT_Auto %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Auto_D_Pulse_2 <- WT_Auto %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_2") %>%
  rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Auto_D_Total_Count <- Auto_D_Pulse_1 %>%
  mutate(P2_No_Centroids = Auto_D_Pulse_2$P2_No_Centroids,
         Total_Centroids = Auto_D_Pulse_1$P1_No_Centroids + Auto_D_Pulse_2$P2_No_Centroids) %>%
  select(-Normalised_No_Centroids)

#calculate total molecules for Turnover Auto WT
Auto_T_Pulse_1 <- WT_Auto %>%
  filter(Data_Set == "Turnover" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

#cannot calculate as Auto compile does not contain pulse 2 turnover data 
Auto_T_Pulse_2 <- WT_Auto %>%
   filter(Data_Set == "Turnover" & Pulse == "Pulse_2") %>%
   rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
 Auto_T_Total_Count <- Auto_T_Pulse_1 %>%
   mutate(P2_No_Centroids = Auto_T_Pulse_2$P2_No_Centroids,
          Total_Centroids = Auto_T_Pulse_1$P1_No_Centroids + Auto_T_Pulse_2$P2_No_Centroids) %>%
   select(-Normalised_No_Centroids)

#combine dynamic and turnover data for Auto compile
Auto_Total_Count <- Auto_D_Total_Count %>%
  bind_rows(Auto_T_Total_Count)


################################################################################
#STEP c: TOTAL MOLECULE COUNT ANALSYIS 

# #plot original vs manual vs auto compile
# setwd(paste(File_Path, Cell_Type, "Molecule_Count", "All_Cell_Lines_Merged", sep="\\"))
# Original_Total_Count <- read_csv("Cleaned_Other_Cell_Lines_Molecule_Count_Compile.csv") %>%
#   filter(Cell_Line == "WT") %>%
#   mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_original"))
  
  
#combine auto and manual data to present
Present_Count <- Auto_Total_Count %>%
  select(-Cloud_Name) %>%
  bind_rows(Manual_Total_Count) %>% 
  select(Cell_Line, Phase, Total_Centroids) %>%
  #bind_rows(Original_Total_Count) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",                          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%                 
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                                         #creates key variable to base colour coding off


#order data to be presented in box plot 
Present_Count$Cell_Line <- factor(Present_Count$Cell_Line, 
                                 levels = c(#"WT_original",
                                            "WT_manual", "WT_auto"), ordered = TRUE) #use names defined in line 51

#put in order for the graphs to appear on the x-axis
Present_Count$Key <- factor(Present_Count$Key, 
                           levels = c(#"Expansion_WT_original", "Steady_State_WT_original",
                                      "Expansion_WT_manual", "Steady_State_WT_manual",
                                      "Expansion_WT_auto", "Steady_State_WT_auto"), ordered = TRUE)


#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri", size = (11)), 
               legend.text = element_text(family = "Calibri", size = (11)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (9)))


#make box plot
count_plot <- ggplot(Present_Count, aes(x = Phase, y = Total_Centroids, fill = Key)) +       #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1) +                                                      #shows box plots for selected cell lines
  scale_fill_manual(values = c(#"#FA9FB5", "#DD3497",                                      #WT_original exp and SS colours 
                               "#FC9272", "#EF3B2C",                                      #WT_manual exp and SS colours 
                               "#BDBDBD", "#737373" )) +                                  #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((0),(290) )) +                                             #adjusts scale so whiskers don't touch the end graph 
  labs(title = "A) Validate Total Molecule Count Analysis",                                                     #sets name of the axes
       x = "Phase",
       y = "Total Xist Count [no. centroids]") +                                         
  theme                                                                                   #adds the theme which defines the font, text size etc


#view the box plot
count_plot

################################################################################
#STEP d: TOTAL MOLECULE COUNT ANALSYIS 
#preform t-test Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("TOTAL_Compare_Both_Phases_WT_vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- "Overall, the automated compile has statistically higher centroid count than the manual compile"
  } else {
    outcome_1 <- "Overall, the new automated compile has statistically lower centroid count than the manual compile"
  }} else {
    outcome_1 <- "Overall, there is no statistical difference in centroid count between the automated compile and the manual compile"
  }


#save p values

setwd(paste(File_Path, Cell_Type, "Molecule_Count", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Count_t-test_p_values.csv", sep="_"))

#return outcome of the stats test
print(outcome_1)


################################################################################
#TRANSCRIPTION DYNAMICS ANALSYIS 
#clean Transcription Dynamics
#CHECK - I think Lisa has included 0min and 60min from turnover data - need to repeat when I get access to the turnover pulse 2

Manual_Transcription <- WT_Manual %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% 
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c("10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 10 ~ 20,                    #new data needs to grouped by times - 20, 40, 60
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

Auto_Transcription <- WT_Auto %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% 
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c( "10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 10 ~ 20,                   #new data needs to grouped by times - 20, 40, 60
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

# #load original transcription dynamics data
# setwd(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines_Merged", sep="\\"))
# Original_Transcription <- read_csv("Cleaned_Other_Cell_Lines_Transcription_Dynamics_Compile.csv") %>%
#   mutate(Cell_Line = "WT_original") %>%
#   filter(Time %in%                                           #selects only time 20, 40 and 60min to present
#            c( "0", "10", "20", "30", "40", "50", "60"))

#load automated data before specific datapoints are removed
setwd(paste(File_Path, Cell_Type, "Turnover", New_Line_Name, sep="\\"))

Auto_Original <- read_csv("WT_auto_Centroids_Compile_ALL_DATA.csv") %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% #only want to include nascent_Xist_dynamics data
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c( #"0",
              "10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 0 ~ 20,                    #new data needs to grouped by times - 20, 40, 60
                          Time == 10 ~ 20,       
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

################################################################################
#STEP b: TRANSCRIPTION DYNAMICS ANALSYIS 
#plot manual vs auto compile

Present_Transcription  <- Manual_Transcription %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_manual"))    %>%
  bind_rows(Auto_Original,
            #Original_Transcription
            ) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",
                           Phase == "Maintenance" ~ "Steady_State")) %>%
  mutate(Key = paste(Phase, Cell_Line, sep="_")) %>%                                  #creates variable to base colour coding in plot
  filter(Time %in%                                                                    #selects only time 20, 40 and 60min to present
           c("20","40", "60"))

#order data to be presented in box plot 
Present_Transcription$Cell_Line <- factor(Present_Transcription$Cell_Line, 
                                    levels = c(#"WT_original",
                                               "WT_manual", "WT_auto"), ordered = TRUE) 

#put in order for the graphs to appear on the x-axis
Present_Transcription$Key <- factor(Present_Transcription$Key, 
                                levels = c(#"Expansion_WT_original", "Steady_State_WT_original",
                                           "Expansion_WT_manual", "Steady_State_WT_manual",
                                           "Expansion_WT_auto", "Steady_State_WT_auto"), ordered = TRUE)

Present_Transcription$Time  <- factor(Present_Transcription$Time,
                             levels = c("20","40", "60"), ordered = TRUE)

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri", size = (11)), 
               legend.text = element_text(family = "Calibri", size = (11)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (9)))

#calculate the y scale for the plot
Min_y = min(Present_Transcription$No_Centroids, na.rm = TRUE)
Max_y = max(Present_Transcription$No_Centroids, na.rm = TRUE)  #adjust value according to graph output as doesn't contain outliers

#make box plot 
transcription_plot <- ggplot(Present_Transcription, aes(x = Time, y = No_Centroids, fill = Key)) + #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                         #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1) +                                                #shows box plots seperated by phase
  scale_fill_manual(values = c(#"#FA9FB5", "#DD3497",                                      #WT_original exp and SS colours 
                               "#FC9272", "#EF3B2C",                                      #WT_manual exp and SS colours 
                               "#BDBDBD", "#737373" )) +                                  #WT_auto exp and SS colours 
  coord_cartesian(ylim = c(Min_y, 90)) +                                        #adjusts so whiskers don't touch the end graph 
  labs(title = "C) Validate Transcription Dynamics Analysis",                              #sets name of the axes
       x = "Time [min]",
       y = "Newly Synthesized Xist Count [no. centroids]") +           
  theme                                                                         #adds the theme which defines the font, text size etc

#view the box plot
transcription_plot

################################################################################
#c: TRANSCRIPTION DYNAMICS ANALSYIS 
#preform t-test Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff_T <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff_T $p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less_T  <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less_T $p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More_T  <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More_T $p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
d <- Stat_Diff_T $p.value
e <- Diff_Less_T $p.value
f <- Diff_More_T $p.value
p_values <- tibble(d, e, f) %>%
  transmute(Comparision = paste("TOTAL_Compare_Both_Phases", New_Line_Name, sep='_'),
            General_Stat_Diff = d,
            Diff_Statistically_Less = e,
            Diff_Statistically_More = f) 

if (Stat_Diff_T$p.value < 0.05) {
  if (Diff_Less_T$p.value < Diff_More_T$p.value){
    outcome_2 <- "Overall, the automated compile has statistically higher transcription rate than the manual compile"
  } else {
    outcome_2 <- "Overall, the new automated compile has statistically lower transcription rate than the manual compile"
  }} else {
    outcome_2 <- "Overall, there is no statistical difference in transcription rate between the automated compile and the manual compile"
  }


#save p values

setwd(paste(File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Transcription_Dynamics_t-test_p_values.csv", sep="_"))

#return outcome of the stats test
print(outcome_1)
print(outcome_2)
