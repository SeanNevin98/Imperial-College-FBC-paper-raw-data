#PROLOG
#Purpose: Plot Figure 1B, Figure 2, and Figure S1, Figure S2, Figure S3, Figure S4.
#Edited by: Sean Nevin

rstudioapi::restartSession()
#clear environment
rm(list = ls())

#load necessary packages
library(tidyverse)
library(readxl)
library(readr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr) #stat_qq_line - qq plot code in this package. 
library(purrr)
library(ggforce)
library(ggview)
library(ggridges) #for ridge plots - overlaid density plots. 

setwd() #set wd to directory containing raw data file. 

# Read in and raw data
raw_data <- read_excel("GitHub Repository - FBC paper raw data.xlsx", sheet = "raw_data")


# Wrangle data - select variables of interest, convert to long format, correct data types
vk_data <- raw_data %>% 
  mutate_all(~ifelse(. == ".", NA, .)) %>% # Replace "." with NA
  dplyr::select(ID_no, vload_PP:vload_Conv) %>%
  melt(id.vars = "ID_no", variable.name = "variable") %>%
  mutate(
    value = as.numeric(value),
    log10_value = log10(value + 1)
  ) %>%
  separate(variable, into = c("variable", "inf_time"), sep = "_", extra = "merge") %>%
  rename(Sample = ID_no)

# Set infection timepoint order
timings <- c("PP", "FP", "FP+4", "FP+7", "FP+14", "Conv")
vk_data$inf_time <- factor(vk_data$inf_time, levels = timings)

str(vk_data)


# PLOTTING SECTION --------------------------------------------------------
# Figure 1B ---------------------------------------------------------------


(Fig1B <- ggplot(data = vk_data, aes(x = inf_time, y = log10_value)) +
  geom_point(
    aes(group = Sample),
    color = "grey40",
    size = 1.5,
    shape = 19
  ) +
  geom_hline(
    yintercept = log10(8),
    colour = "sienna1",
    linetype = 2,
    linewidth = 1
  ) +
  geom_line(
    data = vk_data[!is.na(vk_data$log10_value), ],
    aes(group = Sample),
    linewidth = 0.75,
    alpha = 0.75,
    colour = "grey70"
  ) +
  xlab("Timepoint") +
  ylab("Viral Copy Number (Log10 RNA copies/mL)") +
  scale_x_discrete(limits = timings) +
  theme_classic(base_size = 22) +
  theme(panel.border = element_rect(
    color = "black",
    fill = NA,
    linewidth = 1
  ))
)

# Figure 2 and S1 ---------------------------------------------------------

#Wrangle FBC data
# Create the data set of significant comparisons
# All timepoint means were compared using a multiple mixed-effects model,
# And Dunnet's multiple comparison test in GraphPad PRISM.
# Multiple parameters were corrected for using Bonferroni adjustment.

# Reload in raw data
raw_data <- read_excel("GitHub Repository - FBC paper raw data.xlsx")

# Load in table of statistical comparison test results - analysed in GraphPad
signif_comp <- read_csv("Statistical comparisons Dunnets test with Bonferroni correction.csv")
#Remove NA rows
signif_comp <- signif_comp %>% drop_na(p.signif)

# Create a vector of the timings (x axis groups)
timings <- c("PP","FP","FP+7","FP+14","Conv")  

# Melt data into long format
df_melt <- raw_data %>%
  mutate_all(~ifelse(. == ".", NA, .)) %>% #potentially keep if lines do not pass through NA values
  reshape2::melt(id.vars = c("ID_no")) %>% 
  drop_na(value)

# Function to split string in variable column
separate_and_create_new_cols <- function(data, column_name, separator, left_col_name, right_col_name) {
  data <- separate(data, column_name, into = c(left_col_name, right_col_name), sep = separator, extra = "merge")
  return(data)
}

# Split variable column into variable and infection timepoint
data <- separate_and_create_new_cols(df_melt, "variable", "[_]", "variable", "inf_time")

# Rename column names
colnames(data) = c("Sample", "variable", "inf_time", "value")

# Convert variable names into non-abbreviated terms
data$variable[data$variable == "WBC"]         <- "White Blood Cells"
data$variable[data$variable == "Neutrophils"] <- "Neutrophils"
data$variable[data$variable == "Lymphocytes"] <- "Lymphocytes"
data$variable[data$variable == "Monocytes"]   <- "Monocytes"
data$variable[data$variable == "Eosinophils"] <- "Eosinophils"
data$variable[data$variable == "Basophils"]   <- "Basophils"
data$variable[data$variable == "Platelets"]   <- "Platelets"
data$variable[data$variable == "MPV"]         <- "Mean Platelet Volume"
data$variable[data$variable == "RBC"]         <- "Red Blood Cells"
data$variable[data$variable == "Haem"]        <- "Haemoglobin"
data$variable[data$variable == "Haematocrit"] <- "Haematocrit"  # Not used in the plots, but keeping as is
data$variable[data$variable == "MCV"]         <- "Mean Corpuscular Volume"
data$variable[data$variable == "MCH"]         <- "Mean Corpuscular Haemoglobin"
data$variable[data$variable == "MCHC"]        <- "MCH Concentration"

# Convert value to numeric class
data$value <- as.numeric(data$value)

str(data)
table(data$variable)



# PLOT VIOLIN DATA DISTRIBUTIONS AT EACH TIMEPOINT
# The rectangle scale factor (rect_sf) default is 20,
# Decreasing the value increases the distance between the rectangles

#The asterisk scale factor (text_sf) default is 1,
# Increasing this value increases the distance from the rectangle to the asterisk above

#### Initialising the Pre-requisite functions ####
# This works like a seq() function
special_seq <- function(start, gap, length) {
  
  sequence <- c(start)
  temp <- start
  
  # Loop through sequence
  for (i in 1:length-1) {
    # add the gap length each time
    temp = temp + gap
    sequence <- c(sequence, (temp))
  }
  
  return(sequence[1:length])
  
}


# Creates significance symbols for adjusted p_values if required
add_signif <- function(p_adjusted) {
  
  sig_symbol_list <- c()
  
  for (i in 1:length(p_adjusted)) {
    
    if (p_adjusted[i] < 0.001) {
      sig_symbol_list <- c(sig_symbol_list, "***")
    } else if (p_adjusted[i] < 0.01) {
      sig_symbol_list <- c(sig_symbol_list, "**")
    } else if (p_adjusted[i] < 0.05) {
      sig_symbol_list <- c(sig_symbol_list, "*")
    } else {
      sig_symbol_list <- c(sig_symbol_list, "")
    }
    
    
  }
  return(sig_symbol_list)
}


# This calculates the coordinates of the comparison lines on the final graph
rect_dataMaker <- function(df, all_data , cell, rsf) {
  
  # Isolate the cell you want to use
  df <- df[df$cell_type == cell,]
  # Identify the number of significant comparisons
  sig_no <- nrow(df)
  
  # Isolate data to plot
  all_data <- all_data[all_data$variable == cell,]
  
  
  ### Identifying the 'y' coordinates of the rectangles ###
  max_value <- as.numeric(summary(all_data$value)[6])
  
  # Calculating the relative scale
  
  scale_no <- ceiling(max_value/rsf)
  
  if (max_value < 10) {
    scale_no <- (max_value/rsf)
  }
  
  #range_val <- max(all_data$value, na.rm = TRUE) - min(all_data$value, na.rm = TRUE)
  range_val <- max(all_data$value, na.rm = TRUE) - sort(all_data$value, na.last = NA)[2]
  buffer <- 0.12 * range_val  # 12% of the range
  
  first_value <- max_value + scale_no + buffer
  
  #first_value <- (max_value + scale_no)
  ymin = special_seq(first_value, scale_no, sig_no)
  ymax = as.character(ymin + (scale_no*0.07))
  ymin = as.character(ymin)
  
  
  #### Identifying the 'x' coordinates for the rectangles ###
  #### 
  
  #change timepoints in signif comp to x axis values
  df <- df %>% 
    mutate(group1 = case_when(
      group1 == "Conv" ~ 5),
      group2 = case_when(
        group2 == "PP" ~ 1, 
        group2 == "FP" ~ 2,
        group2 == "FP+7" ~ 3,
        group2 == "FP+14" ~ 4
      ))
  
  df$group1 <- as.numeric(df$group1)
  df$group2 <- as.numeric(df$group2)
  
  xmin <- c()
  xmax <- c()
  
  for (i in 1:nrow(df)) {
    xmin <- c(xmin,min(df$group1[i],df$group2[i]))
    xmax <- c(xmax,max(df$group1[i],df$group2[i]))
  }
  
  # Final dataset of coordinates
  final_df <- data.frame(xmin=xmin, xmax=xmax,
                         ymin=ymin, ymax=ymax)
  
  # Makes sure all coordinates are numeric
  final_df$xmin <- as.numeric(final_df$xmin)
  final_df$xmax <- as.numeric(final_df$xmax)
  final_df$ymin <- as.numeric(final_df$ymin)
  final_df$ymax <- as.numeric(final_df$ymax)
  
  return(final_df)
  
}


# This calculates the coordinates of the significance symbols
text_dataMaker <- function(df,cell,rect_data,tsf) {
  
  # Isolate the cell you want to use
  df <- df[df$cell_type == cell,]
  sig_no <- nrow(df)
  
  x <- c()
  y <- c()
  lab <- c()
  
  for (i in 1:nrow(df)) {
    
    x_value <- (((rect_data$xmax[i] - rect_data$xmin[i])/2) + rect_data$xmin[i])
    y_value <- rect_data$ymin[i] + (tsf*0.3)
    
    if (min(rect_data$ymin) < 5 & min(rect_data$ymin) > 1) {
      y_value <- rect_data$ymin[i] + (tsf*0.03)
    } else if ((min(rect_data$ymin) < 5)) {
      y_value <- rect_data$ymin[i] + (tsf*0.003)
    }
    
    lab_value <- df$p.signif[i]
    
    x <- c(x,x_value)
    y <- c(y,y_value)
    lab <- c(lab,lab_value)
  }
  
  final_data <- data.frame(x=x, y=y, lab=lab)
  return(final_data)
  
  
}



harry_plotter <- function(cell_query, rect_sf=20, text_sf=1) {
  
  # Print the cell type being generated currently...
  print(paste("Making ", cell_query, " plot..."))
  
  
  # Plot for graph that has significant correlations
  if (cell_query %in% signif_comp$cell_type) {
    
    # Generate the data for the comparison line coordinates
    rect_data <- rect_dataMaker(signif_comp,data,cell_query,rect_sf)
    # Generate the data for the significance symbol coordinates
    text_data <- text_dataMaker(signif_comp,cell_query,rect_data,text_sf)
    
    test_plot <- ggplot(data[data$variable == cell_query,], aes(inf_time, value)) + 
      facet_grid(. ~ variable) +
      geom_line(aes(group = Sample), color="gray70") +
      geom_point(aes(group = Sample), color="gray40") +
      geom_violin(trim = FALSE, colour = "black", fill = NA, scale = "area") +
      theme_pubr() +
      ylim(0,NA) +
      xlab("")+
      scale_x_discrete(limits=timings) +
      geom_rect(data=rect_data, inherit.aes=FALSE,
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
      geom_text(data=text_data, inherit.aes = FALSE, size=6,
                aes(label=lab,x=x,y=y)) +
      stat_summary(fun.data=mean_cl_normal, geom='smooth', aes(group=1), 
                   fun.args = list(conf.int = .95), #can change confidence intervals here
                   data = data[data$variable == cell_query,],
                   se=T, color='steelblue4',
                   fill='steelblue4')+
      theme_classic(base_size = 16) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.text.x = element_text(face = "bold"))
    
    
    
  } else {
    
    # Plot for graph that desn't have significant correlations
    test_plot <- ggplot(data[data$variable == cell_query,], aes(inf_time, value)) + 
      facet_grid(. ~ variable) +
      geom_line(aes(group = Sample), color="gray70") +
      geom_point(aes(group = Sample), color="gray40") +
      geom_violin(trim = FALSE, colour = "black", fill = NA, scale = "area") +
      theme_pubr() +
      ylim(0,NA) +
      xlab("")+
      scale_x_discrete(limits=timings) +
      stat_summary(fun.data=mean_cl_normal, geom='smooth', aes(group=1), 
                   fun.args = list(conf.int = .95), #can change confidence intervals here
                   data = data[data$variable == cell_query,],
                   se=T, color='steelblue4',
                   fill='steelblue4')+
      theme_classic(base_size = 16) +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.text.x = element_text(face = "bold"))
    
    
  }
  return(test_plot)
  
}




# Figure 2 ----------------------------------------------------------------

# Plot the FBC data as a ggarrange object
Fig2_violin <- ggpubr::ggarrange(
  
  harry_plotter("White Blood Cells", text_sf = 0.5, rect_sf = 110) +
    geom_hline(yintercept = 4.2, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 11.2, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(1.8, NA) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Neutrophils", text_sf = 0.5, rect_sf = 15) +
    geom_hline(yintercept = 2, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 7.1, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Lymphocytes", text_sf = 0.2) +
    geom_hline(yintercept = 1.1, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 3.6, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(0.4, NA) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Monocytes", text_sf = 0.4) +
    geom_hline(yintercept = 0.3, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 0.9, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Platelets", rect_sf = 15, text_sf = 25) +
    geom_hline(yintercept = 130, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 400, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(60, NA) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Mean Platelet Volume", text_sf = 0.5, rect_sf = 80) +
    geom_hline(yintercept = 7.5, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 11.5, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(5, NA) +
    labs(y = "fL"),
  
  nrow = 2,
  ncol = 3
)

print(Fig2_violin)



# Figure S1
FigS1_violin <- ggpubr::ggarrange(
  
  harry_plotter("Eosinophils", rect_sf = 20, text_sf = 5) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Basophils", rect_sf = 10) +
    labs(y = expression(10^9/L)),
  
  harry_plotter("Red Blood Cells", rect_sf = 30, text_sf = 0.1) +
    geom_hline(yintercept = 3.73, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 5.46, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(3.5, NA) +
    labs(y = expression(10^12/L)),
  
  harry_plotter("Haemoglobin", rect_sf = 25, text_sf = 4) +
    geom_hline(yintercept = 114, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 168, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(70, NA) +
    labs(y = expression(g/L)),
  
  harry_plotter("Mean Corpuscular Volume", rect_sf = 100, text_sf = 2) +
    geom_hline(yintercept = 83.5, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 99.5, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(80, NA) +
    labs(y = expression(fL)),
  
  harry_plotter("Mean Corpuscular Haemoglobin", rect_sf = 50, text_sf = 0.8) +
    geom_hline(yintercept = 27.5, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 33.1, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(19.8, NA) +
    labs(y = expression(pg)),
  
  harry_plotter("MCH Concentration", text_sf = 5, rect_sf = 50) +
    geom_hline(yintercept = 315, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = 350, colour = "sienna1", linetype = 2, linewidth = 0.5) +
    ylim(310, NA) +
    labs(y = expression(g/L)),
  
  ncol = 3,
  nrow = 3
)

print(FigS1_violin)





# Figure S3 ---------------------------------------------------------------
# Visualising data distributions using historgram and QQ plots

# Check normal distributions of each FBC parameter - 
# Perform a histogram, with density overlaid
# Perform a qq plot

# Define plotting order of variables
var_list <- c(
  "White Blood Cells",
  "Neutrophils",
  "Lymphocytes",
  "Monocytes",
  "Platelets",
  "Mean Platelet Volume",
  "Eosinophils",
  "Basophils",
  "Red Blood Cells",
  "Haemoglobin",
  "Mean Corpuscular Volume",
  "Mean Corpuscular Haemoglobin",
  "MCH Concentration"
)

data <- filter(data, variable != "vload")

# Initiate empty lists to collect plots
hist_list <- list()
qq_list   <- list()

for (i in seq_along(var_list)) {
  
  var <- var_list[i]
  
  # Filter dataframe for selected variable
  data_sub <- data %>%
    filter(variable == var)
  
  # Histogram + density
  p_hist <- ggplot(data_sub, aes(x = value)) +
    geom_histogram(aes(y = ..density..),
                   bins = 30,
                   fill = "grey",
                   colour = "black",
                   alpha = 0.6) +
    geom_density(colour = "red", linewidth = 1.2) +
    labs(x = var,
         y = "Density",
         title = paste(var)) +
    theme_classic(base_size = 13.5) +
    theme(plot.title = element_text(face = "bold"))
  
  # QQ plot
  p_qq <- ggplot(data_sub, aes(sample = value)) +
    stat_qq(size = 1.2, colour = "blue") +
    stat_qq_line(colour = "red", linewidth = 1) +
    labs(title = paste(var)) +
    theme_classic(base_size = 13.5)  +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_blank())
  
  # Show in loop
  print(p_hist)
  print(p_qq)
  
  # Store plots in lists
  hist_list[[i]] <- p_hist
  qq_list[[i]]   <- p_qq
}

# Arrange plots
# Figure S3 - histograms of FBC parameters
(FigureS3 <- ggarrange(plotlist = hist_list,
                                   nrow = 5, ncol = 3))

#Not shown in final manuscript - QQ plots of parameters
(individual_qqplots <- ggarrange(plotlist = qq_list,
                                nrow = 5, ncol = 3))
  
(qqplots_mod <- ggpubr::annotate_figure(
    individual_qqplots,
    left = ggpubr::text_grob("Sample Quantiles", rot = 90),
    bottom = ggpubr::text_grob("Theoretical Quantiles")
  )
)


# Figure S4 ---------------------------------------------------------------

#Use ggridge for easier comparison of data distributions. 
# Optional: order timepoints if needed
timepoints <- c("FP","FP+7","FP+14","Conv")
data$inf_time <- factor(data$inf_time, levels = timepoints)
data$variable <- factor(data$variable, levels = var_list)


#filter dataframe to var in varlist
  data_sub <- data %>% 
    filter(variable != "vload") %>% 
    filter(inf_time != "PP")
  
unique(data_sub$variable)
unique(data$inf_time)
  
#use facet_grid to plot each inf timepoint.
(FigureS4 <- density_plot_comparison <- 
  ggplot(data_sub, aes(x = value, y = inf_time)) +
    geom_density_ridges(
      alpha = 0.7,
      quantile_lines = TRUE,
      quantiles = 2) +
    facet_wrap(~variable, 
               scales = "free", 
               nrow = 5, 
               ncol = 3) +  # one row per variable, one column per timepoint
    labs(x = "FBC measurement value",
         y = "Density") + 
    theme_classic(base_size = 13.5) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          strip.text.x = element_text(face = "bold"))
  ) 





# Figure S2 ---------------------------------------------------------------

#Ppts with the greatest decrease in Mean Platelet volume at FP+7 see
#the greatest increase in platelet count at FP+14
raw_data <- read_excel("GitHub Repository - FBC paper raw data.xlsx")

df <- raw_data %>% 
  mutate_all(~ifelse(. == ".", NA, .)) %>% 
  mutate_at(c(1:77), as.numeric)

my_theme <- theme_classic(base_size = 12) +
  theme(panel.border = element_rect(colour = "black",
                                    fill = NA, 
                                    linewidth = 1))

(FigureS2 <- ggplot(data = df, aes(x = `MPV_FP+7`, y = `Platelets_FP+14`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..rr.label.., ..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson", size = 6) +
  my_theme +
  labs(x = "Mean Platelet Volume (MPV) at FP+7", 
       y = "Platelet count at FP+14")
)



# Run linear regression
model <- lm(`Platelets_FP+14` ~ `MPV_FP+7`, data = df)
summary(model)

# Regression coefficients with SE and p-values
coefs <- summary(model)$coefficients

# Confidence intervals
confint_vals <- confint(model)

# R-squared and adjusted R-squared
r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

# Combine into one data frame
results <- cbind(
  Estimate = coefs[, "Estimate"],
  Std_Error = coefs[, "Std. Error"],
  p_value = coefs[, "Pr(>|t|)"],
  CI_lower = confint_vals[, 1],
  CI_upper = confint_vals[, 2]
)

print(round(results, 4))
cat("R-squared:", round(r_squared, 4), "\n")
cat("Adjusted R-squared:", round(adj_r_squared, 4), "\n")


# END OF SCRIPT -----------------------------------------------------------


