rm(list=ls())

#prepare the dataset for control tubes
################################################################################
#size and population relationship analysis on earth mite
library(readxl)
setwd("/Users/user/Documents/MSCA/Bulb mite experiment/Smallegange & Deere 2014/Dataset for density dependent effect")

CountData=read_xlsx("MQ test_ Count and Size data.xlsx", sheet = "CountData")
colnames(CountData)<-c(colnames(CountData)[1:23],"noeggs")
SizeData=read_xlsx("MQ test_ Count and Size data.xlsx", sheet = "SizeData")
# get the data for control
CountData_C<-subset(CountData,treatment=="C")
SizeData_C<-subset(SizeData,treatment=="C")
CountData_C$tube<-as.factor(CountData_C$tube)
SizeData_C$tube<-as.factor(SizeData_C$tube)
###############################################################################

#mark the data before equlibrium and after equlibrium and trend of population dynamic
###############################################################################


#logistic regression of population dynamic
#-------------------------------------------------------------------------------
   #change T3 to T5, T8 for tube 5, 8
   library(nls2)
   # Define the logistic model function
   logistic_model <- function(day, K, N0, r) {
     K / (1 + ((K - N0) / N0) * exp(-r * day))
   }

   # Initial parameter guesses
   start_params <- list(K = max(CountData_CT8$nototal), N0 = min(CountData_CT8$nototal), r = 0.1)

   # Fit the logistic model
   fit <- nls2(nototal ~ logistic_model(day, K, N0, r),
               data = CountData_CT8,
               start = start_params,
               algorithm = "port")

   # View the fit
   summary(fit)

   # Add fitted values to the data frame
   CountData_CT8$fit <- predict(fit, newdata = CountData_CT8$nototal)

   # Plot the time series with the fitted logistic model
   ggplot(CountData_CT8, aes(x = day)) +
     geom_line(aes(y = nototal), color = "blue", size = 1) +
     geom_line(aes(y = fit), color = "red", linetype = "dashed") +
     labs(title = "Population Size with Logistic Model Fit of tube 8", x = "Day", y = "Population Size") +
     theme_minimal()
#-------------------------------------------------------------------------------
   #mark the euqlibrium status
   for (t in c(3,5,8)) {
     rowID<-which(CountData_C$tube==t)
     if (t==3) {
       day_eq<-75
     } 
     if (t==5) {
       day_eq<-61
     }
     else {
       day_eq<-71
     }
     
     for (i in 1:(length(rowID))) {
       ifelse(CountData_C$day[rowID[i]]<=day_eq, 
              CountData_C$Equ[rowID[i]]<-"Pre-equ", 
              CountData_C$Equ[rowID[i]]<-"Post-equ")  
     }
   }
#----------------------------------------------------------------------------   
#continued from line 52
#calculate lambda of the regressed logistic growth over time
   # Extract parameters
   params <- coef(fit)
   # Print parameters
   K <- params["K"]
   N0 <- params["N0"]
   r <- params["r"]
   # Time points for which you want to calculate lambda
   time_points <- 1:365  # Example time points
   
   # Calculate population sizes at different times
   population_sizes <- logistic_model(time_points, K, N0, r)
   
   # Calculate finite growth rate (lambda) for each time point
   lambda_values <- c(NA, population_sizes[-1] / population_sizes[-length(population_sizes)])
   
   # Remove NA values for lambda and corresponding population sizes
   valid_indices <- !is.na(lambda_values)
   time_points <- time_points[valid_indices]
   population_sizes <- population_sizes[valid_indices]
   lambda_values <- lambda_values[valid_indices]
   
   # Create a data frame for plotting
   plot_data <- data.frame(
     Time = time_points,
     Population_Size = population_sizes,
     Lambda = lambda_values
   )
   
   # Plot the relationship between Lambda and Population Size
   ggplot(plot_data, aes(x = time_points, y = Lambda)) +
     geom_point(color = "blue") +
     geom_line(color = "red") +
     labs(
       title = "Relationship Between Lambda and Population Size",
       x = "Population Size",
       y = "Finite Growth Rate (Lambda)"
     ) +
     theme_minimal()
   #----------------------------------------------------------------------------
   #Calculate the lambda of observed data in the experiment and assign it to CountDate_C
       for (tube in c(3,5,8)) { #t= 5， 8 for tube 5 and 8
       rowID<-which(CountData_C$tube==tube) 
       sub_pop <- data.frame(
       Day = CountData_C$day[rowID],  
       Population_Size = CountData_C$nototal[rowID] 
       )
     all_days <- seq(min(sub_pop$Day), max(sub_pop$Day), by = 1)  # Creates a sequence of continuous days
     # Perform linear interpolation using approx()
     interpolated_population <- approx(x = sub_pop$Day, y = sub_pop$Population_Size, xout = all_days)
     
     # Create a new data frame with interpolated values
     interpolated_data <- data.frame(
       day = interpolated_population$x,
       Population_Size = interpolated_population$y
     )
     
     # View the interpolated data
     #print(interpolated_data)
     # Calculate finite growth rate (lambda)
     lambda_values <- c(NA, interpolated_data$Population_Size[-1] / interpolated_data$Population_Size[-length(interpolated_data$Population_Size)])
     # Add lambda values to the interpolated data frame
     interpolated_data$Lambda <- lambda_values
     interpolated_data$tube<-as.factor(tube)
     CountData_C$tube<-as.numeric(CountData_C$tube)
     for (i in 1:dim(interpolated_data[1])) {
       day<-interpolated_data$day[i]
       tube<-interpolated_data$tube[i]
       row_IDD<- which(CountData_C$day==day & CountData_C$tube==tube)
       CountData_C$Lambda[row_IDD]<-interpolated_data$Lambda[i]
       }
       }
   CountData_C$tube<-as.factor(CountData_C$tube)
#--------------------------------------------------------------------------------  
#analyze relationshio between lambda and population size 
  # Remove rows with missing Lambda values
     data_clean <- CountData_C %>%
       filter(!is.na(Lambda))
   
     skew_normal_model <- function(m,n, nototal, xi, omega, alpha) {
       phi <- dnorm((nototal - xi) / omega)
       # Standard normal cumulative distribution function
       Phi <- pnorm(alpha * (nototal - xi) / omega)
       # Compute the Skew Normal PDF
       m* (2 / omega) * phi * Phi+n
     }
   
  plot(data_clean$nototal, data_clean$Lambda)   
   lines (1:800,skew_normal_model(70,1, 1:800,200,46.2,-0.05))
   
     # Initial parameter guesses
     start_params <- list(m=70 , n=1, xi=200, omega=46.2, alpha=-0.05)
     fit_skew_normal <- nls2(Lambda ~ skew_normal_model(m,n, nototal, xi, omega, alpha),
                       data = data_clean,
                       start = start_params,
                       algorithm = "port")
     
     # View the fit
     summary(fit_skew_normal)
     plot(predict(fit_skew_normal, 1:800))
     # Add fitted values to the data frame
     data_clean$fit <- predict(fit_skew_normal, newdata = data_clean$nototal)
     #!!!!!need to be tested
     data_clean$fit <- skew_normal_model(60,1, data_clean$nototal,200,46.2,-0.05)

     
     # Plot the time series with the fitted logistic model
     ggplot(data_clean, aes(x = nototal)) +
       geom_point(aes(y = Lambda), color = "blue", size = 1) +
       geom_line(aes(y = fit), color = "red", linetype = "dashed") +
       labs(title = "Lambda vs population size", x = "Population size", y = "Lambda") +
       theme_minimal()
     
################################################################################

#plot population structure dynamic across tubes
##################################################################
setwd("/Users/user/Documents/MSCA/Bulb mite experiment/Smallegange & Deere 2014/Dataset for density dependent effect/figures")
tiff("Population structure dynamics.tif", width = 4, height = 3, units = "in", res = 600)
# Load necessary libraries
library(tidyverse)
# Reshape the data to long format
long_data <- CountData_C %>%
  #pivot_longer(cols = c(noadults, nojuveniles, noeggs),  # Specify all stage group columns here
  pivot_longer(cols = c( noeggs, nojuveniles, noadults),  # Specify all stage group columns here
               names_to = "stage_group", 
               values_to = "count")
# Plot totalno vs day, with totalno composed of individuals from each stage group, and tube as a factor
long_data$stage_group<-factor(long_data$stage_group,levels=c("noeggs","nojuveniles","noadults"))
ggplot(long_data, aes(x = factor(day), y = count, fill = stage_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ tube) + 
  scale_x_discrete(breaks = unique(long_data$day)) +  # Show only 5 evenly spaced ticks
  labs(title = "Total Number of Individuals vs Day",
       x = "Day",
       y = "Total Number of Individuals",
       fill = "Stage Group") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),  # X-axis title size
        axis.title.y = element_text(size = 12),  # Y-axis title size
        axis.text.x  = element_text(size = 10,angle = 45, hjust = 1),  # X-axis text size
        axis.text.y  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))  # Y-axis text size)

dev.off()
###########################################################################

#plot the correlation between size of individuals from different stage groups with size of stage groups
#############################################################################

#merge countdata and size data
merged_data <- left_join(SizeData_C, CountData_C, by = c("day", "tube"))
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- cor(x, y, use="pairwise.complete.obs") # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * abs(Cor)) # Resize the text by level of correlation
}

merged_data_Preequ<-subset(merged_data, merged_data$Equ=="Pre-equ")
merged_data_Postequ<-subset(merged_data,merged_data$Equ=="Post-equ")

al <-merged_data[, c(9:20, 42, 43, 44, 41)] 
al <-merged_data_Preequ[, c(9:15, 42, 43, 44, 41)] 
al <-merged_data_Postequ[, c(9:20, 42, 43, 44, 41)] 

pairs(al,
      #col=merged_data$tube, 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth, na.rm=T) # Smoothed regression lines)


al <-merged_data[, c(42,9)]  #F vs noadults
pairs(al,
      #col=merged_data$trend, 
      col=merged_data$Equ, 
      #col=merged_data$tube, 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth, na.rm=T) # Smoothed regression lines)



al_log <- data.frame(noadults = log(al$noadults), F = log(al$F))
pairs(
  al_log,
  col=merged_data$Equ, 
  #col=merged_data$Equ, 
  #col=merged_data$tube, 
  upper.panel = panel.cor,
  # Correlation panel
  lower.panel = panel.smooth,
  na.rm = T
) # Smoothed regression lines)


############################################################################3

# relationshiop at specific tube
#############################################################################
#Tube 3
merged_data_T3<-subset(merged_data,tube==3)
al <-merged_data_T3[, c(9:20, 42, 43, 44, 41,46)] 

al<-merged_data_T3[, c(42, 9)]  #F vs noadults
al<-merged_data_T3[, c(42, 10)] # FM vs noadults
al<-merged_data_T3[, c(42, 11)] # SM vs noaudlts
al<-merged_data_T3[, c(42, 16)] # Qtritos vs noaudlts
al<-merged_data_T3[, c(42, 17)] # Qprotos vs noaudlts
al<-merged_data_T3[, c(42, 19)] # deutos vs noaudlts

#select one al from above
pairs(al,
      col=merged_data_T3$Equ, 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth, na.rm=T) # Smoothed regression lines)


#Tube 5
merged_data_T5<-subset(merged_data,tube==5)
al <-merged_data_T5[, c(9:20, 42, 43, 44, 41)] 
# Remove columns that are completely NA
al <- al[, colSums(is.na(al)) < nrow(al)] 

al<-merged_data_T5[, c(42, 9)]  #F vs noadults
al<-merged_data_T5[, c(42, 10)] # FM vs noadults
al<-merged_data_T5[, c(42, 11)] # SM vs noaudlts
al<-merged_data_T5[, c(42, 16)] # Qtritos vs noaudlts
al<-merged_data_T5[, c(42, 17)] # Qprotos vs noaudlts
al<-merged_data_T5[, c(42, 19)] # deutos vs noaudlts

#select one al from above
pairs(al,
      col=merged_data_T5$Equ, 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth, na.rm=T) # Smoothed regression lines)

#Tube 8
merged_data_T8<-subset(merged_data,tube==8)
al <-merged_data_T8[, c(9:20, 42, 43, 44, 41)] 
al<-merged_data_T8[, c(42, 9)]  #F vs noadults
al<-merged_data_T8[, c(42, 10)] # FM vs noadults
al<-merged_data_T8[, c(42, 11)] # SM vs noaudlts
al<-merged_data_T8[, c(42, 16)] # Qtritos vs noaudlts
al<-merged_data_T8[, c(42, 17)] # Qprotos vs noaudlts
al<-merged_data_T8[, c(42, 19)] # deutos vs noaudlts

#select one al from above
pairs(al,
      col=merged_data_T8$Equ, 
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth, na.rm=T) # Smoothed regression lines)

################################################################################

# Calculate mean and standard deviation for each combination of day and noadults
################################################################################
#For F vs noadults with st-dev as legend
#--------------------------------------------------------------------------------
#for tube 3,5, 8,  change merged_data_T3, to ..._5, ..._8 
al<-merged_data_T8[, c(4, 42, 9)]  #F vs noadults
summary_al <- al %>%
  group_by(day, noadults) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = noadults, y = mean_F, color = sd_F)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)
  
  # Customizing the x-axis scale (setting limits and breaks)
  scale_x_continuous(
    name = "Number of Adults (noadults)",  # Custom x-axis label
    limits = c(0, 80),  # Set limits from 0 to 20
    breaks = seq(0, 80, by = 10)  # Breaks every 5 units
  ) +
  
  # Customizing the y-axis scale (setting limits and breaks)
  scale_y_continuous(
    name = "Mean of F (mean_F)",  # Custom y-axis label
    limits = c(500, 900),  # Set limits from 0 to 100
    breaks = seq(500, 900, by = 100)  # Breaks every 20 units
  ) +
  
  # Customizing the color legend (scale for sd_F)
  scale_color_gradient(
    low = "blue", high = "red",  # Color gradient from blue to red
    limits = c(0, 200),  # Set limits for sd_F from 0 to 30
    breaks = seq(0, 200, by = 20),  # Breaks for the legend at intervals of 10
    name = "Standard Deviation (sd_F)"  # Legend title
  ) +
  
  labs(title = "Scatter Plot of mean_F vs noadults with Custom Scales") +
  theme_minimal()
#-----------------------------------------------------------------------------
##For F vs noadults with equalibrium status color coded
#-----------------------------------------------------------------------------
#for tube 3,5, 8,  change merged_data_T3, to ..._5, ..._8 
al<-merged_data_T8[, c(4, 42,46,9)]  #F vs noadults
summary_al <- al %>%
  group_by(day, noadults, Equ) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = noadults, y = mean_F, color = Equ)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs noadults with Custom Scales") +
  theme_minimal()
##-------------------------------------------------------------------------------
#For F vs nototal with st-dev as legend
#--------------------------------------------------------------------------------
#for tube 3,5, 8,  change merged_data_T3, to ..._5, ..._8 
al<-merged_data_T3[, c(4, 41, 9)]  #F vs noadults
summary_al <- al %>%
  group_by(day, nototal) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = nototal, y = mean_F, color = sd_F)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)
  
  # Customizing the x-axis scale (setting limits and breaks)
  scale_x_continuous(
    name = "Population size (nototal)",  # Custom x-axis label
    limits = c(0, 600),  # Set limits from 0 to 20
    breaks = seq(0, 600, by = 100)  # Breaks every 5 units
  ) +
  
  # Customizing the y-axis scale (setting limits and breaks)
  scale_y_continuous(
    name = "Mean of F (mean_F)",  # Custom y-axis label
    limits = c(500, 900),  # Set limits from 0 to 100
    breaks = seq(500, 900, by = 100)  # Breaks every 20 units
  ) +
  
  # Customizing the color legend (scale for sd_F)
  scale_color_gradient(
    low = "blue", high = "red",  # Color gradient from blue to red
    limits = c(0, 200),  # Set limits for sd_F from 0 to 30
    breaks = seq(0, 200, by = 20),  # Breaks for the legend at intervals of 10
    name = "Standard Deviation (sd_F)"  # Legend title
  ) +
  
  labs(title = "Scatter Plot of mean_F vs nototal with Custom Scales") +
  theme_minimal()
#-----------------------------------------------------------------------------
##For F vs nototal with equalibrium status color coded
#-----------------------------------------------------------------------------
#for tube 3,5, 8,  change merged_data_T3, to ..._5, ..._8 
al<-merged_data_T3[, c(4, 41,46,9)]  #F vs noadults
summary_al <- al %>%
  group_by(day, nototal, Equ) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = nototal, y = mean_F, color = Equ)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs nototal with Custom Scales") +
  theme_minimal()
##-------------------------------------------------------------------------------

al<-merged_data_T8[, c(4, 42,46,9)]  #F vs noadults
summary_al <- al %>%
  group_by(day, noadults, Equ) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = noadults, y = mean_F, color = Equ)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  
  
  labs(title = "Scatter Plot of mean_F vs noadults with Custom Scales") +
  theme_minimal()

#--------------------------------------------------------------------------

#For FM
#for tube 3,5, 8,  change merged_data_T3, to ..._5, ..._8 
al<-merged_data_T3[, c(4, 42, 10)]  #FM vs noadults
summary_al <- al %>%
  group_by(day, noadults) %>%
  summarize(
    mean_FM = ifelse(all(is.na(FM)), NA, mean(FM, na.rm = TRUE)),  # Mean of FM
    sd_FM = ifelse(all(is.na(FM)), NA, sd(FM, na.rm = TRUE))       # Standard deviation of FM
  )

# Load the ggplot2 package
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_al, aes(x = noadults, y = mean_FM, color = sd_FM)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)
  
  # Customizing the x-axis scale (setting limits and breaks)
  scale_x_continuous(
    name = "Number of Adults (noadults)",  # Custom x-axis label
    limits = c(0, 80),  # Set limits from 0 to 20
    breaks = seq(0, 80, by = 10)  # Breaks every 5 units
  ) +
  
  # Customizing the y-axis scale (setting limits and breaks)
  scale_y_continuous(
    name = "Mean of FM (mean_FM)",  # Custom y-axis label
    limits = c(500, 900),  # Set limits from 0 to 100
    breaks = seq(500, 900, by = 100)  # Breaks every 20 units
  ) +
  
  # Customizing the color legend (scale for sd_FM)
  scale_color_gradient(
    low = "blue", high = "red",  # Color gradient from blue to red
    limits = c(0, 200),  # Set limits for sd_FM from 0 to 30
    breaks = seq(0, 200, by = 20),  # Breaks for the legend at intervals of 10
    name = "Standard Deviation (sd_FM)"  # Legend title
  ) +
  
  labs(title = "Scatter Plot of mean_FM vs noadults with Custom Scales") +
  theme_minimal()


###############################################################################
#statistical test with raw readings
StaTest<-merged_data[,c(4,8,9,42,41,45,46)] 

# Load lme4 package
library(lme4)
#------------------------------------------------------------------------------
# F vs noadults
#------------------------------------------------------------------------------
# Fit a mixed-effects model with tube as random effects, day as fixed
model21 <- lmer(F ~ noadults + day + (1 | tube) , data=StaTest)
model22 <- lmer(F ~ noadults + noadults^2 + day  + (1 | tube) , data=StaTest) 
model23 <- lmer(F ~ noadults + noadults^2 + noadults^3 +day  + (1 | tube) , data=StaTest) 
# Fit a mixed-effects model with tube as random effects, day as fixed, equ as fixed
model31 <- lmer(F ~ noadults + day +Equ + (1 | tube) , data=StaTest) 
model32 <- lmer(F ~ noadults + noadults^2 + day +Equ + (1 | tube) , data=StaTest) 
model33 <- lmer(F ~ noadults + noadults^2 + noadults^3 + day +Equ + (1 | tube) , data=StaTest) 
# Fit a mixed-effects model with tube as random effects, day as random, equ as fixed
model41 <- lmer(F ~ noadults  +Equ + (1 | tube) + (1|day), data=StaTest) 
model42 <- lmer(F ~ noadults + noadults^2  +Equ + (1 | tube) + (1|day), data=StaTest) 
model43 <- lmer(F ~ noadults + noadults^2 + noadults^3 + Equ + (1 | tube)+ (1|day) , data=StaTest) 
AIC(model21, model22,model23, model31, model32, model33, model41, model42, model43)
#------------------------------------------------------------------------------
# F vs nototal
#------------------------------------------------------------------------------
# Fit a mixed-effects model with tube as random effects, day as fixed
model21 <- lmer(F ~ nototal + day + (1 | tube) , data=StaTest)
model22 <- lmer(F ~ nototal + nototal^2 + day  + (1 | tube) , data=StaTest) 
model23 <- lmer(F ~ nototal + nototal^2 + nototal^3 +day  + (1 | tube) , data=StaTest) 
# Fit a mixed-effects model with tube as random effects, day as fixed, equ as fixed
model31 <- lmer(F ~ nototal + day +Equ + (1 | tube) , data=StaTest) 
model32 <- lmer(F ~ nototal + nototal^2 + day +Equ + (1 | tube) , data=StaTest) 
model33 <- lmer(F ~ nototal + nototal^2 + nototal^3 + day +Equ + (1 | tube) , data=StaTest) 
# Fit a mixed-effects model with tube as random effects, day as random, equ as fixed
model41 <- lmer(F ~ nototal  +Equ + (1 | tube) + (1|day), data=StaTest) 
model42 <- lmer(F ~ nototal + nototal^2  +Equ + (1 | tube) + (1|day), data=StaTest) 
model43 <- lmer(F ~ nototal + nototal^2 + nototal^3 + Equ + (1 | tube)+ (1|day) , data=StaTest) 
AIC(model21, model22,model23, model31, model32, model33, model41, model42, model43)

#------------------------------------------------------------------------------
#mean_F vs noadults
#------------------------------------------------------------------------------
#statistical test with mean body size
summary_StaTest <- StaTest %>%
  group_by(day, tube, Equ, noadults) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_StaTest, aes(x = noadults, y = mean_F, color = Equ)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs noadults with Custom Scales") +
  theme_minimal()

ggplot(summary_StaTest, aes(x = noadults, y = mean_F, color = sd_F)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs noadults with Custom Scales") +
  theme_minimal()


# Fit a mixed-effects model with tube as random effects, day as fixed
model21 <- lmer(mean_F ~ noadults + day + (1 | tube) , data=summary_StaTest)
model22 <- lmer(mean_F ~ noadults + noadults^2 + day  + (1 | tube) , data=summary_StaTest) 
model23 <- lmer(mean_F ~ noadults + noadults^2 + noadults^3 +day  + (1 | tube) , data=summary_StaTest) 
# Fit a mixed-effects model with tube as random effects, day as fixed, equ as fixed
model31 <- lmer(mean_F ~ noadults + day +Equ + (1 | tube) , data=summary_StaTest) 
model32 <- lmer(mean_F ~ noadults + noadults^2 + day +Equ + (1 | tube) , data=summary_StaTest) 
model33 <- lmer(mean_F ~ noadults + noadults^2 + noadults^3 + day +Equ + (1 | tube) , data=summary_StaTest) 
# Fit a mixed-effects model with tube as random effects, day as random, equ as fixed
model41 <- lmer(mean_F ~ noadults  +Equ + (1 | tube) + (1|day), data=summary_StaTest) 
model42 <- lmer(mean_F ~ noadults + noadults^2  +Equ + (1 | tube) + (1|day), data=summary_StaTest) 
model43 <- lmer(mean_F ~ noadults + noadults^2 + noadults^3 + Equ + (1 | tube)+ (1|day) , data=summary_StaTest) 
AIC(model21, model22,model23, model31, model32, model33, model41, model42, model43)

#-------------------------------------------------------------------------------
#mean_F vs nototal
#------------------------------------------------------------------------------
#statistical test with mean body size
summary_StaTest <- StaTest %>%
  group_by(day, tube, Equ, nototal) %>%
  summarize(
    mean_F = ifelse(all(is.na(F)), NA, mean(F, na.rm = TRUE)),  # Mean of F
    sd_F = ifelse(all(is.na(F)), NA, sd(F, na.rm = TRUE))       # Standard deviation of F
  )
library(ggplot2)
# Scatter plot with custom axis and color scales
ggplot(summary_StaTest, aes(x = nototal, y = mean_F, color = Equ)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs nototal with Custom Scales") +
  theme_minimal()

ggplot(summary_StaTest, aes(x = nototal, y = mean_F, color = sd_F)) +
  geom_point(size = 3) +  # Scatter plot points
  geom_smooth(method = "loess", se = FALSE, color = "black") +  # Smooth trend line (loess)+
  labs(title = "Scatter Plot of mean_F vs nototal with Custom Scales") +
  theme_minimal()


# Fit a mixed-effects model with tube as random effects, day as fixed
model21 <- lmer(mean_F ~ nototal + day + (1 | tube) , data=summary_StaTest)
model22 <- lmer(mean_F ~ nototal + nototal^2 + day  + (1 | tube) , data=summary_StaTest) 
model23 <- lmer(mean_F ~ nototal + nototal^2 + nototal^3 +day  + (1 | tube) , data=summary_StaTest) 
# Fit a mixed-effects model with tube as random effects, day as fixed, equ as fixed
model31 <- lmer(mean_F ~ nototal + day +Equ + (1 | tube) , data=summary_StaTest) 
model32 <- lmer(mean_F ~ nototal + nototal^2 + day +Equ + (1 | tube) , data=summary_StaTest) 
model33 <- lmer(mean_F ~ nototal + nototal^2 + nototal^3 + day +Equ + (1 | tube) , data=summary_StaTest) 
# Fit a mixed-effects model with tube as random effects, day as random, equ as fixed
model41 <- lmer(mean_F ~ nototal  +Equ + (1 | tube) + (1|day), data=summary_StaTest) 
model42 <- lmer(mean_F ~ nototal + nototal^2  +Equ + (1 | tube) + (1|day), data=summary_StaTest) 
model43 <- lmer(mean_F ~ nototal + nototal^2 + nototal^3 + Equ + (1 | tube)+ (1|day) , data=summary_StaTest) 
AIC(model21, model22,model23, model31, model32, model33, model41, model42, model43)



# Asssumption check
#STEP1-homoscedasticity check  (Residuals vs Fitted values plot)
#you want to see residuals have constant variance across fitted values (homoscedasticity), i.e. no clear pattern or funnel shape.
library(ggplot2)
# Extract fitted values and residuals
fitted_vals <- fitted(model)
residuals_vals <- resid(model)
# Create a residuals vs fitted plot

ggplot(data = NULL, aes(x = fitted_vals, y = residuals_vals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs Fitted") +
  theme_minimal()


#STEP2-Normality of Residuals check (Q-Q Plot)
#  Residuals should ideally fall along the reference line.
ggplot(data = NULL, aes(sample = residuals_vals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Residuals") +
  theme_minimal()

# Another Normality check (Histogram of residuals)
ggplot(data = NULL, aes(x = residuals_vals)) +
  geom_histogram(binwidth = 10, fill = "blue", alpha = 0.7, color = "black") +
  labs(x = "Residuals", y = "Count", title = "Histogram of Residuals") +
  theme_minimal()

#STEP3-Check for Influential Observations (Leverage and Cook’s Distance)
# Calculate Cook's distance for influential points
#This will help you identify any points that have undue influence on your model.
cooksd <- cooks.distance(model)
# Plot Cook's distance
plot(cooksd, type = "h", ylab = "Cook's Distance", main = "Cook's Distance Plot")
abline(h = 4/(nrow(StaTest)-length(fixef(model))), col = "red")  # Threshold line for large Cook's distance

#STEP4- Check for Random Effect Independence
# Plot residuals vs random effects
ranef_tube <- model.frame(model)$tube
ranef_day<-model.frame(model)$day
residuals_vals <- resid(model)
REI<-data.frame("ranef_tube"=ranef_tube,
               "ranef_day"=ranef_day,
               "residual_vals"=residuals_vals)
ggplot(data = REI, aes(x = ranef_tube, y = residuals_vals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Random Effects (tube)", y = "Residuals", title = "Residuals vs Random Effects (tube)") +
  theme_minimal()
ggplot(data = REI, aes(x = ranef_day, y = residuals_vals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Random Effects (tube)", y = "Residuals", title = "Residuals vs Random Effects (tube)") +
  theme_minimal()




model <- lmer(body_size ~ population_density + (1 | tube) + (1 | day/tube), data = StaTest)


# Summary of the model
summary(model)
