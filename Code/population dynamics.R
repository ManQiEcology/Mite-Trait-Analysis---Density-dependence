rm(list=ls())

library(here)
library(readxl)
library(nls2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Metrics)     # for rmse()
library(yardstick)
library(tidyverse)
library(patchwork)
library(forcats)

setwd(here())

################################################################################
#prepare the dataset for control tubes
################################################################################
#load the abundance data
CountData=read_xlsx("Data/MQ test_ Count and Size data.xlsx", sheet = "CountData")
colnames(CountData)<-c(colnames(CountData)[1:23],"noeggs")

# get the data for control
CountData_C<-subset(CountData,treatment=="C")
CountData_C$tube<-as.factor(CountData_C$tube)

###############################################################################
#mark the data before equlibrium and after equlibrium and trend of population dynamic
###############################################################################

#logistic regression of population dynamic
#-------------------------------------------------------------------------------
# Define the logistic model function
logistic_model <- function(day, K, N0, r) {
  K / (1 + ((K - N0) / N0) * exp(-r * day))
}


#fit the logistic model and plot the population dynamic
# Create an empty list to store plots
day_eq<-NULL
plot_list <- list()
fit_model<-list()
CountData_C_extd <- NULL
   for (tubeID in c(3,5,8)) {
     CountData_CTtube<-subset(CountData,tube==tubeID)
     
     # Initial parameter guesses
     start_params <- list(K = max(CountData_CTtube$nototal), N0 = min(CountData_CTtube$nototal), r = 0.1)
     
     # Fit the logistic model
     N_fit <- nls2(nototal ~ logistic_model(day, K, N0, r),
                 data = CountData_CTtube,
                 start = start_params,
                 algorithm = "port")
     fit_model[[tubeID]]<-N_fit
    
      # View the fit
      summary(N_fit)
      # Add fitted values to the data frame
      CountData_CTtube$N_fit <- predict(N_fit, newdata = CountData_CTtube$nototal) 
      CountData_C_extd<-rbind(CountData_C_extd,CountData_CTtube)
      #get the day that population achieve equalibrium (95% of carrying capacity)
      params <- coef(N_fit)
      K <- params["K"]
      target_N <- 0.99 * K
      day_eq_tubeID <- CountData_CTtube$day[min(which(CountData_CTtube$N_fit>=target_N))]
      day_eq <- c(day_eq,day_eq_tubeID)
      
      # Plot the time series with the fitted logistic model
      p<-ggplot(CountData_CTtube, aes(x = day)) +
         geom_line(aes(y = nototal), color = "blue", size = 1) +
         geom_line(aes(y = N_fit), color = "red", linetype = "dashed") +
        labs(title = NULL, 
             x = "Day", 
             y = ifelse(tubeID==3, "Population abundance", "")) +
        ylim(0,750)+
        theme(
          panel.background = element_blank(),      # Remove panel background
          panel.grid.minor = element_blank(),      # Remove minor grid lines
          panel.grid.major = element_line(color = "grey80"))  # Optional: keep major grid lines
      
      plot_list[[tubeID]] <- p
   }

names(day_eq)<-c("tube_3","tube_5","tube_8")

# Combine into a 1x3 layout
combined_plot_for_population_dynamics <- plot_list[[3]] + plot_list[[5]] + plot_list[[8]]+
plot_annotation(tag_levels = "a")  # Automatically labels "a", "b", "c"
# Print the layout

ggsave("Results/Population_dynamics_fitted_with_logistic_model.tiff", 
       combined_plot_for_population_dynamics, width = 6, height = 3, dpi = 600, compression = "lzw")
  
#-------------------------------------------------------------------------------
   #mark the equilibrium status on CountData_C_extd
   for (tubeID in c(3,5,8)) {
     rowID<-which(CountData_C_extd$tube==tubeID)
     
     #get the threshold of abundance beyond which equilibrium was achieved
     day_eq_tubeID<-day_eq[paste("tube_",tubeID,sep = "")]
     
     for (i in 1:(length(rowID))) {
       ifelse(CountData_C_extd$day[rowID[i]]<=day_eq_tubeID, 
              CountData_C_extd$Equ[rowID[i]]<-"Pre-equ", 
              CountData_C_extd$Equ[rowID[i]]<-"Post-equ")  
     }
   }

#----------------------------------------------------------------------------   
#calculate lambda of the regressed logistic growth over time

for (tubeID in c(3,5,8)) {
  
  # Extract parameters
  N_fit<-fit_model[[tubeID]]
  params <- coef(N_fit)
  K <- params["K"]
  N0 <- params["N0"]
  r <- params["r"]
  
  # Time points for which you want to calculate lambda
  time_points <- 1:365  # Example time points
  
  # Calculate population sizes at different times
  Abundance_regrs <- logistic_model(time_points, K, N0, r)
  
  # Calculate finite growth rate (lambda) for each time point
  lambda_regrs<- tail(Abundance_regrs, -1) / head(Abundance_regrs, -1)
  
  # Plot the relationship between Lambda and Population Size
  if (tubeID==3) {
  plot(head(Abundance_regrs,-1), 
       lambda_regrs, 
       type="l",
       xlab = expression("Population size at time " * t),
       ylab = expression("Realized " * lambda * " at time " * t),
       col="grey90")
    text(300, 1.08, "tube 3", col="grey90")
  }
  else {
  lines(head(Abundance_regrs,-1), 
        lambda_regrs,
        col=ifelse(tubeID==5,"grey60","grey30")) 
    text(x=317,y=ifelse(tubeID==5,1.075,1.070), labels=tubeID,
         col=ifelse(tubeID==5,"grey60","grey30"))
  }
}
  
#----------------------------------------------------------------------------
#Calculate the lambda of observed data in the experiment and assign it to CountDate_C
interpolated_data<-NULL
       for (tubeID in c(3,5,8)) { #t= 5， 8 for tube 5 and 8
       rowID<-which(CountData_C_extd$tube==tubeID) 
       
       sub_pop <- data.frame(
       Day = CountData_C_extd$day[rowID],  
       Population_Size = CountData_C_extd$nototal[rowID] 
       )
     all_days <- seq(min(sub_pop$Day), max(sub_pop$Day), by = 1)  # Creates a sequence of continuous days
     # Perform linear interpolation using approx()
     interpolated_population_size <- approx(x = sub_pop$Day, y = sub_pop$Population_Size, xout = all_days)
     
     # Create a new data frame with interpolated values and realized lambda
     interpolated_data_tubeID <- data.frame(
       day = interpolated_population_size$x,
       population_Size = interpolated_population_size$y,
       lambda = c(tail(interpolated_population_size$y,-1)/
                    head(interpolated_population_size$y,-1),
                  NA),
       tube = tubeID
     )
     interpolated_data<-rbind(interpolated_data,interpolated_data_tubeID)
       }
#add realized lambda to original dataset: CountData_C_extd
CountData_C_extd$lambda<-NULL
  for (i in 1:dim(interpolated_data[1])) {
    day <- interpolated_data$day[i]
    tubeID <- interpolated_data$tube[i]
    row_IDD <- which(CountData_C_extd$day==day & CountData_C_extd$tube==tubeID)
    CountData_C_extd$lambda[row_IDD]<-interpolated_data$lambda[i]
  }

#--------------------------------------------------------------------------------  
#analyze relationship between population growth rate and population size 

#calculate r with lambda
CountData_C_extd$r<-log(CountData_C_extd$lambda)

plot(data_clean$nototal, data_clean$r)

#Remove rows with missing Lambda values
data_clean <- CountData_C_extd %>%
       filter(!is.na(lambda))

# Define outlier condition — adjust based on what you found
outlier_condition <- data_clean$r > 0.4

#remove the outlier for regression model fitting
data_model <- data_clean[!outlier_condition, ]

# Bin into deciles (you can increase this to 20 for more granularity)
data_model <- data_model %>%
  mutate(bin = cut(nototal, breaks = 20)) %>%
  group_by(bin) %>%
  mutate(bin_count = n()) %>%
  ungroup()

# Compute meaningful weights: more weight for sparse bins
max_count <- max(data_model$bin_count)
data_model <- data_model %>%
  mutate(weight = max_count / bin_count)

# Replicate rows according to weight
data_rep <- data_model[rep(1:nrow(data_model), round(data_model$weight)),]

  # Model 1: Ricker Model (with only negative density dependent effect)
  model_Ricker<- nlsLM(r ~ d + a * (nototal+c)*exp(-b*(nototal+c)),
                        data = data_rep,
                        start = list(a = 0.005, b=0.01, c = -100, d = -0.05),
                       #lower = c(a = 0.003, b = 0.008, c = -150, d = -0.03),
                       #upper = c(a = 0.007, b = 0.012, c = -50, d = -0.07),
                       control = nls.lm.control(maxiter = 500))
  
  # Model 2: Modified Beverton-Halt Model(with Allee effect and negative density dependent effect)
  model_Mberverton_halt <- nlsLM(
    r ~ e + (a * (nototal)^(d - 1)) / (1 + b * (nototal)^d),
    data = data_rep_sub,
    start = list(a = 0.0006, b = 0.00001, d = 3, e = -0.15),
    #lower = c(a = 0.0004, b = 0.000008, d = 2.4, e = -0.19),
    #upper = c(a = 0.0008, b = 0.000012, d = 3.6, e = -0.11),
    control = nls.lm.control(maxiter = 500)
  )
  
  
  # Predict r from all models
  data_clean <- data_clean %>%
    mutate(pred_Ricker = predict(model_Ricker, newdata=data_clean),
           pred_M_BH   = predict(model_Mberverton_halt, newdata=data_clean))
  
  # Define RMSE manually
  rmse <- function(actual, predicted) {
    sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  }
  
  # Apply it to your models
  rmse_Ricker <- rmse(log(data_clean$lambda), data_clean$pred_Ricker)
  rmse_M_BH    <- rmse(log(data_clean$lambda), data_clean$pred_M_BH)
  
 
  
  print(paste("RMSE - Ricker: ", round(rmse_Ricker, 4)))
  print(paste("RMSE - Modified Beverton Halt:", round(rmse_M_BH, 4)))
  
  data_clean$tube<-as.factor(data_clean$tube)
  # Plot results
  r_vs_N <- ggplot(data_clean, aes(x = nototal, color = tube)) +
    geom_point(aes(y = r), alpha = 0.4) +
    geom_line(aes(y = pred_Ricker), color = "blue", linetype = "dashed") +
    geom_line(aes(y = pred_M_BH), color = "green") +
    geom_point(data = subset(data_clean, r > 0.4 & nototal < 200), 
               aes(x = nototal, y = r), 
               color = "black", size = 3, shape = 1)+
    labs(title =NULL,
         y = "Population growth rate", 
         x = "Population size",
         fill= "Population ID") +
    theme_minimal()+   # Use a minimal base theme
    theme(
      legend.position = c(0.9,0.7),
      panel.background = element_blank(),         # Remove panel background
      panel.grid.minor = element_blank(),         # Remove minor grid lines
      panel.grid.major = element_blank(),         # Optionally remove major grid
      axis.line = element_line(color = "black"),  # Draw x and y axes
      axis.ticks = element_line(color = "black")  # Show axis ticks
    )
  
  ggsave("Results/r_vs_N.tiff", 
         r_vs_N, width = 3, height = 2, dpi = 600, compression = "lzw")
  
################################################################################
#plot population structure dynamic across tubes
##################################################################

# Reshape the data to long format
long_data <- CountData_C_extd %>%
  pivot_longer(cols = c( noeggs, nojuveniles, noadults),  # Specify all stage group columns here
               names_to = "stage_group", 
               values_to = "count")
  
# Calculate total count per tube and day
long_data <- long_data %>%
  group_by(tube, day) %>%
  mutate(total_count = sum(count, na.rm = TRUE),
         proportion = count / total_count) %>%
  ungroup()

#replace stage_group names with the ones more reader friendly 
long_data$stage_group <- fct_recode(long_data$stage_group,
                                    "Eggs" = "noeggs",
                                    "Juveniles" = "nojuveniles",
                                    "Adults" = "noadults")

# Plot totalno vs day, with totalno composed of individuals from each stage group, and tube as a factor
# Prepare a unique dataset for N_fit line (one per day & tube)
line_data <- long_data %>%
  select(day, tube, N_fit) %>%
  distinct()
# Plot
count_stage_group <- ggplot() +
  # Stacked bar by stage group
  geom_bar(data = long_data,
           aes(x = day, y = count, fill = stage_group),
           stat = "identity", position = "stack", alpha = 1) +
  
  # Smooth dashed line for total fit
  geom_line(data = line_data,
            aes(x = day, y = N_fit),
            color = "black", linetype = "dashed", size = 0.9) +
  
  # Facet by tube
  facet_wrap(~ tube, scales = "free_x") +
  
  # Labels and styling
  labs(
    x = NULL,
    y = "Number of individuals \n from each stage group",
    fill = "Stage group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),         # Remove panel background
    panel.grid.minor = element_blank(),         # Remove minor grid lines
    panel.grid.major = element_blank(),         # Optionally remove major grid
    axis.line = element_line(color = "black"),  # Draw x and y axes
    axis.ticks = element_line(color = "black")  # Show axis ticks
  )



proportion_stage_group1<- 
  ggplot()+
  geom_bar(data = long_data,
           aes(x = day, y = proportion, fill = stage_group),
           stat = "identity", position = "stack", alpha = 1) +
  facet_wrap(~ tube,scales = "free_x") + 
  labs(title = NULL,
       x = "Day",
       y = "Proportion \n of each stage group",
       fill = "Stage group") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(),
    legend.position = "top",
    panel.background = element_blank(),         # Remove panel background
    panel.grid.minor = element_blank(),         # Remove minor grid lines
    panel.grid.major = element_blank(),         # Optionally remove major grid
    axis.line = element_line(color = "black"),  # Draw x and y axes
    axis.ticks = element_line(color = "black")  # Show axis ticks
    )
proportion_stage_group2<- 
  ggplot()+
  geom_line(data = long_data,
           aes(x = day, y = proportion, colour  = stage_group)) +
  facet_wrap(~ tube,scales = "free_x") + 
  labs(title = NULL,
       x = "Day",
       y = "Proportion \n of each stage group",
       color = "Stage group") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(),
    legend.position = c(0.95,0.91),
    panel.background = element_blank(),         # Remove panel background
    panel.grid.minor = element_blank(),         # Remove minor grid lines
    panel.grid.major = element_blank(),         # Optionally remove major grid
    axis.line = element_line(color = "black"),  # Draw x and y axes
    axis.ticks = element_line(color = "black")  # Show axis ticks
  )

# combine the count and proportion plot
population_structure1 <- count_stage_group / proportion_stage_group1+  # Stack vertically
plot_annotation(tag_levels = "a")  # Automatically labels "a", "b", "c"
# Save to file
population_structure2 <- count_stage_group / proportion_stage_group2+  # Stack vertically
  plot_annotation(tag_levels = "a")  # Automatically labels "a", "b", "c"
# Save to file
ggsave("Results/population_structure1.tiff", 
       population_structure1, width = 7, height = 5, dpi = 600, compression = "lzw")
ggsave("Results/population_structure2.tiff", 
       population_structure2, width = 7, height = 5, dpi = 600, compression = "lzw")

###########################################################################

