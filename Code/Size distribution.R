################################################################################
# Density-dependent body-size decomposition using 20 bins along N/K
################################################################################
library(dplyr)
library(tidyr)
library(e1071)
library(ggplot2)
library(moments)

library(readxl)
library(purrr)
library(minpack.lm)
library(e1071)
library(patchwork)
library(here)
library(gratia)



# Define relevant size columns
stage_cols <- c(
  "F", "FM", "SM", "tritos", "protos", "larvae", "eggs",
  "Qtritos", "Qprotos", "Qlarvae", "deutos", "Qdeutos"
)

analyze_body_size_distribution <- function(data, group_by_tube = FALSE) {
  # Define grouping variables
  grouping_vars <- if (group_by_tube) c("day", "tube") else "day"
  
  # Summarize metrics
  summary_df <- data %>%
    select(all_of(c(grouping_vars, "q", stage_cols))) %>%
    pivot_longer(cols = all_of(stage_cols), names_to = "stage", values_to = "size") %>%
    filter(!is.na(size)) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarize(
      mean_size = mean(size, na.rm = TRUE),
      cv_size   = sd(size, na.rm = TRUE) / mean(size, na.rm = TRUE),
      skew_size = skewness(size, na.rm = TRUE),
      q   = mean(q, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_df$tube<-as.factor(summary_df$tube)
  
  library(mgcv)
  
  #analyse relationship between mean body size and relative density
  gam_mean <- gam(
    mean_size ~ s(q) + tube,
    data = summary_df,
    method = "REML"
  )
  
  #analyse relationship between cv of body size and relative density
  gam_cv <- gam(
    cv_size ~ s(q) + tube,
    data = summary_df,
    method = "REML"
  )
  #analyse relationship between skewness of body size and relative density
  gam_skew <- gam(
    skew_size ~ s(q) + tube,
    data = summary_df,
    method = "REML"
  )
  
 #predict the fitted mean body size 
  newdat_size <- data.frame(
    q = seq(min(summary_df$q, na.rm = TRUE),
            max(summary_df$q, na.rm = TRUE),
            length.out = 200),
    tube = factor(levels(summary_df$tube)[1], levels = levels(summary_df$tube))
  )
  
  pred_size <- predict(
    gam_mean,
    newdata = newdat_size,
    se.fit = TRUE
  )
  
  newdat_size$fit <- pred_size$fit
  newdat_size$lwr <- pred_size$fit - 1.96 * pred_size$se.fit
  newdat_size$upr <- pred_size$fit + 1.96 * pred_size$se.fit
  
  #predict the fitted CV of body size 
  newdat_cv <- data.frame(
    q = seq(min(summary_df$q, na.rm = TRUE),
            max(summary_df$q, na.rm = TRUE),
            length.out = 200),
    tube = factor(levels(summary_df$tube)[1], levels = levels(summary_df$tube))
  )
  
  pred_cv <- predict(
    gam_cv,
    newdata = newdat_cv,
    se.fit = TRUE
  )
  
  newdat_cv$fit <- pred_cv$fit
  newdat_cv$lwr <- pred_cv$fit - 1.96 * pred_cv$se.fit
  newdat_cv$upr <- pred_cv$fit + 1.96 * pred_cv$se.fit
  
  
  #predict the fitted CV of body size 
  newdat_skew <- data.frame(
    q = seq(min(summary_df$q, na.rm = TRUE),
            max(summary_df$q, na.rm = TRUE),
            length.out = 200),
    tube = factor(levels(summary_df$tube)[1], levels = levels(summary_df$tube))
  )
  
  pred_skew <- predict(
    gam_skew,
    newdata = newdat_skew,
    se.fit = TRUE
  )
  
  newdat_skew$fit <- pred_skew$fit
  newdat_skew$lwr <- pred_skew$fit - 1.96 * pred_skew$se.fit
  newdat_skew$upr <- pred_skew$fit + 1.96 * pred_skew$se.fit
  
  
  # Set color aesthetics if tube included
  aes_color <- if (group_by_tube) aes(color = as.factor(tube)) else NULL
  
  # Plot 1: Mean size
  p1 <- ggplot(summary_df, aes(x = q, y = mean_size)) +
    geom_point(aes_color, alpha = 0.6) +
    geom_ribbon(
      data = newdat_size,
      aes(
        x = q,
        ymin = lwr,
        ymax = upr
      ),
      inherit.aes = FALSE,
      alpha = 0.2,
      fill = "grey70"
    )+
    geom_line(
      data = newdat_size,
      aes(
        x = q,
        y = fit
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1
    )+
    theme_minimal() +
    theme(
      legend.position = c(0.4, 1.2),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = unit(2, "mm"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
      ) + 
    labs(
      title = "",
      y =expression("Mean ("~mu~"m)"),
      color = if (group_by_tube) "Tube" else NULL
    )+
    scale_x_continuous(
      breaks = seq(0, 1.6, by = 0.2))
  
  # Plot 2: CV
  p2 <- ggplot(summary_df, aes(x = q, y = cv_size)) +
    geom_point(aes_color, alpha = 0.6) +
    geom_ribbon(
      data = newdat_cv,
      aes(
        x = q,
        ymin = lwr,
        ymax = upr
      ),
      inherit.aes = FALSE,
      alpha = 0.2,
      fill = "grey70"
    )+
    geom_line(
      data = newdat_cv,
      aes(
        x = q,
        y = fit
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1
    )+
    theme_minimal() +
    theme(
    legend.position = if (group_by_tube) "none" else NULL,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.4),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(2, "mm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()) + 
    labs(
      title = "",
      y = "CV",
      color = if (group_by_tube) "Tube" else NULL
    )+  
    scale_x_continuous(
      breaks = seq(0, 1.6, by = 0.2))
  
  # Plot 3: Skewness
  p3 <- ggplot(summary_df, aes(x = q, y = skew_size)) +
    geom_point(aes_color, alpha = 0.6) +
    geom_ribbon(
      data = newdat_skew,
      aes(
        x = q,
        ymin = lwr,
        ymax = upr
      ),
      inherit.aes = FALSE,
      alpha = 0.2,
      fill = "grey70"
    )+
    geom_line(
      data = newdat_skew,
      aes(
        x = q,
        y = fit
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1
    )+
    theme_minimal() +
    theme (
      legend.position = if (group_by_tube) "none" else NULL,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black"),
      axis.ticks.length = unit(2, "mm")
    ) +
    labs(
      title = "",
      x = "Relative density (N/K)", y = "Skewness",
      color = if (group_by_tube) "Tube" else NULL)+
    scale_x_continuous(
      breaks = seq(0, 1.6, by = 0.2))
  
  return(list(
    summary = summary_df,
    plot_mean = p1,
    plot_cv = p2,
    plot_skew = p3
  ))
  
}

plot_distr_metric<-analyze_body_size_distribution(merged_data,TRUE) 
plot_distr1 <- (
  plot_distr_metric[["plot_mean"]] /
    plot_distr_metric[["plot_cv"]] /
    plot_distr_metric[["plot_skew"]]
) +
  plot_annotation(tag_levels = 'a')

ggsave ("Results/body size distribution metrics vs relative population size.tiff", 
        plot_distr1, width = 6, height = 4, dpi = 600, compression = "lzw")



#################################EXTRA CODE####################################
plot_body_size_distribution_by_day <- function(data, 
                                               stage_cols, 
                                               group_by_tube = FALSE,
                                               bins = 30,
                                               day_range = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(moments)  # for skewness
  
  # Step 1: Reshape to long format
  long_data <- data %>%
    select(day, tube, nototal, all_of(stage_cols)) %>%
    pivot_longer(cols = all_of(stage_cols), names_to = "stage", values_to = "size") %>%
    filter(!is.na(size))
  
  # Step 2: Filter days if specified
  if (!is.null(day_range)) {
    long_data <- long_data %>%
      filter(day >= day_range[1], day <= day_range[2])
  }
  
  # Step 3: Summarize for annotation
  summary_stats <- long_data %>%
    group_by(day) %>%
    summarize(
      mean_nototal = mean(nototal, na.rm = TRUE),
      se_nototal   = sd(nototal, na.rm = TRUE) / sqrt(n()),
      mean_size    = mean(size, na.rm = TRUE),
      cv_size      = sd(size, na.rm = TRUE) / mean(size, na.rm = TRUE),
      skew_size    = skewness(size, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("Day = ", day,
                          "\nPop size = ", round(mean_nototal, 1), " ± ", round(se_nototal, 1),
                          "\nMean = ", round(mean_size, 2),
                          "\nCV = ", round(cv_size, 2),
                          "\nSkew = ", round(skew_size, 2)))
  
  # Step 4: Build plot
  p <- ggplot(long_data, aes(x = size)) +
    facet_wrap(~ day, scales = "fixed", ncol = 6) +
    theme_minimal() +
    theme(strip.text = element_blank()) +
    labs(x = "Body size", y = "Count", title = "") +
    geom_text(data = summary_stats, 
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.2,
              size = 3.1)
  
  xrange <- c(0,950)
  yrange <- c(0, 30)
  
  # Step 5: Histogram logic
  if (group_by_tube) {
    p <- p +
      geom_histogram(aes(fill = as.factor(tube)), binwidth = 0.5, alpha = 0.6, position = "identity") +
      scale_fill_discrete(name = "Tube")
    
  } else {
    p <- p + 
      scale_x_continuous(limits = xrange) +
      scale_y_continuous(limits = yrange)+
      geom_histogram(breaks = seq(xrange[1], xrange[2], length.out = bins+1),
                     fill = "grey60", color = "white")
  }
  
  return(p)
}

long_data <- data %>%
  select(day, tube, nototal, all_of(stage_cols)) %>%
  pivot_longer(cols = all_of(stage_cols), names_to = "stage", values_to = "size") %>%
  filter(!is.na(size))

long_data$stage<-factor(long_data$stage, levels = c("F", "FM", "SM","Qtritos", "Qprotos", "Qdeutos", "Qlarvae",
                                                    "tritos", "deutos", "protos","larvae","eggs"))
ggplot(long_data,aes(x = size, fill = stage)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.4, color = NA) +
  labs(x = "Body size", y = "Count", fill = "Stage")


  plot_distr_201 <- plot_body_size_distribution_by_day(
    data = merged_data,
    stage_cols = stage_cols,
    bins = 30,
    group_by_tube = FALSE,
    day_range = c(200, 204)
  )
  
  plot_distr_26 <- plot_body_size_distribution_by_day(
    data = merged_data,
    stage_cols = stage_cols,
    bins = 30,
    group_by_tube = FALSE,
    day_range = c(25, 26)
  )
  
  plot_distr_6 <- plot_body_size_distribution_by_day(
    data = merged_data,
    stage_cols = stage_cols,
    bins = 30,
    group_by_tube = FALSE,
    day_range = c(5, 6)
  )
  
  
  plot_distr_merge <- (
    plot_distr_6|
      plot_distr_26|
      plot_distr_201
  ) +
    plot_annotation(tag_levels = 'a')
  
  ggsave ("Results/body size distribution change.tiff", 
          plot_distr_merge, width = 6, height = 4, dpi = 600, compression = "lzw")
  
  
  