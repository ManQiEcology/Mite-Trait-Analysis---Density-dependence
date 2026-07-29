##########################################3
#find the best combination of window length and step size
##########################################

library(dplyr)
library(tidyr)
library(purrr)
library(mgcv)
library(ggplot2)

################################################################################
# Candidate grid
################################################################################

candidate_grid <- tidyr::expand_grid(
  window_n = 3:8,
  step_n   = 1:7
) %>%
  filter(step_n < window_n)   # step must be smaller than window

################################################################################
# Helper: score one (window_n, step_n) combination
################################################################################
# Assumes you already have:
#   - count_fit
#   - size_long
#   - stage_cols
#   - overall_pools
#   - make_tube_windows()
#   - decompose_window_transition()
#
# If your helper function names differ, just replace them below.

score_window_combo <- function(window_n,
                               step_n,
                               count_fit,
                               size_long,
                               stage_cols,
                               overall_pools,
                               B = 300,
                               Nsim = 2000) {
  
  # Build moving windows for each tube
  window_pairs <- count_fit %>%
    arrange(tube, day) %>%
    group_split(tube) %>%
    map_dfr(~ make_tube_windows(.x, window_n = window_n, step_n = step_n))
  
  if (is.null(window_pairs) || nrow(window_pairs) < 5) {
    return(tibble(
      window_n = window_n,
      step_n = step_n,
      n_transitions = nrow(window_pairs),
      mean_ribbon_width = NA_real_,
      mean_dev_expl = NA_real_,
      mean_rmse = NA_real_,
      mean_roughness = NA_real_
    ))
  }
  
  # Run decomposition for each transition
  decomp_results <- pmap_dfr(
    window_pairs,
    ~ decompose_window_transition(
      tube_i      = ..1,
      from_start  = ..2,
      from_end    = ..3,
      to_start    = ..4,
      to_end      = ..5,
      count_fit   = count_fit,
      size_long   = size_long,
      stage_cols  = stage_cols,
      overall_pools = overall_pools,
      B           = B,
      Nsim        = Nsim
    )
  )
  
  if (is.null(decomp_results) || nrow(decomp_results) == 0) {
    return(tibble(
      window_n = window_n,
      step_n = step_n,
      n_transitions = nrow(window_pairs),
      mean_ribbon_width = NA_real_,
      mean_dev_expl = NA_real_,
      mean_rmse = NA_real_,
      mean_roughness = NA_real_
    ))
  }
  
  # Long format for panel d-f
  decomp_long <- decomp_results %>%
    pivot_longer(
      cols = c(structure_share_abs, plasticity_share_abs),
      names_to = "component",
      values_to = "share"
    ) %>%
    mutate(
      component = recode(
        component,
        structure_share_abs = "structure",
        plasticity_share_abs = "plasticity"
      ),
      component = factor(component, levels = c("structure", "plasticity")),
      metric = factor(metric, levels = c("Mean", "CV", "Skewness"))
    ) %>%
    filter(is.finite(q_mid), is.finite(share))
  
  # Fit GAMs and calculate smoothness metrics
  fit_one <- function(dat) {
    if (nrow(dat) < 10 || n_distinct(dat$q_mid) < 5) {
      return(tibble(
        dev_expl = NA_real_,
        rmse = NA_real_,
        mean_ribbon_width = NA_real_,
        roughness = NA_real_
      ))
    }
    
    # small basis dimension to avoid overfitting
    k_use <- min(5, max(3, floor(n_distinct(dat$q_mid) / 2)))
    
    m <- gam(
      share ~ s(q_mid, k = k_use) + tube,
      data = dat,
      method = "REML"
    )
    
    pred_obs <- predict(m, newdata = dat)
    rmse <- sqrt(mean((dat$share - pred_obs)^2, na.rm = TRUE))
    
    grid <- data.frame(
      q_mid = seq(min(dat$q_mid, na.rm = TRUE),
                  max(dat$q_mid, na.rm = TRUE),
                  length.out = 200),
      tube = factor(levels(dat$tube)[1], levels = levels(dat$tube))
    )
    
    pr <- predict(m, newdata = grid, se.fit = TRUE)
    
    tibble(
      dev_expl = summary(m)$dev.expl,
      rmse = rmse,
      mean_ribbon_width = mean(2 * 1.96 * pr$se.fit, na.rm = TRUE),
      roughness = mean(abs(diff(dat$share[order(dat$q_mid)])), na.rm = TRUE)
    )
  }
  
  gam_stats <- decomp_long %>%
    group_by(metric, component) %>%
    group_modify(~ fit_one(.x)) %>%
    ungroup()
  
  tibble(
    window_n = window_n,
    step_n = step_n,
    n_transitions = nrow(window_pairs),
    n_points = nrow(decomp_results),
    mean_ribbon_width = mean(gam_stats$mean_ribbon_width, na.rm = TRUE),
    mean_dev_expl = mean(gam_stats$dev_expl, na.rm = TRUE),
    mean_rmse = mean(gam_stats$rmse, na.rm = TRUE),
    mean_roughness = mean(gam_stats$roughness, na.rm = TRUE)
  )
}

################################################################################
# Run the grid search
################################################################################

set.seed(123)

search_results <- pmap_dfr(
  list(candidate_grid$window_n, candidate_grid$step_n),
  ~ score_window_combo(
    window_n = ..1,
    step_n = ..2,
    count_fit = count_fit,
    size_long = size_long,
    stage_cols = stage_cols,
    overall_pools = overall_pools,
    B = 300,      # smaller B for search speed
    Nsim = 2000
  )
)

################################################################################
# Rank candidate combinations
################################################################################

search_ranked <- search_results %>%
  arrange(
    mean_ribbon_width,          # smallest ribbon width first
    desc(mean_dev_expl),        # higher explained deviance is better
    mean_rmse,                  # lower RMSE is better
    mean_roughness              # lower raw roughness is better
  )

print(search_ranked)

best_combo <- search_ranked %>% slice(1)
print(best_combo)

################################################################################
# Optional: visualize the search surface
################################################################################

ggplot(search_results, aes(x = window_n, y = step_n, fill = mean_ribbon_width)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", na.value = "grey90") +
  labs(
    x = "Window length (censuses)",
    y = "Step size (censuses)",
    fill = "Mean ribbon\nwidth"
  ) +
  theme_classic()

ggsave("Results/optimal combination of window length and step size.tiff", 
       Fig_2, width = 7, height = 5, dpi = 600, compression = "lzw")

################################################################################
# After choosing the best combo, rerun the decomposition
# with a larger B (for final figure)
################################################################################

# Example:
# final_window_n <- best_combo$window_n
# final_step_n   <- best_combo$step_n
#
# Then rerun the decomposition with:
# B = 1000 or more