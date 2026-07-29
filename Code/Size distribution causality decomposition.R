################################################################################
# Tube-specific moving-window decomposition of body-size distributions
# ----------------------------------------------------------------------
# For each tube:
#   1) order censuses by day
#   2) build overlapping windows of consecutive censuses
#   3) compare window i with window i+1
#   4) decompose the change in mean / CV / skewness into:
#        - developmental-stage composition
#        - within-stage body-size plasticity
#   5) smooth contribution curves with GAMs for plotting
################################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(minpack.lm)
library(e1071)
library(mgcv)
library(patchwork)
library(here)

setwd(here())

################################################################################
# 1) Read data
################################################################################

wb <- "Data/MQ test_ Count and Size data.xlsx"

count_raw <- read_xlsx(wb, sheet = "CountData") %>%
  filter(treatment == "C")

size_raw <- read_xlsx(wb, sheet = "SizeData") %>%
  filter(treatment == "C")

################################################################################
# 2) Standardize stage columns
################################################################################

stage_cols <- c(
  "eggs",
  "larvae", "Qlarvae",
  "protos", "Qprotos",
  "tritos", "Qtritos",
  "deutos", "Qdeutos",
  "F", "FM", "SM"
)

count_df <- count_raw %>%
  transmute(
    tube    = factor(tube),
    day     = as.numeric(day),
    nototal = as.numeric(nototal),
    
    eggs    = `noeggs...15`,
    larvae  = nolarvae,
    Qlarvae = noQL,
    protos  = noprotos,
    Qprotos = noQP,
    tritos  = notritos,
    Qtritos = noQT,
    deutos  = nodeutos,
    Qdeutos = noQD,
    F       = noF,
    FM      = noFM,
    SM      = noSM
  )

stopifnot(all(stage_cols %in% names(count_df)))

size_long <- size_raw %>%
  select(tube, day, all_of(stage_cols)) %>%
  mutate(
    tube = factor(tube),
    day  = as.numeric(day)
  ) %>%
  pivot_longer(
    cols = all_of(stage_cols),
    names_to = "stage",
    values_to = "size"
  ) %>%
  filter(!is.na(size)) %>%
  mutate(stage = factor(stage, levels = stage_cols))

################################################################################
# 3) Fit tube-specific logistic trajectories to obtain relative density q = N/K
################################################################################

logistic_model <- function(day, K, N0, r) {
  K / (1 + ((K - N0) / N0) * exp(-r * day))
}

fit_logistic_by_tube <- function(dat) {
  dat <- dat %>% arrange(day)
  
  start_params <- list(
    K  = max(dat$nototal, na.rm = TRUE),
    N0 = max(min(dat$nototal, na.rm = TRUE), 1),
    r  = 0.1
  )
  
  fit <- nlsLM(
    nototal ~ logistic_model(day, K, N0, r),
    data = dat,
    start = start_params,
    control = nls.lm.control(maxiter = 500)
  )
  
  pars <- coef(fit)
  K_hat <- unname(pars["K"])
  
  dat <- dat %>%
    mutate(
      N_fit = predict(fit, newdata = dat),
      K_hat = K_hat,
      q     = nototal / K_hat
    )
  
  list(fit = fit, data = dat, coef = pars)
}

fit_list <- split(count_df, count_df$tube) %>%
  lapply(fit_logistic_by_tube)

count_fit <- bind_rows(lapply(fit_list, `[[`, "data")) %>%
  mutate(tube = factor(tube))

tube_params <- tibble(
  tube = names(fit_list),
  K    = map_dbl(fit_list, ~ unname(coef(.x$fit)["K"])),
  N0   = map_dbl(fit_list, ~ unname(coef(.x$fit)["N0"])),
  r    = map_dbl(fit_list, ~ unname(coef(.x$fit)["r"]))
)

################################################################################
# 4) Helper functions for counts, pools, and synthetic populations
################################################################################

calc_metrics <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(c(mean = NA_real_, cv = NA_real_, skew = NA_real_))
  }
  
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  c(
    mean = m,
    cv   = ifelse(is.na(m) || m == 0, NA_real_, s / m),
    skew = e1071::skewness(x, type = 2, na.rm = TRUE)
  )
}

get_window_counts <- function(df_tube, idx, stage_cols) {
  counts <- colSums(df_tube[idx, stage_cols, drop = FALSE], na.rm = TRUE)
  counts <- as.numeric(counts)
  names(counts) <- stage_cols
  counts
}

get_window_pools <- function(size_long, tube_i, days_vec, stage_cols, fallback_pools) {
  tmp <- size_long %>%
    filter(tube == tube_i, day %in% days_vec) %>%
    group_by(stage) %>%
    summarise(pool = list(size), .groups = "drop")
  
  out <- setNames(vector("list", length(stage_cols)), stage_cols)
  
  for (st in stage_cols) {
    hit <- tmp$pool[tmp$stage == st]
    if (length(hit) > 0 && length(hit[[1]]) > 0) {
      out[[st]] <- hit[[1]]
    } else {
      out[[st]] <- fallback_pools[[st]]
      if (is.null(out[[st]])) out[[st]] <- numeric(0)
    }
  }
  
  out
}

allocate_counts <- function(counts, Nsim = 2000) {
  counts <- as.numeric(counts)
  names(counts) <- stage_cols
  counts <- pmax(counts, 0)
  
  if (sum(counts, na.rm = TRUE) == 0) {
    out <- rep(0, length(counts))
    names(out) <- stage_cols
    return(out)
  }
  
  props <- counts / sum(counts, na.rm = TRUE)
  expected <- props * Nsim
  
  alloc <- floor(expected)
  remainder <- Nsim - sum(alloc)
  
  if (remainder > 0) {
    frac <- expected - alloc
    idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
    alloc[idx] <- alloc[idx] + 1
  }
  
  names(alloc) <- stage_cols
  alloc
}

sample_population_metrics <- function(counts, pools, Nsim = 2000) {
  alloc <- allocate_counts(counts, Nsim = Nsim)
  x <- unlist(lapply(stage_cols, function(st) {
    n <- alloc[[st]]
    pool <- pools[[st]]
    if (is.na(n) || n <= 0 || length(pool) == 0) return(numeric(0))
    #randomly sampling individuals from size distribution of a given developmental stage
    sample(pool, size = n, replace = TRUE) 
  }), use.names = FALSE)
  calc_metrics(x)
}

simulate_scenario <- function(counts, pools, B = 1000, Nsim = 2000) {
  sims <- replicate(B, sample_population_metrics(counts, pools, Nsim = Nsim))
  rowMeans(sims, na.rm = TRUE)
}

################################################################################
# 5) Overall fallback pools across all censuses
################################################################################

overall_pools <- size_long %>%
  group_by(stage) %>%
  summarise(pool = list(size), .groups = "drop") %>%
  { setNames(.$pool, .$stage) }

################################################################################
# 6) Build overlapping windows within each tube
################################################################################

window_n <- 5   # value of window_n was determined by the code "Detecting suitable combination of window length and step size.R"
step_n <- 2     # value of step_n was determined by the code "Detecting suitable combination of window length and step size.R"

make_tube_windows <- function(df, window_n = 5, step_n = 5) {
  df <- df %>% arrange(day)
  n <- nrow(df)
  
  if (n < window_n + step_n) return(NULL)
  
  map_dfr(seq(1, n - window_n - step_n + 1, by = step_n), function(i) {
    from_idx <- i:(i + window_n - 1)
    to_idx   <- (i + step_n):(i + step_n + window_n - 1)
    
    tibble(
      tube = df$tube[1],
      from_start = i,
      from_end   = i + window_n - 1,
      to_start    = i + step_n,
      to_end      = i + step_n + window_n - 1,
      day_from_mid = mean(df$day[from_idx], na.rm = TRUE),
      day_to_mid   = mean(df$day[to_idx], na.rm = TRUE),
      q_from       = mean(df$q[from_idx], na.rm = TRUE),
      q_to         = mean(df$q[to_idx], na.rm = TRUE)
    )
  })
}

window_pairs <- count_fit %>%
  arrange(tube, day) %>%
  group_split(tube) %>%
  map_dfr(~ make_tube_windows(.x, window_n = window_n, step_n = step_n))

################################################################################
# 7) Counterfactual decomposition for one window transition
################################################################################

decompose_window_transition <- function(tube_i,
                                        from_start, from_end,
                                        to_start, to_end,
                                        count_fit,
                                        size_long,
                                        stage_cols,
                                        overall_pools,
                                        B = 1000,
                                        Nsim = 2000) {
  
  df_tube <- count_fit %>%
    filter(tube == tube_i) %>%
    arrange(day)
  
  if (nrow(df_tube) < to_end) return(NULL)
  
  from_idx <- from_start:from_end
  to_idx   <- to_start:to_end
  
  from_counts <- get_window_counts(df_tube, from_idx, stage_cols)
  to_counts   <- get_window_counts(df_tube, to_idx, stage_cols)
  
  from_days <- df_tube$day[from_idx]
  to_days   <- df_tube$day[to_idx]
  
  from_pools <- get_window_pools(size_long, tube_i, from_days, stage_cols, overall_pools)
  to_pools   <- get_window_pools(size_long, tube_i, to_days, stage_cols, overall_pools)
  
  # baseline = window i
  m0  <- simulate_scenario(from_counts, from_pools, B = B, Nsim = Nsim)
  
  # structure-only = stage composition from window i+1, body sizes from window i
  ms  <- simulate_scenario(to_counts, from_pools, B = B, Nsim = Nsim)
  
  # plasticity-only = stage composition from window i, body sizes from window i+1
  mp  <- simulate_scenario(from_counts, to_pools, B = B, Nsim = Nsim)
  
  # full = window i+1
  msp <- simulate_scenario(to_counts, to_pools, B = B, Nsim = Nsim)
  
  # symmetric decomposition
  S <- 0.5 * ((ms - m0) + (msp - mp))
  P <- 0.5 * ((mp - m0) + (msp - ms))
  I <- msp - m0 - S - P
  
  total <- msp - m0
  abs_sum <- abs(S) + abs(P) + abs(I)
  
  q_mid <- mean(c(
    mean(df_tube$q[from_idx], na.rm = TRUE),
    mean(df_tube$q[to_idx],   na.rm = TRUE)
  ), na.rm = TRUE)
  
  tibble(
    tube = tube_i,
    day_from_mid = mean(from_days, na.rm = TRUE),
    day_to_mid   = mean(to_days, na.rm = TRUE),
    q_mid = q_mid,
    metric = c("Mean", "CV", "Skewness"),
    baseline = as.numeric(m0),
    structure_only = as.numeric(ms),
    plasticity_only = as.numeric(mp),
    full = as.numeric(msp),
    total_change = as.numeric(total),
    structure_effect = as.numeric(S),
    plasticity_effect = as.numeric(P),
    interaction_effect = as.numeric(I),
    structure_share_abs = ifelse(abs_sum == 0, NA_real_, 100 * abs(S) / abs_sum),
    plasticity_share_abs = ifelse(abs_sum == 0, NA_real_, 100 * abs(P) / abs_sum),
    interaction_share_abs = ifelse(abs_sum == 0, NA_real_, 100 * abs(I) / abs_sum)
  )
}

set.seed(123)

decomp_results <- pmap_dfr(
  window_pairs,
  ~ decompose_window_transition(
    tube_i = ..1,
    from_start = ..2,
    from_end   = ..3,
    to_start   = ..4,
    to_end     = ..5,
    count_fit = count_fit,
    size_long = size_long,
    stage_cols = stage_cols,
    overall_pools = overall_pools,
    B = 1000,
    Nsim = 2000
  )
)

decomp_results$metric <- factor(
  decomp_results$metric,
  levels = c("Mean", "CV", "Skewness")
)

################################################################################
# Panel d-f: decomposition of structure vs plasticity by metric
################################################################################
library(dplyr)
library(mgcv)
library(ggplot2)
library(patchwork)

# long format
decomp_long <- decomp_results %>%
  pivot_longer(
    cols = c(structure_share_abs, plasticity_share_abs),
    names_to = "component",
    values_to = "share"
  ) %>%
  mutate(
    component = recode(
      component,
      structure_share_abs  = "Structure",
      plasticity_share_abs = "Plasticity"
    ),
    component = factor(component, levels = c("Structure", "Plasticity")),
    metric = factor(metric, levels = c("Mean", "CV", "Skewness")),
    tube = factor(tube)
  )

# one GAM per metric x component
fit_one <- function(df) {
  gam(
    share ~ s(q_mid, k = 5) + tube,
    data = df,
    method = "REML"
  )
}

gam_list <- list(
  Mean = list(
    Structure  = fit_one(filter(decomp_long, metric == "Mean", component == "Structure")),
    Plasticity = fit_one(filter(decomp_long, metric == "Mean", component == "Plasticity"))
  ),
  CV = list(
    Structure  = fit_one(filter(decomp_long, metric == "CV", component == "Structure")),
    Plasticity = fit_one(filter(decomp_long, metric == "CV", component == "Plasticity"))
  ),
  Skewness = list(
    Structure  = fit_one(filter(decomp_long, metric == "Skewness", component == "Structure")),
    Plasticity = fit_one(filter(decomp_long, metric == "Skewness", component == "Plasticity"))
  )
)

# prediction helper
make_pred_df <- function(model, df) {
  q_seq <- seq(
    min(df$q_mid, na.rm = TRUE),
    max(df$q_mid, na.rm = TRUE),
    length.out = 200
  )
  
  ref_tube <- levels(df$tube)[1]
  
  newdat <- data.frame(
    q_mid = q_seq,
    tube = factor(ref_tube, levels = levels(df$tube))
  )
  
  pred <- predict(model, newdata = newdat, se.fit = TRUE)
  
  newdat$fit <- pred$fit
  newdat$lwr <- pred$fit - 1.96 * pred$se.fit
  newdat$upr <- pred$fit + 1.96 * pred$se.fit
  
  newdat
}

# build prediction data
pred_df <- bind_rows(
  lapply(names(gam_list), function(m) {
    bind_rows(
      lapply(names(gam_list[[m]]), function(comp) {
        df_m <- filter(decomp_long, metric == m, component == comp)
        make_pred_df(gam_list[[m]][[comp]], df_m) %>%
          mutate(metric = m, component = comp)
      })
    )
  })
) %>%
  mutate(
    metric = factor(metric, levels = c("Mean", "CV", "Skewness")),
    component = factor(component, levels = c("Structure", "Plasticity"))
  )

make_panel <- function(metric_name, ylab, show_x = FALSE, show_legend = FALSE) {
  p <- ggplot(
    decomp_long %>% filter(metric == metric_name),
    aes(x = q_mid, y = share, color = component)
  ) +
    geom_point(alpha = 0.35, size = 1.4) +
    geom_ribbon(
      data = pred_df %>% filter(metric == metric_name),
      aes(x = q_mid, ymin = lwr, ymax = upr, fill = component),
      inherit.aes = FALSE,
      alpha = 0.15,
      color = NA
    ) +
    geom_line(
      data = pred_df %>% filter(metric == metric_name),
      aes(x = q_mid, y = fit, color = component),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    scale_x_continuous(breaks = seq(0, 1.6, by = 0.2)) +
    labs(
      x = if (show_x) "Relative density (N/K)" else NULL,
      y = ylab,
      color = NULL,
      fill = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.title.x = if (show_x) element_text() else element_blank(),
      legend.position = if (show_legend) c(0.4, 1.08) else "none",
      legend.justification = c(0.5, 0),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
  
  p
}


p_d <- make_panel("Mean", "Contribution (%)\n(Mean)", show_x = FALSE, show_legend = TRUE)
p_e <- make_panel("CV", "Contribution (%)\n(CV)", show_x = FALSE, show_legend = FALSE)
p_f <- make_panel("Skewness", "Contribution (%)\n(Skewness)", show_x = TRUE, show_legend = FALSE)


Fig_2<-(plot_distr1|(p_d / p_e / p_f))+
  plot_annotation(tag_levels = "a")  # Automatically labels "a", "b", "c"
ggsave("Results/population level body size distribution and contributors.tiff", 
       Fig_2, width = 7, height = 5, dpi = 600, compression = "lzw")


