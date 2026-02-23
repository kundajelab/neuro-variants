#############################################################
## 1.  USER-EDITABLE PATHS ##################################
#############################################################

infile <- "/Users/amarderstein/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/chrombpnet_variant_effects/output/pred/results/scores_sfari.txt"

#############################################################
## 2.  LOAD  +  PREP  DATA ##################################
#############################################################
library(data.table)
library(ggplot2)
library(zoo)

dnm <- fread(infile)[
  , .(snp_id,
      status = factor(Pheno, levels = c("control", "case")),
      FLARE  = FLARE_fb,
      # FLARE  = FLARE_heart,
      phylop)
][!is.na(FLARE)]

setorder(dnm, -FLARE, snp_id)  # enforce reproducibility

cat("Loaded", nrow(dnm), "variants (",
    sum(dnm$status == "case"),    "cases, ",
    sum(dnm$status == "control"), "controls)\n")

#############################################################
## 3.  SWEEP OVER FLARE THRESHOLDS ##########################
#############################################################

cut_vec <- seq(0.999, 0.90, by = -0.001)

if ("Proband_ID" %in% names(dnm)) {
  n_cases <- length(unique(dnm[status == "case", Proband_ID]))
} else {
  n_cases <- 1902
}

results <- lapply(cut_vec, function(qc) {
  thr <- quantile(dnm$FLARE, qc, na.rm = TRUE)
  dnm[, high := FLARE >= thr]
  
  tb <- table(dnm$high, dnm$status)
  if (!all(c("TRUE", "FALSE") %in% rownames(tb))) return(NULL)
  
  Ocase <- tb["TRUE", "case"]
  Octrl <- tb["TRUE", "control"]
  Ncase <- sum(tb[, "case"])
  Nctrl <- sum(tb[, "control"])
  
  Ecase <- Ncase * (Octrl / Nctrl)
  excess <- Ocase - Ecase
  
  fish <- fisher.test(tb, alternative = "greater")
  var_excess <- Ocase + Ecase
  ci <- excess + c(-1, 1) * 1.96 * sqrt(var_excess)
  
  data.table(
    cutoff_pct = qc * 100,
    FLARE_threshold = thr,
    Ocase = Ocase,
    Octrl = Octrl,
    expected_case = Ecase,
    excess = excess,
    percent_probands = 100 * excess / n_cases,
    percent_lo       = 100 * ci[1] / n_cases,
    percent_hi       = 100 * ci[2] / n_cases,
    p_fisher         = fish$p.value
  )
})

summary_tab <- rbindlist(results)
print(summary_tab)

# Save enrichment sweep
# fwrite(summary_tab, "flare_quantile_sweep.tsv", sep = "\t")

# ggsave("flare_threshold_curve.pdf", width = 6, height = 4.5)

#############################################################
## 4.  CHOOSE OPTIMAL FLARE THRESHOLD #######################
#############################################################

chosen_row <- summary_tab[
  !is.na(p_fisher) & p_fisher < 0.05,
  .SD[which.max(percent_probands)]
]

if (nrow(chosen_row) == 0) {
  cat("\nNo threshold met significance; falling back to 95th percentile.\n")
  chosen_row <- summary_tab[cutoff_pct == 95]
}

print(chosen_row)

chosen_thresh <- chosen_row$FLARE_threshold
chosen_pct    <- chosen_row$cutoff_pct

cat(sprintf(
  "\n>>> Selected FLARE threshold = %.4f  (top %.2f%% of variants)\n",
  chosen_thresh, chosen_pct
))

ggplot(summary_tab, aes(cutoff_pct, percent_probands)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = percent_lo, ymax = percent_hi), alpha = 0.2) +
  labs(x = "Top X% of FLARE-scored DNMs",
       y = "% of probands with excess high-FLARE DNMs",
       title = "Enrichment of high-FLARE DNMs in cases vs. controls") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_abline(slope=0,intercept=0,lty='dashed',col='red') +
  geom_point(aes(x = chosen_row$cutoff_pct,y=chosen_row$percent_probands),col="steelblue1",lty='dashed',size=rel(5))


dnm[, high_flare := FLARE >= chosen_thresh]

#############################################################
## 5.  PHYLOP-STRATIFIED EXPECTATION & EXCESS ###############
#############################################################

dnm[, phy_bin := cut(
  phylop,
  breaks = quantile(phylop, probs = seq(0, 1, 0.1), na.rm = TRUE),
  include.lowest = TRUE, ordered_result = TRUE)]

ctrl_rate <- dnm[status == "control",
                 .(ctrl_high = sum(high_flare), ctrl_n = .N),
                 by = phy_bin][,
                               bin_rate := ctrl_high / ctrl_n]

case_size <- dnm[status == "case", .N, by = phy_bin]
exp_tbl   <- merge(case_size, ctrl_rate, by = "phy_bin", all.x = TRUE)
exp_tbl[is.na(bin_rate), bin_rate := 0]
exp_tbl[, exp_case_high := N * bin_rate]

expected  <- sum(exp_tbl$exp_case_high)
observed  <- dnm[status == "case" & high_flare == TRUE, .N]
excess    <- observed - expected

var_delta <- observed + expected
ci_delta  <- excess + c(-1, 1) * 1.96 * sqrt(var_delta)

prop_excess <- 100 * excess      / n_cases
prop_lo     <- 100 * ci_delta[1] / n_cases
prop_hi     <- 100 * ci_delta[2] / n_cases

cat("\n===== PHYLOP-STRATIFIED BURDEN =====\n")
cat("Observed high-FLARE DNMs in cases :", observed, "\n")
cat("Expected (phyloP-matched controls) :", formatC(expected, digits = 3), "\n")
cat("Excess                              :", formatC(excess,   digits = 3), "\n")
cat(sprintf("%% of probands with excess DNM    : %.2f %%  (95%% CI %.2f–%.2f %%)\n",
            prop_excess, prop_lo, prop_hi))
cat("====================================\n")
