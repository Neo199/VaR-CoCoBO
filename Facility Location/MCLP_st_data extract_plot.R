# Extract TRUE BO Coverage (Excluding Initial Samples)
# For all trials and all BO instances

library(dplyr)

# Assuming you have loaded:
# load('results_MCLP/stress_test/stress_test_complete.RData')

# Initialize storage
n_trials <- 5
n_bo_instances <- 10
n_init <- 5  # First 5 are initial samples

# Create dataframe to store results
bo_true_results <- data.frame(
  trial_id = integer(),
  bo_instance = integer(),
  best_overall = numeric(),
  best_from_bo = numeric(),
  best_iteration_overall = integer(),
  best_iteration_bo = integer()
)

# Extract for each trial and each BO instance
for (trial in 1:n_trials) {
  
  cat("\n==============================================\n")
  cat("TRIAL", trial, "\n")
  cat("==============================================\n\n")
  
  for (inst in 1:n_bo_instances) {
    
    # Get coverage history
    coverage_hist <- stress_test_results$all_trial_results[[trial]]$bo_results[[inst]]$coverage_history
    
    # Overall best (including initial samples)
    best_overall <- max(coverage_hist)
    best_iter_overall <- which.max(coverage_hist)
    
    # Best from BO iterations only (positions 6-250)
    coverage_bo_only <- coverage_hist[(n_init + 1):length(coverage_hist)]
    best_from_bo <- max(coverage_bo_only)
    best_iter_bo <- which.max(coverage_bo_only) + n_init  # Adjust for actual position
    
    # Store results
    bo_true_results <- rbind(bo_true_results, data.frame(
      trial_id = trial,
      bo_instance = inst,
      best_overall = best_overall,
      best_from_bo = best_from_bo,
      best_iteration_overall = best_iter_overall,
      best_iteration_bo = best_iter_bo
    ))
    
    cat(sprintf("  Instance %2d: Best Overall = %6.0f (iter %3d) | Best from BO = %6.0f (iter %3d)\n",
                inst, best_overall, best_iter_overall, best_from_bo, best_iter_bo))
  }
  
  # Trial summary
  cat("\n  Trial Summary:\n")
  cat(sprintf("    Best overall (all instances): %.0f\n", 
              max(bo_true_results$best_overall[bo_true_results$trial_id == trial])))
  cat(sprintf("    Best from BO (all instances): %.0f\n", 
              max(bo_true_results$best_from_bo[bo_true_results$trial_id == trial])))
  cat(sprintf("    Mean BO coverage: %.0f ± %.0f\n",
              mean(bo_true_results$best_from_bo[bo_true_results$trial_id == trial]),
              sd(bo_true_results$best_from_bo[bo_true_results$trial_id == trial])))
}

# ---------------------------------------------------------
# AGGREGATE STATISTICS
# ---------------------------------------------------------

cat("\n\n==============================================\n")
cat("AGGREGATE STATISTICS ACROSS ALL TRIALS\n")
cat("==============================================\n\n")

# Summary by trial
trial_summary <- bo_true_results %>%
  group_by(trial_id) %>%
  summarise(
    mean_best_overall = mean(best_overall),
    sd_best_overall = sd(best_overall),
    mean_best_bo = mean(best_from_bo),
    sd_best_bo = sd(best_from_bo),
    max_best_bo = max(best_from_bo),
    min_best_bo = min(best_from_bo)
  )

print(trial_summary)

# ---------------------------------------------------------
# UPDATED COMPARISON WITH GLPK AND GA
# ---------------------------------------------------------

cat("\n\n==============================================\n")
cat("CORRECTED COMPARISON: GLPK vs GA vs BO\n")
cat("==============================================\n\n")

# Get GLPK and GA results from original summary
original_summary <- stress_test_results$summary_df

# Create corrected comparison
corrected_comparison <- data.frame(
  trial_id = 1:n_trials,
  glpk_coverage = original_summary$glpk_coverage,
  ga_coverage = original_summary$ga_coverage,
  bo_mean_coverage = trial_summary$mean_best_bo,
  bo_sd_coverage = trial_summary$sd_best_bo,
  bo_max_coverage = trial_summary$max_best_bo,
  bo_min_coverage = trial_summary$min_best_bo
)

# Calculate gaps
corrected_comparison$ga_gap <- 100 * (corrected_comparison$glpk_coverage - corrected_comparison$ga_coverage) / corrected_comparison$glpk_coverage
corrected_comparison$bo_gap_mean <- 100 * (corrected_comparison$glpk_coverage - corrected_comparison$bo_mean_coverage) / corrected_comparison$glpk_coverage

print(corrected_comparison)

# Overall statistics
cat("\n\nOVERALL STATISTICS:\n")
cat(sprintf("GLPK Mean Coverage: %.0f ± %.0f\n", 
            mean(corrected_comparison$glpk_coverage), 
            sd(corrected_comparison$glpk_coverage)))
cat(sprintf("GA Mean Coverage: %.0f ± %.0f (Gap: %.2f%% ± %.2f%%)\n", 
            mean(corrected_comparison$ga_coverage), 
            sd(corrected_comparison$ga_coverage),
            mean(corrected_comparison$ga_gap),
            sd(corrected_comparison$ga_gap)))
cat(sprintf("BO Mean Coverage: %.0f ± %.0f (Gap: %.2f%% ± %.2f%%)\n", 
            mean(corrected_comparison$bo_mean_coverage), 
            sd(corrected_comparison$bo_mean_coverage),
            mean(corrected_comparison$bo_gap_mean),
            sd(corrected_comparison$bo_gap_mean)))

# ---------------------------------------------------------
# CHECK: Did BO improve beyond initial samples?
# ---------------------------------------------------------

cat("\n\n==============================================\n")
cat("DID BO IMPROVE BEYOND INITIAL SAMPLES?\n")
cat("==============================================\n\n")

improved_count <- sum(bo_true_results$best_from_bo > bo_true_results$best_overall)
total_runs <- nrow(bo_true_results)

if (improved_count == 0) {
  cat("❌ BO NEVER improved beyond initial samples in ANY run!\n")
  cat(sprintf("   0 out of %d runs showed improvement\n", total_runs))
} else {
  cat(sprintf("✓ BO improved in %d out of %d runs (%.1f%%)\n", 
              improved_count, total_runs, 100 * improved_count / total_runs))
}

# Check if best from BO equals best overall (meaning it came from initial samples)
initial_sample_best <- sum(bo_true_results$best_from_bo == bo_true_results$best_overall)
cat(sprintf("\nRuns where best came from initial samples: %d out of %d (%.1f%%)\n",
            initial_sample_best, total_runs, 100 * initial_sample_best / total_runs))

# ---------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------

cat("\n\nSaving corrected results...\n")

write.csv(bo_true_results, "bo_true_results_detailed.csv", row.names = FALSE)
write.csv(corrected_comparison, "corrected_comparison_summary.csv", row.names = FALSE)

cat("\nFiles saved:\n")
cat("  - bo_true_results_detailed.csv\n")
cat("  - corrected_comparison_summary.csv\n")

# ---------------------------------------------------------
# RETURN RESULTS
# ---------------------------------------------------------

# Return as list for further analysis
results <- list(
  detailed = bo_true_results,
  trial_summary = trial_summary,
  corrected_comparison = corrected_comparison
)

cat("\n==============================================\n")
cat("EXTRACTION COMPLETE!\n")
cat("==============================================\n")
