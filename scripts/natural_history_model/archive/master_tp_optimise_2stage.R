# Master file to call fn_markov_short.R
rm(list = ls())
library(tidyverse)

source("/data/WIPH-CanDetect/HealthEco/route.R")
source("/data/WIPH-CanDetect/HealthEco/scripts/fn_tp_optimise_2stage.R")

dt <- readRDS(file.path(wd, "pop_ugi_symp_action_20250708.rds"))

d4s <- dt %>% filter(stage != "") %>% mutate(
  age_group3 = case_when(age_index < 60 ~ "<60",
                         age_index >= 60 & age_index < 80 ~ "60-79",
                         age_index >=80 ~ ">=80"
                         ),
  age_group3 = factor(age_group3, levels = c("<60", "60-79", ">=80")),
  month = floor(time2diag/30.5),
  diagnosed_stage = if_else(stage %in% c("1", "2"), 3, 4)
) %>% select(e_patid, month, diagnosed_stage, site, age_group3)


rt <- estimate_parameters(d4s,
                    initial_params = c(0.2, 0.1, 0.15),
                    lower_params = c(0.001, 0.001, 0.001),
                    upper_params = c(0.999, 0.999, 0.999)
)


# By age and site ---------------------------------------------------------

d4s <- d4s

# Set optimization bounds
lower_params <- c(0.001, 0.001, 0.001)
upper_params <- c(0.999, 0.999, 0.999)

# Initialize results dataframe
results <- data.frame(
  site = character(),
  age_group3 = character(),
  n_patients = integer(),
  tp1 = numeric(),
  tp2 = numeric(),
  tp3 = numeric(),
  likelihood = numeric(),
  convergence = integer(),
  stringsAsFactors = FALSE
)

# Get unique combinations
combinations <- d4s %>%
  distinct(site, age_group3) %>%
  arrange(site, age_group3)

print(paste("Processing", nrow(combinations), "combinations"))

# Process each combination
for (i in 1:nrow(combinations)) {
  site_val <- combinations$site[i]
  age_val <- combinations$age_group3[i]
  
  print(paste("Processing:", site_val, "-", age_val))
  
  # Filter data for this combination
  subset_data <- d4s %>%
    filter(site == site_val, age_group3 == age_val)
  
  # Skip if too few patients
  if (nrow(subset_data) < 10) {
    print(paste("Skipping - only", nrow(subset_data), "patients"))
    next
  }
  
  # Get suggested starting parameters
  initial_params <- suggest_starting_parameters(subset_data)
  
  # Run optimization
  result <- estimate_parameters(
    patient_data = subset_data,
    initial_params = initial_params,
    lower_params = lower_params,
    upper_params = upper_params,
    save_interval = NA,
    filename_base = paste0("opt_", site_val, "_", age_val),
    save_final = FALSE
  )
  
  # Store results
  results <- rbind(results, data.frame(
    site = site_val,
    age_group3 = age_val,
    n_patients = nrow(subset_data),
    tp1 = result$final_params[1],
    tp2 = result$final_params[2],
    tp3 = result$final_params[3],
    likelihood = result$final_likelihood,
    convergence = result$convergence_code,
    stringsAsFactors = FALSE
  ))
}

# Display results
print("Optimization Results:")
print(results)

# Save results
write.csv(results, file.path(output, "optimization_results.csv"), row.names = FALSE)





