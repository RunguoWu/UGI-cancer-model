# Old validation functions
# For archive only

pred_d4s <- NULL

for (st in sites) {

  optimized_params <- as.numeric(params[params$site==st, grepl("tp", colnames(params))])

  sub_d4s <- subset(d4s, site==st)

  sub_summary <- sub_d4s %>%
    group_by(month, diagnosed_stage) %>%
    summarise(.groups = 'drop')

  sub_summary$correct_prob_diag <- 0
  sub_summary$correct_prob_b4diag <- 0
  sub_summary$most_likely_state <- 0
  sub_summary$most_likely_state_1monthb4 <- 0


  for (i in 1:nrow(sub_summary)) {

    month_diagnosis <- sub_summary$month[i]
    diagnosed_stage <- sub_summary$diagnosed_stage[i]

    # in the month of observed diagnosis
    correct_prob_diag <- predict_stage_at_diagnosis(optimized_params,
                                     month_diagnosis = month_diagnosis,
                                     diagnosed_stage = diagnosed_stage,
                                     before_diag = FALSE
    )

    most_likely_state <- predict_most_likely_state_index(optimized_params,
                                                         month_diagnosis = month_diagnosis,
                                                         diagnosed_stage = diagnosed_stage,
                                                         before_diag = FALSE
    )

    # 1 month before observed diagnosis month
    correct_prob_b4diag <- predict_stage_at_diagnosis(optimized_params,
                                                    month_diagnosis = month_diagnosis,
                                                    diagnosed_stage = diagnosed_stage,
                                                    before_diag = TRUE
    )

    most_likely_state_1monthb4 <- predict_most_likely_state_index(optimized_params,
                                                         month_diagnosis = month_diagnosis,
                                                         diagnosed_stage = diagnosed_stage,
                                                         before_diag = TRUE
    )

    sub_summary[i, "correct_prob_diag"] <- correct_prob_diag
    sub_summary[i, "correct_prob_b4diag"] <- correct_prob_b4diag
    sub_summary[i, "most_likely_state"] <- most_likely_state
    sub_summary[i, "most_likely_state_1monthb4"] <- most_likely_state_1monthb4

  }

  sub_d4s <- sub_d4s %>% left_join(sub_summary, by = c("month", "diagnosed_stage"))

  pred_d4s <- rbind(pred_d4s, sub_d4s)

}

d4s2 <- d4s %>% left_join(pred_d4s[, c("e_patid", "correct_prob_diag", "correct_prob_b4diag",
                                      "most_likely_state", "most_likely_state_1monthb4"
                                      )])


create_validation_dotplot(d4s2)
create_accuracy_summary(d4s2)
create_stage_comparison(d4s2)

table(d4s2$diagnosed_stage+4==d4s2$most_likely_state)

table(d4s2$diagnosed_stage==d4s2$most_likely_state_1monthb4)

# the most_likely_state result shows that
# misprediction only happen between diagnosed and undiagnosed at the same state
#

rt <- estimate_starting_distribution(d4s, params)
starting_stage_probs <- rt$stom

P <- create_transition_matrix(
  optimized_params[1], optimized_params[2], optimized_params[3],
  optimized_params[4], optimized_params[5], optimized_params[6], optimized_params[7]
)

simulate_single_patient(P, starting_stage_probs)