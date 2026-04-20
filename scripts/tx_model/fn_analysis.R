# Plot flexible parametric survival model
plot_stpm2 <- function(fit, by_stage=0, max_time=8, dt, add_legend=TRUE){
  
  time_seq <- seq(0, max_time, length.out=100)  # More points for smoother curves
  
  if (by_stage==0){
    
    pred_mat <- matrix(NA, nrow = length(time_seq), ncol = nrow(dt))
    
    for(i in 1:length(time_seq)){
      # Create newdata with the time variable set
      newdata <- dt
      newdata$fu_diag2cens <- time_seq[i]
      pred_mat[i, ] <- predict(fit, newdata=newdata, type="surv")
    }
    avg_pred <- rowMeans(pred_mat, na.rm=TRUE)
    
    km <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, dt)

    # Create data frames for ggplot
    km_df <- data.frame(time=km$time, surv=km$surv, type="Observed")
    pred_df <- data.frame(time=time_seq, surv=avg_pred, type="Predicted")
    
    p <- ggplot() +
      geom_step(data=km_df, aes(x=time, y=surv, color="Overall", linetype="Observed"), size=1) +
      geom_line(data=pred_df, aes(x=time, y=surv, color="Overall", linetype="Predicted"), size=1) +
      labs(x="Time", y="Survival") +
      theme_bw() +
      ylim(0, 1)
    
  } else {
    dt1 <- subset(dt, stage_imp=="1")
    dt2 <- subset(dt, stage_imp=="2")
    dt3 <- subset(dt, stage_imp=="3")
    dt4 <- subset(dt, stage_imp=="4")
    
    # Predictions for each stage - matrices for all individuals
    pred_mat1 <- matrix(NA, nrow = length(time_seq), ncol = nrow(dt1))
    pred_mat2 <- matrix(NA, nrow = length(time_seq), ncol = nrow(dt2))
    pred_mat3 <- matrix(NA, nrow = length(time_seq), ncol = nrow(dt3))
    pred_mat4 <- matrix(NA, nrow = length(time_seq), ncol = nrow(dt4))
    
    for(i in 1:length(time_seq)){
      newdata1 <- dt1
      newdata1$fu_diag2cens <- time_seq[i]
      pred_mat1[i, ] <- predict(fit, newdata=newdata1, type="surv")
      
      newdata2 <- dt2
      newdata2$fu_diag2cens <- time_seq[i]
      pred_mat2[i, ] <- predict(fit, newdata=newdata2, type="surv")
      
      newdata3 <- dt3
      newdata3$fu_diag2cens <- time_seq[i]
      pred_mat3[i, ] <- predict(fit, newdata=newdata3, type="surv")
      
      newdata4 <- dt4
      newdata4$fu_diag2cens <- time_seq[i]
      pred_mat4[i, ] <- predict(fit, newdata=newdata4, type="surv")
    }
    
    avg_pred1 <- rowMeans(pred_mat1, na.rm = T)
    avg_pred2 <- rowMeans(pred_mat2, na.rm = T)
    avg_pred3 <- rowMeans(pred_mat3, na.rm = T)
    avg_pred4 <- rowMeans(pred_mat4, na.rm = T)
    
    # KM curves
    km1 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, dt1)
    km2 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, dt2)
    km3 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, dt3)
    km4 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, dt4)
    
    # Create data frames for ggplot
    plot_data <- rbind(
      data.frame(time=km1$time, surv=km1$surv, stage="Stage 1", type="Observed"),
      data.frame(time=km2$time, surv=km2$surv, stage="Stage 2", type="Observed"),
      data.frame(time=km3$time, surv=km3$surv, stage="Stage 3", type="Observed"),
      data.frame(time=km4$time, surv=km4$surv, stage="Stage 4", type="Observed"),
      data.frame(time=time_seq, surv=avg_pred1, stage="Stage 1", type="Predicted"),
      data.frame(time=time_seq, surv=avg_pred2, stage="Stage 2", type="Predicted"),
      data.frame(time=time_seq, surv=avg_pred3, stage="Stage 3", type="Predicted"),
      data.frame(time=time_seq, surv=avg_pred4, stage="Stage 4", type="Predicted")
    )
    
    # Color scheme
    stage_colors <- c("Stage 1"="forestgreen", "Stage 2"="blue", 
                      "Stage 3"="orange", "Stage 4"="red")
    
    p <- ggplot(plot_data, aes(x=time, y=surv, color=stage, linetype=type)) +
      geom_line(size=1) +
      scale_color_manual(values=stage_colors) +
      scale_linetype_manual(values=c("Observed"=3, "Predicted"=1)) +
      labs(x="Time", y="Survival", color="Stage", linetype="Type") +
      theme_bw() +
      ylim(0, 1) +
      xlim(0, max_time)
    
    if (!add_legend) {
      p <- p + theme(legend.position="none")
    }
  }
  
  return(p)
}


# Function to calculate stage distribution under partial interventions
partial_intervention_summary <- function(data, 
                                         pct_img = 0, 
                                         pct_2ww = 0, 
                                         seed = NULL) {
  
  # Check that percentages sum to <= 1
  if (pct_img + pct_2ww > 1) {
    stop("Sum of pct_img and pct_2ww cannot exceed 1")
  }
  
  pct_observed = 1 - pct_img - pct_2ww
  
  if (!is.null(seed)) set.seed(seed)
  
  data %>%
    group_by(site) %>%
    mutate(
      # Randomly assign intervention status (img, 2ww, or observed)
      intervention_group = sample(c("img", "2ww", "observed"), 
                                  n(), 
                                  replace = TRUE, 
                                  prob = c(pct_img, pct_2ww, pct_observed))
    ) %>%
    mutate(
      # Assign probabilities based on intervention group
      pred_stage1 = case_when(
        intervention_group == "img" ~ X1_img,
        intervention_group == "2ww" ~ X1_2ww,
        intervention_group == "observed" ~ as.numeric(stage_imp == 1)
      ),
      pred_stage2 = case_when(
        intervention_group == "img" ~ X2_img,
        intervention_group == "2ww" ~ X2_2ww,
        intervention_group == "observed" ~ as.numeric(stage_imp == 2)
      ),
      pred_stage3 = case_when(
        intervention_group == "img" ~ X3_img,
        intervention_group == "2ww" ~ X3_2ww,
        intervention_group == "observed" ~ as.numeric(stage_imp == 3)
      ),
      pred_stage4 = case_when(
        intervention_group == "img" ~ X4_img,
        intervention_group == "2ww" ~ X4_2ww,
        intervention_group == "observed" ~ as.numeric(stage_imp == 4)
      ),
      # Calculate predicted 5-year survival for each patient
      pred_surv = case_when(
        intervention_group == "img" ~ X1_img * diag_s1 + X2_img * diag_s2 + 
          X3_img * diag_s3 + X4_img * diag_s4,
        intervention_group == "2ww" ~ X1_2ww * diag_s1 + X2_2ww * diag_s2 + 
          X3_2ww * diag_s3 + X4_2ww * diag_s4,
        intervention_group == "observed" ~ surv_cur
      )
    ) %>%
    summarise(
      pct_img = pct_img,
      pct_2ww = pct_2ww,
      pct_observed = pct_observed,
      # Observed proportions at each stage
      obs_stage1 = mean(stage_imp == 1),
      obs_stage2 = mean(stage_imp == 2),
      obs_stage3 = mean(stage_imp == 3),
      obs_stage4 = mean(stage_imp == 4),
      # Mean stage probabilities under partial interventions
      pred_stage1 = mean(pred_stage1),
      pred_stage2 = mean(pred_stage2),
      pred_stage3 = mean(pred_stage3),
      pred_stage4 = mean(pred_stage4),
      # 5-year survival
      obs_surv = mean(surv_cur),
      pred_surv = mean(pred_surv),
      # Number in each group
      n_img = sum(intervention_group == "img"),
      n_2ww = sum(intervention_group == "2ww"),
      n_observed = sum(intervention_group == "observed"),
      n_total = n(),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = matches("(obs|pred)_stage\\d"),
      names_to = c("type", "stage"),
      names_pattern = "(.*)_stage(.*)",
      values_to = "proportion"
    ) %>%
    pivot_wider(
      names_from = type,
      values_from = proportion
    ) %>%
    mutate(
      dif_stage = pred - obs,
      dif_surv = pred_surv -obs_surv 
    ) %>% 
    rename(obs_stage = obs,
           pred_stage = pred) %>% 
    select(site, stage, n_total, pct_img, pct_2ww, pct_observed, 
           obs_stage, pred_stage, dif_stage, obs_surv, pred_surv, dif_surv)
}


