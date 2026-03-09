############################################################
# Section 4.4: Local interpretation – Shapley values
# Description:
#   Compute and visualise per-variable Shapley contributions for a single
#   target point, relative to its nearest decision-boundary counterfactual.
############################################################

#############################
# 1) Select target and CF   #
#############################

{
tdp <- 3
biplot_dataT[tdp, ]
#biplot_dataT %>% filter(pred_col == "red")

{
  # Target row (features only)
  pred_data <- biplot_dataT[tdp, 1:num_vars]
  pred_info <- biplot_dataT[tdp, ]

  # Choose counterfactual boundary in X-space
  #  - single-point search:
  # boundary <- biplot_search$B_boundary_x_value[1, 1:num_vars]
  boundary <- biplot_search_T$B_boundary_x_value[1, 1:num_vars]
  #  - direct projection alternative:
  # boundary <- proj_boundary
  #  - multi-point search alternative:
  # boundary <- pc.cf[nearest_pc, ]

  # Difference vector (X): boundary - data
  data_to_boundary <- as.numeric(boundary - pred_data)
}

#############################################
# 2) Run Shapley and collect contributions  #
#############################################

{
  shapley_store <- Run_Shapley(
    model_use, num_vars,
    test_data = pred_data, B_boundary = boundary,   # B_boundary_x also acceptable
    individual = TRUE, tdp = 1,                     # can pass specific SPD index here
    svm_u = svm_u, lda = lda, rpart_used = rpart_used,
    rounding = rounding
  )$shapley_store

  # Single-point (sdp) contribution row
  shapley_cause <- shapley_store[1, ]
}

#############################################
# 3) Predictions and ordering of variables  #
#############################################

{
  # Current point prediction
  pred_value <- pred_function(
    model_use = model_use,
    model_select = model_select,
    rounding = rounding,
    new_data = pred_data
  )

  # Boundary point prediction
  pred_boundary <- pred_function(
    model_use = model_use,
    model_select = model_select,
    rounding = rounding,
    new_data = boundary
  )

  # Class and order: determines decreasing order for plotting
  pred_use <- as.numeric(ifelse(pred_value >= cut_off, 1, 0))     # predicted class (0/1)
  dec_order <- ifelse(pred_use == 1, TRUE, FALSE)                 # TRUE => decreasing
  dist_order <- order(shapley_cause, decreasing = dec_order)      # if class 1, negatives come first
}

#############################################
# 4) Assemble plotting data (no logic change)
#############################################

{
  # Base frame: data value, delta to boundary, and shapley contribution
  df <- as.data.frame(cbind(as.numeric(pred_data), data_to_boundary, shapley_cause))
  df <- df %>% dplyr::rename(pred_data = V1)

  # Standardised distance option (kept as in original; default FALSE)
  df$vec_cause <- df$data_to_boundary
  std_dist <- FALSE
  if (std_dist == TRUE) {
    df$vec_cause <- df$data_to_boundary / X_sd
  }

  # Add names and flags
  df <- cbind(varnames = var_names, dist_order, df)
  df <- df %>% dplyr::mutate(pred_use = pred_use)

  # Label whether each variable supports or contradicts the current class
  df <- df %>% dplyr::mutate(
    Contribute = dplyr::case_when(
      pred_use == 0 & shapley_cause >= 0 ~ "Supports",     # in class 0, increasing causes supports
      pred_use == 0 & shapley_cause <  0 ~ "Contradicts",
      pred_use == 1 & shapley_cause >= 0 ~ "Contradicts",  # in class 1, increasing causes contradicts
      pred_use == 1 & shapley_cause <  0 ~ "Supports",
      TRUE ~ "Green"
    )
  )

  # Apply order
  df <- df[dist_order, , drop = FALSE]
  df <- df %>% dplyr::mutate(varnames = factor(varnames, levels = varnames))  # keep order after sorting

  # Pretty y-axis labels (value -> delta)
  df <- df %>% dplyr::mutate(
    varnames_p = paste0(varnames, ": ", round(pred_data, 2), "-> ", round(data_to_boundary, 2))
  )
  df <- df %>% dplyr::mutate(varnames_p = factor(varnames_p, levels = varnames_p))
}

########################
# 5) Build the graphic #
########################

{
  # Method 1 (length by change) retained but commented in original
  # ggplot1 <- ggplot(df, aes(fill = Contribute, colour = Contribute, x = vec_cause, y = varnames_p)) +
  #   scale_fill_manual(values = c("Supports" = "darkblue", "Contradicts" = "grey")) +
  #   scale_color_manual(values = c("Supports" = "darkblue", "Contradicts" = "grey")) +
  #   geom_col() + theme_bw() +
  #   geom_label(aes(label = round(shapley_cause, rounding), x = 0), colour = "white", show.legend = FALSE) +
  #   labs(
  #     x = paste0("Change to Boundary - ", ifelse(std_dist == TRUE, "Standardised", "Original")),
  #     y = "Variable contribution order",
  #     title = "Distance to Boundary: Contribution & order from Shapley value",
  #     subtitle = paste("ID and class: ", tdp, ",", class_data$class[tdp],
  #                      "; Predicted Class and score: ", class_data$pred_use[tdp], ",",
  #                      round(class_data$pred_value[tdp], 3),
  #                      "; Boundary score: ", round(pred_boundary, 3))
  #   )

  # Method 2 (length by probability change) — used
  shapley_plot <- ggplot(df, aes(fill = Contribute, colour = Contribute,
                                                   x = shapley_cause, y = varnames_p)) +
    scale_fill_manual(name = "Contribution",
                               values = c("Supports" = "darkblue", "Contradicts" = "grey")) +
    scale_color_manual(name = "Contribution",
                                values = c("Supports" = "darkblue", "Contradicts" = "grey")) +
    geom_col() + theme_bw() +
    geom_label(aes(label = round(shapley_cause, rounding),
                                     x = round(shapley_cause, rounding)),
                        colour = "white", show.legend = FALSE) +
    theme(legend.position = "bottom") +
    labs(
      x = "[-> change to reach boundary] : Impact on prediction of data values",
      y = "Variable contribution order",
      title = "Shapley Contribution Plot",
      subtitle = paste(
        "ID and class: ", tdp, ",", pred_info$class,
        "; Predicted Class and score: ", pred_info$pred_use, ",", round(pred_info$pred_value, 3),
        "; Boundary score: ", round(pred_boundary, 3)
      )
    )
}

print(shapley_plot)

}
###################################
# 6) Investigate outcome
## not needed ###
###################################

biplot_dataT[tdp,]
biplot_grid_Test$Z.st[tdp,] %*% tV[1:2,] + X_center
Z_st_subset[tdp,] %*% tV[1:2,] + X_center
Bz_subset[tdp,] %*% tV[1:2,] + X_center
vec_to_boundary_Z[tdp,] %*% tVr

summary(train_data)
vec_to_boundary[tdp,]
X_sd
vec_to_boundary_sd[tdp,]
sort(abs(vec_to_boundary_sd[tdp,]))

p <- 3.4 * -2.77662 + 5.4 * 0.13186 + 6.2 * -0.06556 + 7.5234 
p <- (3.4-0.59*0) * -2.77662 + (5.4+0.16*1) * 0.13186 + (6.2+0.14*1) * -0.06556 + 7.5234 
1/ (1+exp(-p))
summary(model_use)

###################################
# 7) Update a modified data point #

## allows you to see the new data with the select variable updated to the CF value
## Only update the variables to the CF value necessary to change the outcome
## therefore, the variables in blue bars


## plot: 
## plot the updated value onto the rotated biplot. 
## Because updated CF is not on the biplot plane anymore, it is not always in the correct biplot prediction region
##     but should be closer visually
###################################

{
  final_data <- rotate_data[1, 1:num_vars]
  # final_data <- biplot_data[tdp, 1:(num_vars + 2)]   # original alternative

  # Select which factors to change:
  # (Edit indices here as needed)
  var_names
  change_values <- c(5, 2, 4)
  final_data[change_values] <- boundary[change_values]

  pred_value <- pred_function(
    model_use = model_use,
    model_select = model_select,
    rounding = rounding,
    new_data = final_data
  )

  print(final_data)
  print(pred_value)
}

# plot the updated value

point_inter_proj(
  X = final_data,
  V = V,
  plot = TRUE,
  proj = pc.values,
  col_x = "darkgrey",
  cex_x = 1.5,
  standardise = standardise_use,
  CVA_method = CVA_method
)

#########################
####### END #############
#########################