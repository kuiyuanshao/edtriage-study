########################### Calibration Plot ###########################

calibration.plot <- function(data, lab, title, df, xlim, ylim, min, interval){
  lines_x_obs <- rep(data$predRate, 2)
  lines_x_pred <- lines_x_obs + interval
  lines_y_obs <- c(min - data$obsRate * -(abs(min) - (min - ylim[1])), rep(min, 10))
  lines_y_pred <- c(min - (1 - data$obsRate) * -(abs(min) - (min - ylim[1])), rep(min, 10))
  lines_coord <- data.frame(lines_x = c(lines_x_obs, lines_x_pred),
                            lines_y = c(lines_y_obs, lines_y_pred),
                            colour = rep(c("obs", "pred"), each = 20),
                            group = c(rep(LETTERS[1:10], 2), rep(LETTERS[11:20], 2)))
  text_x <- c(data$predRate, data$predRate + interval)
  text_y <- c(min - data$obsRate * -(abs(min) - (min - ylim[1])), 
              min - (1 - data$obsRate) * -(abs(min) - (min - ylim[1])))
  texts <- c(round(data$obsRate * data$obsNo), round((1 - data$obsRate) * data$obsNo))
  texts_coord <- data.frame(text_x = text_x, text_y = text_y,
                            texts = texts, colour = rep(c("obs", "pred"), each = 10))
  
  ggplot_object <- ggplot(data, aes(x = predRate, y = obsRate)) + 
    geom_errorbar(aes(ymin = obsRate_LCL, ymax = obsRate_UCL), width = interval) + 
    xlab(paste("Predicted probability of", lab)) + 
    ylab(paste("Observed probability of", lab)) + 
    geom_line() + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    ggtitle(paste("Calibration Plot for", title)) + 
    theme_bw() + 
    scale_y_continuous(paste("Observed probability of", lab),
                       limits = ylim,
                       expand = c(0, 0),
                       breaks = seq(0, ylim[2], by = ylim[2] / 5),
                       minor_breaks = NULL,
                       sec.axis = sec_axis(~.*(max(data$obsNo) * 5 + max(data$obsNo) * 1.2), 
                                           name = "Number of ED presentations",
                                           breaks = c(-(max(data$obsNo) * 5 + max(data$obsNo) * 1.2)*abs(min), 
                                                      -(max(data$obsNo) * 5 + max(data$obsNo) * 1.2)*(min - ylim[1]), 
                                                      (max(data$obsNo) * 5 + max(data$obsNo) * 1.2) * ylim[2]),
                                           labels = c(0, max(data$obsNo), max(data$obsNo) * 5))) +
    scale_x_continuous(expand = c(0, 0),
                       limits = xlim) + 
    geom_line(data = lines_coord, aes(x = lines_x, y = lines_y, group = group, colour = colour)) + 
    theme(legend.position = "None") + 
    geom_text(data = texts_coord, aes(x = text_x, y = text_y, label = texts, colour = colour), size = 2, hjust = "right")
    
    
  return (ggplot_object)
}


########################### ROC Curve ###########################

roc.curve <- function(response, predictor, varname, multi = F, modelnames = NULL){
  if (multi == F){
    roc_obj <- roc(response, predictor)
    roc_ci <- ci.se(roc_obj, specificities = seq(0, 1, l = 50))
    roc_ci_df <- data.frame(x = as.numeric(rownames(roc_ci)),
                                lower = roc_ci[, 1],
                                upper = roc_ci[, 3])
  
    roc_curve <- ggroc(roc_obj, legacy.axes=TRUE) + 
      theme_minimal() + 
      geom_abline(slope = 1, 
                  linetype = "dashed", 
                  color = "grey") + 
      coord_equal() + 
      geom_ribbon(
        data = roc_ci_df,
        aes(x = 1 - x, ymin = lower, ymax = upper),
        fill = "black",
        alpha = 0.2,
        inherit.aes = F) + 
      ggtitle(paste0("ROC Curve for ", varname))
    return (roc_curve)
  }else{
    formu <- reformulate(sprintf("predictor%s", 
                                 paste0("[[", 1:length(predictor), "]]")), 
                         "response")
    roc_list <- roc(formu)
    roc_list_ci <- lapply(roc_list, ci.se, specificities = seq(0, 1, l = 50))
    
    roc_ci_df <- lapply(roc_list_ci, function(ciobj) 
      data.frame(x = as.numeric(rownames(ciobj)),
                 lower = ciobj[, 1],
                 upper = ciobj[, 3]))
    
    roc_curve <- ggroc(roc_list,legacy.axes=TRUE) + 
      theme_minimal() + 
      geom_abline(slope = 1, 
                  linetype = "dashed", 
                  color = "grey") + 
      coord_equal() 
    for (i in 1:length(predictor)){
      roc_curve <- roc_curve + 
        geom_ribbon(
          data = roc_ci_df[[i]],
          aes(x = 1 - x, ymin = lower, ymax = upper),
          fill = i + 1,
          alpha = 0.2,
          inherit.aes = F)
    }
    roc_curve <- roc_curve + 
      scale_colour_manual(values = (1:length(predictor))[-1], 
                          labels = modelnames, name = "Models")
    return (roc_curve)
  }
}






