library("readxl")
library("dplyr")
library("lubridate")
library("BMR")
library("ggplot2")


set.seed(42)

# Load data
file_path <- "C:\\Users\\Halil Tunçer\\Desktop\\data_combined.xlsx"
usa <- read_excel(file_path)

# Preprocess data
usa <- na.omit(usa)
usa$Date <- as.Date(usa$Date, format = "%Y.%m.%d")
usa$CPI <- log(usa$CPI)
usa$Industrial_Production <- log(usa$Industrial_Production)
usa$Interest_Rate <- lag(usa$Interest_Rate, 2)
usa$CPI <- lag(usa$CPI, 2)
usa$Industrial_Production <- lag(usa$Industrial_Production, 2)
colnames(usa) <- c("Date", "return", "dvd_yield", "pb_ratio", "pe_ratio", "ipi", "cpi", "rate")
usa <- na.omit(usa)

# Define subsets
subset1 <- usa[usa$Date <= as.Date("2010-12-31"), ]
subset2 <- usa[usa$Date >= as.Date("2011-01-01") & usa$Date <= as.Date("2023-12-31"), ]
subset3 <- usa

# Function to compute DA
compute_DA <- function(actual_returns, predicted_returns) {
  # Compute P hat (actual hit ratio)
  P_hat <- ifelse(actual_returns * predicted_returns > 0, 1, 0)
  
  # Compute proportions of positive actual and predicted returns
  P_r <- mean(ifelse(actual_returns > 0, 1, 0))
  P_r_hat <- mean(ifelse(predicted_returns > 0, 1, 0))
  
  # Compute P hat star (expected hit ratio under independence)
  P_hat_star <- P_r * P_r_hat + (1 - P_r) * (1 - P_r_hat)
  
  # Compute the variances
  var_P_hat <- var(P_hat)
  var_P_hat_star <- P_hat_star * (1 - P_hat_star)
  
  # Compute DA
  DA <- (mean(P_hat) - P_hat_star) / sqrt(var_P_hat - var_P_hat_star)
  
  return(DA)
}

# Function to compute EP statistics
compute_EP_statistics <- function(y_hat, y) {
  T <- length(y) # Number of observations
  p_hat_y <- 0.5 * (1 + sum(sign(y_hat)) / T) # Estimate of p_hat_y
  var_y <- var(y) # Variance of y
  
  # Compute EP variance estimate
  V_EP <- 4 / T^2 * p_hat_y * (1 - p_hat_y) * sum((y - mean(y))^2)
  
  # Compute EP statistic
  A_T <- mean(sign(y_hat) * y)
  B_T <- (mean(sign(y_hat)) * mean(y))
  EP_statistic <- (A_T - B_T) / sqrt(V_EP)
  
  return(EP_statistic)
}

perform_analysis1 <- function(data, window_size) {
  forecast_results <- data.frame()
  actual_values <- data.frame()
  dates <- data.frame()
  
  # Iterate over the data with a moving window
  for (i in 1:(nrow(data) - window_size)) {
    # Extract the current window
    
    current_window <- data[i:(i + window_size - 1), ]
    
    # Build BVAR model for the current window
    bvar_obj <- new(bvarm)
    bvar_obj$build(data.matrix(current_window[,-1]),
                   TRUE, # constant
                   1)    # lags
    prior <- c(1,1,1,1,1,1,1)
    bvar_obj$prior(prior, # prior mean value
                   1,     # var_type
                   1,     # decay_type
                   0.2,   # HP1
                   0.2,   # HP2
                   10^3,  # HP3
                   1.0)   # HP4
    bvar_obj$gibbs(10000)
    
    # Generate forecast for the next period
    forecast_result <- forecast(bvar_obj, shocks = TRUE, var_names = colnames(data)[-1], 
                                back_data = 5, save = FALSE, periods = 1)
    
    # Store the forecast in the list
    forecast_result <- forecast_result[1]
    forecast <- forecast_result$forecast_mean[1, 1]
    forecast_results <- rbind(forecast_results, forecast)
    
    actual_value <- data[(i + window_size - 1), "return"]
    actual_values <- rbind(actual_values, actual_value)
    
    date <- data[(i+ window_size - 1), "Date"]
    dates <- rbind(dates, date)
  }
  
  df_moving <- cbind(dates, forecast_results, actual_values)
  df_moving$Date <- as.Date(df_moving$Date, format = "%Y-%m-%d")
  colnames(df_moving) <- c("Date", "Forecast", "Actual")
  df_moving <- df_moving[-nrow(df_moving), ]
  
  n_obs <- nrow(df_moving)
  
  # Calculate the number of correctly predicted signs (positive, negative)
  correct_pos <- sum(df_moving$Forecast > 0 & df_moving$Actual > 0)
  correct_neg <- sum(df_moving$Forecast < 0 & df_moving$Actual < 0)
  
  # Calculate the ratio of correctly predicted signs (positive, negative)
  ratio_correct_pos <- correct_pos / n_obs
  ratio_correct_neg <- correct_neg / n_obs
  
  # Calculate other forecast performance metrics
  # Mean Absolute Error (MAE)
  mae <- mean(abs(df_moving$Forecast - df_moving$Actual))
  
  # Root Mean Squared Error (RMSE)
  rmse <- sqrt(mean((df_moving$Forecast - df_moving$Actual)^2))
  
  # Mean Absolute Percentage Error (MAPE)
  mape <- mean(abs((df_moving$Actual - df_moving$Forecast) / df_moving$Actual)) * 100
  
  # Directional Accuracy (DA)
  da <- (correct_pos + correct_neg) / n_obs
  
  # Mean Forecast Error (MFE)
  mfe <- mean(df_moving$Forecast - df_moving$Actual)
  
  # Mean Absolute Scaled Error (MASE) - Assuming no seasonal component
  mase <- mae / mean(abs(diff(df_moving$Actual)))
  
  # Compute DA and EP statistics
  DA <- compute_DA(df_moving$Actual, df_moving$Forecast)
  EP_statistic <- compute_EP_statistics(df_moving$Forecast, df_moving$Actual)
  
  # Plot actual vs forecast
  plot_actual_vs_forecast <- ggplot(df_moving, aes(x = Date)) +
    geom_line(aes(y = Actual, color = "Actual")) +
    geom_line(aes(y = Forecast, color = "Forecast")) +
    labs(title = paste("Actual vs Forecast"),
         y = "Value", color = "Series") +
    theme_minimal()
  
  return(list(Ratio_Positive = ratio_correct_pos,
              Ratio_Negative = ratio_correct_neg,
              MAE = mae,
              RMSE = rmse,
              MAPE = mape,
              Percent_correct = da,
              DA = DA,
              MFE = mfe,
              MASE = mase,
              EP_Statistic = EP_statistic,
              plot = plot_actual_vs_forecast
  ))
}





# Perform analysis for each subset
subset1_stats <- perform_analysis1(subset1, 60)
subset1_stats_plot <- subset1_stats[11]
subset1_stats <- subset1_stats[-11]
subset1_stats
subset2_stats <- perform_analysis1(subset2, 60)
subset2_stats_plot <- subset2_stats[11]
subset2_stats <- subset2_stats[-11]
subset2_stats
subset3_stats <- perform_analysis1(subset3, 60)
subset3_stats_plot <- subset3_stats[11]
subset3_stats <- subset3_stats[-11]
subset3_stats

whole_variables <- data.frame(whole_sample = subset3_stats, before2011 = subset1_stats , after2011 = subset2_stats)
whole_variables <- t(whole_variables)
whole_variables <- data.frame(
  Column1 = whole_variables[1:10, 1],
  Column2 = whole_variables[11:20, 1],
  column3 = whole_variables[21:30,1]
)
colnames(whole_variables) <- c("1991-2023", "1991-2011", "2011-2023")
rownames(whole_variables) <- c("Ratio_Positive", "Ratio_Negative", "MAE", "RMSE", "MAPE","Percent_correct", "DA", "MFE","MASE", "EP_Statistic")
whole_variables





# Function to perform analysis for different subsets
perform_analysis_subset <- function(data, window_size) {
  # Initialize a list to store results
  subset_results <- list()
  
  # Perform analysis for usa[,2:5]
  fin_only <- perform_analysis1(data[, c(1,2:5)], window_size)
  subset_results$fin_only <- fin_only
  
  # Perform analysis for usa[,2,6:8]
  macro_only <- perform_analysis1(data[, c(1,2,6:8)], window_size)
  subset_results$macro_only <- macro_only
  
  return(subset_results)
}




# Perform analysis for each subset for the whole sample
subset_results <- perform_analysis_subset(usa, 60)
print(subset_results$fin_only)
subset_results_fin_only_plot <- subset_results$fin_only[11]
subset_results$fin_only <- subset_results$fin_only[-11]
print(subset_results$macro_only)
subset_results_macro_only_plot <- subset_results$macro_only[11]
subset_results$macro_only <- subset_results$macro_only[-11]

wholesample <- data.frame(wholesample_fin = subset_results$fin_only, whole_Sample_macro = subset_results$macro_only)
wholesample <- t(wholesample)
wholesample <- data.frame(
  Column1 = wholesample[1:10, 1],
  Column2 = wholesample[11:20, 1]
)
colnames(wholesample) <- c("WholeSampleFinOnly", "WholeSampleMacroOnly")
rownames(wholesample) <- c("Ratio_Positive", "Ratio_Negative", "MAE", "RMSE", "MAPE","Percent_correct", "DA", "MFE","MASE", "EP_Statistic")
wholesample


# Perform analysis for each subset for subset1 (Before 2011)
subset_results_before2010 <- perform_analysis_subset(subset1, 60)
print(subset_results_before2010$fin_only)
print(subset_results_before2010$macro_only)

subset_results_before2010_fin_only_plot <- subset_results_before2010$fin_only[11]
subset_results_before2010$fin_only <- subset_results_before2010$fin_only[-11]
subset_results_before2010$fin_only

subset_results_before2010_macro_only_plot <- subset_results_before2010$macro_only[11]
subset_results_before2010$macro_only <- subset_results_before2010$macro_only[-11]
subset_results_before2010$macro_only
before2010 <- data.frame(before2010fin_only = subset_results_before2010$fin_only, before2010macro_only = subset_results_before2010$macro_only)
before2010 <- t(before2010)
before2010 <- data.frame(
  Column1 = before2010[1:10, 1],
  Column2 = before2010[11:20, 1]
)
colnames(before2010) <- c("FinOnly_before2010", "MacroOnly_before2010")
rownames(before2010) <- c("Ratio_Positive", "Ratio_Negative", "MAE", "RMSE", "MAPE","Percent_correct", "DA", "MFE","MASE", "EP_Statistic")
before2010


# Perform analysis for each subset for subset2 (After 2010)
subset_results_after2010 <- perform_analysis_subset(subset2, 60)
print(subset_results_after2010$fin_only)
print(subset_results_after2010$macro_only)

subset-results_after2010_fin_only_plot <- subset_results_after2010$fin_only[11]
subset_results_after2010$fin_only <- subset_results_after2010$fin_only[-11]
subset_results_after2010$fin_only

subset-results_after2010_macro_only_plot <- subset_results_after2010$macro_only[11]
subset_results_after2010$macro_only <- subset_results_after2010$macro_only[-11]
subset_results_after2010$macro_only
after2010 <- data.frame( after2010fin_only = subset_results_after2010$fin_only, after2010macro_only = subset_results_after2010$macro_only)
after2010 <- t(after2010)
after2010 <- data.frame(
  Column1 = after2010[1:10, 1],
  Column2 = after2010[11:20, 1]
)
colnames(after2010) <- c("FinOnly_after2010", "MacroOnly_after2010")
rownames(after2010) <- c("Ratio_Positive", "Ratio_Negative", "MAE", "RMSE", "MAPE","Percent_correct", "DA", "MFE","MASE", "EP_Statistic")
after2010


bütünperiyod <- cbind(whole_variables[1], wholesample[1], wholesample[2])
bütünperiyod
doksanbir11 <- cbind(whole_variables[2], before2010[1], before2010[2] )
doksanbir11
onbir23 <- cbind(whole_variables[3], after2010[1], after2010[2])
onbir23


library(xtable)

# Assume df1, df2, and df3 are your dataframes
# Create a LaTeX table for each dataframe
table1 <- xtable(bütünperiyod)
table2 <- xtable(doksanbir11)
table3 <- xtable(onbir23)

# Print LaTeX code for each table
print(table1, include.rownames = TRUE)
print(table2, include.rownames = TRUE)
print(table3, include.rownames = TRUE)

