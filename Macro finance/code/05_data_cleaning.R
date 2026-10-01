# Author: Jean de Biolley #

# The following code cleans and standardizes historical monthly data from the Kenneth French
# 30 Industry Portfolios and Fama-French 3 Factors datasets, computing net monthly excess returns
# over the risk-free rate from July 1927 onward.


library(tidyverse)
library(stringr)
library(dplyr)

all_folders <- c(
  "30_Industry_Portfolios_CSV",
  "F-F_Research_Data_Factors_CSV"
)

process_french_data <- function(folder_name) {
  # Construct input and output file paths dynamically
  file_name <- str_replace(folder_name, "(?i)_csv$", ".csv")
  input_file_path <- file.path("raw data", folder_name, file_name)
  output_file_path <- file.path("input data", file_name)

  if (!file.exists(input_file_path)) {
    warning(paste("File not found, skipping:", input_file_path))
    return(NULL)
  }

  # Set date boundaries and regular expressions based on sampling frequency
  is_daily <- str_detect(tolower(folder_name), "daily")

  if (is_daily) {
    start_date <- 19700101
    end_date <- 20251231
    date_regex <- "^\\d{8}$"
  } else {
    start_date <- 192707
    end_date <- 202603
    date_regex <- "^\\d{6}$"
  }

  is_factor_file <- str_detect(folder_name, "Factors")

  # Load the raw text lines from the data file
  all_lines <- readLines(input_file_path, warn = FALSE)

  if (is_factor_file) {
    # Identify the header line using character matching for factor labels
    header_idx <- which(
      str_detect(all_lines, "Mkt-RF") &
        str_detect(all_lines, "SMB") &
        str_detect(all_lines, "HML")
    )
    if (length(header_idx) == 0) {
      warning(paste("Anchor not found in:", file_name, "- Skipping."))
      return(NULL)
    }
    anchor_idx <- header_idx[1]
    start_idx <- anchor_idx
  } else {
    # Identify the target table segment label for industry structures
    anchor_string <- if (is_daily) {
      "Average Value Weighted Returns -- Daily"
    } else {
      "Average Value Weighted Returns -- Monthly"
    }

    anchor_idx <- which(str_detect(all_lines, fixed(anchor_string)))
    if (length(anchor_idx) == 0) {
      warning(paste("Anchor not found in:", file_name, "- Skipping."))
      return(NULL)
    }
    anchor_idx <- anchor_idx[1]
    start_idx <- anchor_idx + 1
  }

  # Locate the terminal blank line separating the sub-table matrices
  is_blank <- str_trim(all_lines) == ""
  blank_lines_after_anchor <- which(is_blank & seq_along(all_lines) > start_idx)

  if (length(blank_lines_after_anchor) == 0) {
    end_idx <- length(all_lines)
  } else {
    end_idx <- min(blank_lines_after_anchor) - 1
  }

  # Slice the raw text to encapsulate the specific dataset window
  table_lines <- all_lines[start_idx:end_idx]

  # Parse the subset lines and reformat values to decimals
  df <- read_csv(paste(table_lines, collapse = "\n"),
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  )
  clean_df <- df %>%
    rename(date = 1) %>%
    mutate(date = str_trim(as.character(date))) %>%
    filter(str_detect(date, date_regex)) %>%
    mutate(date = as.integer(date)) %>%
    filter(date >= start_date & date <= end_date) %>%
    mutate(across(-date, ~ {
      val <- as.numeric(.)
      val[val == -99.99 | val == -999] <- 0
      return(val / 100)
    }))

  # Create the export directory if missing and write out the clean dataset
  dir.create(dirname(output_file_path), showWarnings = FALSE, recursive = TRUE)
  write_csv(clean_df, output_file_path)

  cat("Successfully processed and saved:", file_name, "\n")
}

# Run the batch transformation routine over target directories
walk(all_folders, process_french_data)

# Import the formatted CSV records
X30_Industry_Portfolios <- read_csv("input data/30_Industry_Portfolios.csv", show_col_types = FALSE)
FF_Research_Data_Factors <- read_csv("input data/F-F_Research_Data_Factors.csv", show_col_types = FALSE)

# Compute portfolio returns relative to the risk-free asset rate
calculate_excess <- function(portfolio_df, factors_df) {
  # Isolate the risk-free rate vector from the asset factor structure
  rf_data <- factors_df %>%
    dplyr::select(date, matches("(?i)^rf$")) %>%
    dplyr::rename(rf = 2)

  # Transform columns by mapping matches and subtracting the risk-free benchmarks
  portfolio_df %>%
    left_join(rf_data, by = "date") %>%
    mutate(across(-c(date, rf), ~ . - rf)) %>%
    select(-rf)
}

# Construct the excess asset matrix structure
X30_monthly_r_excess <- calculate_excess(X30_Industry_Portfolios, FF_Research_Data_Factors)

rm(calculate_excess)
