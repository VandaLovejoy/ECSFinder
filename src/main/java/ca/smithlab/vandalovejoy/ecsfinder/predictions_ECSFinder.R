#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(caret)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript script.R <input_file> <output_file> <model_dir>")
}
input_file  <- args[1]
output_file <- args[2]
model_dir   <- args[3]
# Load the pretrained Random Forest model
model_path <- file.path(model_dir, "final_rf_model.rds")
if (!file.exists(model_path)) {
  stop("Model file not found at: ", model_path)
}
model_rf <- readRDS(model_path)

# Read and validate input data
test_data <- read.csv(input_file, stringsAsFactors = FALSE)
if (nrow(test_data) == 0) {
  stop("Input file is empty or improperly formatted.")
}

# Required columns
required_cols <- c(
  "name_file", "min_energy", "pseudo_energy", "log_min_evalue",
  "covarying_bp", "MPI", "average_MFE_sample", "sd_sample",
  "zscore", "sci"
)
missing_cols <- setdiff(required_cols, names(test_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Prepare features
x_test <- test_data[required_cols]
# Handle infinite and missing values
x_test$log_min_evalue[is.infinite(x_test$log_min_evalue)] <- -1e6
x_test$covarying_bp[is.na(x_test$covarying_bp)]       <- 0

# Prediction function
predict_probs <- function(model, newdata) {
  if ("train" %in% class(model) || "randomForest" %in% class(model)) {
    predict(model, newdata = newdata, type = "prob")[, 2]
  } else {
    stop("Unsupported model class: ", paste(class(model), collapse = ", "))
  }
}

# Generate predictions
predicted_probs <- predict_probs(model_rf, x_test)
threshold <- 0.431
predicted_class <- ifelse(predicted_probs >= threshold, "TP", "FP")

# Compile results
predictions <- data.frame(
  name_file               = test_data$name_file,
  Predicted_Probabilities = predicted_probs,
  Predicted_Class         = predicted_class,
  stringsAsFactors        = FALSE
)

# Write output (create or append)
write_header <- !file.exists(output_file)
write.table(
  predictions,
  file      = output_file,
  sep       = ",",
  row.names = FALSE,
  col.names = write_header,
  quote     = FALSE,
  append    = !write_header
)
message(if (write_header) "File created and data written." else "Data appended to existing file.")

# Display predictions
print(predictions)
