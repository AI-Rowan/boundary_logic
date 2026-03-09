############################################################
# Section 2.1: Data Input Preparation
# Description:
#   Script to prepare the dataset for modeling. It sets up the number of
#   independent variables, ensures class columns are numeric, stores
#   variable names, and labels the columns consistently for downstream
#   analysis.
#
# Notes:
#   * This script assumes `data` is already loaded in the environment.
#   * removed multiclass "original class" variable in step 3
#   * Data engineering:
        # If you want to apply any data engineering, the procedure must be stored and reversible, to allow return to original variables
        # this is needed to interpret changes needed ito original variables understood by the user
        # engineering such as converting categorical to numeric or removing outliers etc can be done not to retain it all the data
        # need to apply it to future test data to be explained
        # Therefore, if standardisation is required, it is done later and process of reversing it is in the code
        # standardise <- T used later is for visualisation and distance calc, not to improve the model prediction accuracy needed in data engineering step

############################################################


#############################
# 1) Select dataset and variables to include: names "data"
#############################

# "2 Data File.R"

#############################
# 2) Set Number of Variables
#############################

# Number of independent variables in the dataset
# Example: PIMA = 8; Iris = 4; Projection example = 2
num_vars <- ncol(data) - 1

##############################################
# 3) Prepare Class Columns for Multiclass Use
##############################################

# For multiclass problems:
#   * Keep the original class in column (num_vars+2) # >> removed this feature, remove multiclass prediction feature
#   * Set the prediction class in column (num_vars+1)
#   * Convert class column to numeric
#   * Assign Original_class for later use
data_transform <- data
data_transform[, num_vars + 1] <- as.numeric(data_transform[, num_vars + 1])

#>> removed: store original class
# data_transform[, num_vars + 2] <- as.numeric(data_transform[, num_vars + 1]) 
# if LDA keep the multi class classification
# if(lda == F) data_transform[,num_vars+1] <- ifelse(data_transform[,num_vars+1]==tdp_class,1,0)
# final comparison will still be two class, so after fitting model, class and prediction will be made 0/1

## Now for all data types, including multiclass, ensure target variable set to 1
data_transform[,num_vars+1] <- ifelse(data_transform[,num_vars+1]==tdp_class,1,0)


#############################
# 4) Store Variable Names
#############################
# Capture the feature names for reference
var_names <- colnames(data_transform)[1:num_vars]

#############################
# 5) Rename Columns
#############################
# Assign consistent column names:
#   * Feature columns = var_names
#   * class = model class (numeric): Required for 2.2 where the model is fitted based on class ~ variables
#   * Original_class = retained true class
#>> removed: colnames(data_transform) <- c(var_names, "class", "Original_class")
colnames(data_transform) <- c(var_names, "class")

##########################################################
# 6) data pre-process and split in test and train
##########################################################

## first standardise the data: will need to convert final answer back 
if(first_standardise  == T)  {
  data_sd <- apply(data, 2, "sd")
  data_center <- apply(data, 2, "mean")
  data[,1:num_vars] <- scale(data[,1:num_vars], center = T, scale = T)
}

## split data into test and train
{
  set.seed(121)
  train_idx <- sample(nrow(data_transform), 4/5* nrow(data_transform))
  train_data <- data_transform[train_idx, ] #train
  test_data <- data_transform[-train_idx, ] #test: data to be explained
}

rm(data_transform)

#########################
####### END #############
#########################
