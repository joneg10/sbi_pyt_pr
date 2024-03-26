from create_training_set_class import TrainingSet
import pandas as pd

###############################################
### Creation of formatting training dataset ###
###############################################

# Select path where PDB files are
paths = TrainingSet("../scpdb_files_2000")

# Get the training set
training_set = paths.get_formated_set()

# Fill NaN values with 0

training_set.fillna(0, inplace=True)


# Set "is_lbs" as the last column

training_set["is_lbs"] = training_set.pop("is_lbs")

# Save the training set to a parquet file
output_file = '../trainingSet_20240325.parquet'
 
training_set.to_parquet(output_file)
