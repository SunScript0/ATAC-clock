import sys, os
import time
import shelve, pickle
import numpy as np
import pandas as pd
from sklearn.model_selection import LeavePGroupsOut, GridSearchCV
from sklearn.metrics import mean_squared_error, median_absolute_error
from scipy.stats import pearsonr

# Settings
n_jobs=10
inner_groups_out = 1

# Read command line arguments
path = sys.argv[1] # Base path of ncv
# path = "/scratch/fmorandi/ChromAcc-clock/clocks/parallel/2023-02-03_14-38_tpm"

start_time = time.time()

# Load common dataset
dataset = shelve.open(path + "/dataset", "r")
for k in dataset:
    globals()[k]=dataset[k]
    print("Loaded", k)
dataset.close() 

# Define inner CV
cv_inner = LeavePGroupsOut(inner_groups_out)

# Optimize hyperparameters by GridSearchCV
scoring = "neg_median_absolute_error"
gridCV = GridSearchCV(pipe, param_grid, scoring=scoring, n_jobs=n_jobs, cv=cv_inner, verbose=1, refit=True)
gridCV.fit(data, labels, groups=groups)

# Save gridcv object
gridcv_file = open(path + "/final_fitted_gridcv", "wb") 
pickle.dump(gridCV, gridcv_file)
gridcv_file.close()

final_coefs = pd.DataFrame({
	"mu": gridCV.best_estimator_["scaler"].mean_,
	"var": gridCV.best_estimator_["scaler"].var_,
	"coef": gridCV.best_estimator_["regressor"].coef_})
final_coefs.index = data.columns
final_coefs.loc["intercept", :] = [np.NaN, np.NaN, gridCV.best_estimator_["regressor"].intercept_]
final_coefs.to_csv(path + "/final_coefs.tsv", sep="\t")