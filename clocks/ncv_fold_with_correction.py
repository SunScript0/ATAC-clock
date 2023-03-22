import sys, os
import time
import shelve, pickle
import numpy as np
import pandas as pd
from sklearn.model_selection import LeavePGroupsOut, GridSearchCV
from sklearn.metrics import mean_squared_error, median_absolute_error
from scipy.stats import pearsonr
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import LinearRegression, RANSACRegressor, HuberRegressor

# Settings
n_jobs=10
inner_groups_out = 1

# Read command line arguments
path = sys.argv[1] # Base path of ncv
i = int(sys.argv[2]) # Fold index
# path = "/scratch/fmorandi/ChromAcc-clock/clocks/parallel/2023-01-28_13-48"
# i=0

# Create directory for outputs
outpath = path + "/fold_" + str(i)
os.makedirs(outpath, exist_ok=True)

start_time = time.time()

log = open(outpath + "/.log", "w")
sys.stdout = log
print(path)
print("Fold =", i)

# Load common dataset
dataset = shelve.open(path + "/dataset", "r")
for k in dataset:
    globals()[k]=dataset[k]
    print("Loaded", k)
dataset.close() 

# Outer partition
train_data = data.iloc[outer_split[i][0], :]
train_ccomp = data_ccomp_only.iloc[outer_split[i][0], :]
train_labels = labels.iloc[outer_split[i][0]]
train_groups = groups[outer_split[i][0]]
test_data = data.iloc[outer_split[i][1], :]
test_ccomp = data_ccomp_only.iloc[outer_split[i][1], :]
test_labels = labels.iloc[outer_split[i][1]]
test_groups = groups[outer_split[i][1]]

# Correct for cell composition
if ccomp_regressor == "RANSAC":
    mlm = MultiOutputRegressor(RANSACRegressor())
elif ccomp_regressor == "Huber":
    mlm = MultiOutputRegressor(HuberRegressor(max_iter=500))
elif ccomp_regressor == "lm":
    mlm = MultiOutputRegressor(LinearRegression())
mlm.fit(train_ccomp.join(train_labels), train_data)
crct_coefs = pd.DataFrame(index = data.columns, columns= np.append(data_ccomp_only.columns.values, "Age"))
for ind, lm in enumerate(mlm.estimators_):
    if ccomp_regressor == "RANSAC":
        crct_coefs.iloc[ind, :] = lm.estimator_.coef_
    else:
        crct_coefs.iloc[ind, :] = lm.coef_
train_ccomp_effect = train_ccomp.values @ crct_coefs.iloc[:, 0:-1].T.values
test_ccomp_effect = test_ccomp.values @ crct_coefs.iloc[:, 0:-1].T.values
train_data = train_data - train_ccomp_effect
test_data = test_data - test_ccomp_effect
crct_coefs.to_csv(outpath + "/ccomp_crct_coefs.tsv", sep="\t")

# Define inner CV
cv_inner = LeavePGroupsOut(inner_groups_out)

# Optimize hyperparameters by GridSearchCV
gridCV = GridSearchCV(pipe, param_grid, scoring=scoring, n_jobs=n_jobs, cv=cv_inner, verbose=1, refit=True)
gridCV.fit(train_data, train_labels, groups=train_groups)

# Predictions
preds = gridCV.predict(test_data)
mse = mean_squared_error(test_labels, preds)
mae = median_absolute_error(test_labels, preds)
r, _ = pearsonr(test_labels, preds)

# Save predictions
out = test_labels.to_frame()
out["Preds"] = preds
out["Fold"] = i
out.to_csv(outpath + "/preds.tsv", sep="\t")

# Save gridcv object
gridcv_file = open(outpath + "/fitted_gridcv", "wb") 
pickle.dump(gridCV, gridcv_file)
gridcv_file.close()

# Print performance to log
print("MSE =", mse, sep="\t")
print("MAE =", mae, sep="\t")
print("r =", r, sep="\t")
print("Finished in", time.time() - start_time)

log.close()