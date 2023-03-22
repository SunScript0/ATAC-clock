import numpy as np
import pandas as pd
import time
import subprocess
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.covariance import EllipticEnvelope
from sklearn.model_selection import KFold, GroupKFold, LeavePGroupsOut, GridSearchCV
from sklearn.linear_model import ElasticNet, ElasticNetCV, LinearRegression
from sklearn.metrics import mean_squared_error, median_absolute_error, r2_score
from scipy.stats import normaltest, kurtosis, skew, pearsonr, spearmanr
from sklearn.feature_selection import SelectFdr, SelectKBest
import warnings
import multiprocessing as mp
import pickle

        
def normalize(df, method, log_transform, row_info=None):
	# Normalizes count matrices with various methods and with/out log transform
	# df is a pandas dataframe with peaks as rows, samples as columns
    # pass row info with column for feature length ("Length") if needed by the norm method
    if method == "cpm":
        sample_counts = np.sum(df, axis = 0)
        df_cpm = df.copy()
        df_cpm = 1e6 * df / sample_counts
        if log_transform:
            df_cpm = np.log1p(df_cpm)
        return df_cpm
    if method == "rpkm":
        sample_counts = np.sum(df, axis = 0)
        df_rpkm = df.copy()
        df_rpkm = 1e6 * df / sample_counts
        try:
            peak_lengths = row_info["Length"]
        except:
            print("Need to provide row_info with feature lengths")
            return
        df_rpkm = 1000 * (df_rpkm.T / peak_lengths).T
        if log_transform:
            df_rpkm = np.log1p(df_rpkm)
        return df_rpkm
    if method == "tpm":
        try:
            peak_lengths = row_info["Length"]
        except:
            print("Need to provide row_info with feature lengths")
            return
        df_tpm = df.copy()
        df_tpm = 1000 * (df.T / peak_lengths).T
        sample_counts = np.sum(df_tpm, axis = 0)
        df_tpm = 1e6* df_tpm / sample_counts
        if log_transform:
            df_tpm = np.log1p(df_tpm)
        return df_tpm
    if method == "like deseq":
        df_deseq = df.copy()
        df_log = df_deseq.apply(np.log)
        pseudo_ref = df_log.mean(axis=1)
        good_ind = pseudo_ref.index[~np.isneginf(pseudo_ref)]
        ratio_to_ref = df_log.loc[good_ind].sub(pseudo_ref[good_ind], axis=0)
        correction_factors = ratio_to_ref.median(axis=0).apply(np.exp)
        # plt.hist(correction_factors)
        df_deseq = df_deseq.divide(correction_factors, axis='columns')
        if log_transform:
            df_deseq = np.log1p(df_deseq)
        return df_deseq
    if method == "like edger":
        df.to_csv("tmp_edger_counts.csv")
        subprocess.run(["Rscript", "tmm_edgeR.R", "tmp_edger_counts.csv", "tmm_edger_cpm.csv"])
        tmm_norm = pd.read_csv("tmm_edger_cpm.csv")
        df_tmm = df.copy()
        df_tmm = tmm_norm
        return df_tmm
    else:
        raise ValueError("Normalization method unknown")
        
def pca_outlier_detection(data, snames, ncomps, plot_boundary, cont):
	# Detect outliers by elliptic envelope on PCA
	# Optionally plots if using less than 2 PCs
	# data is numpy/pandas matrix with samples as rows
	# snames is sample names
	# ncomps is the number of components used for the PCA
	# plot_boundary is to produce a plot (2 PCs only)
	# cont is the fraction of all samples expected to be outliers
    plt_points = 100
    plt_margin = 10
    
    if ncomps > 2 and plot_boundary:
        raise ValueError("Cannot plot boundary when using more than 2 components")
        
    pca = PCA(n_components=ncomps)
    data_PCs = pca.fit_transform(data)
    #print("Explained variance ratio: " + str(pca.explained_variance_ratio_))
    #print("Cumulative explained variance ratio: " + str(np.sum(pca.explained_variance_ratio_)))

    cov = EllipticEnvelope(contamination=cont).fit(data_PCs)
    outliers = snames[cov.predict(data_PCs) < 0]
    
    if plot_boundary:
        fig, axs = plt.subplots(1, 1, figsize=(5, 5))
        axs.scatter(data_PCs[:, 0], data_PCs[:, 1], s=10)
        
        minX = np.min(data_PCs[:, 0]) - plt_margin
        maxX = np.max(data_PCs[:, 0]) + plt_margin
        minY = np.min(data_PCs[:, 1]) - plt_margin
        maxY = np.max(data_PCs[:, 1]) + plt_margin

        xx1, yy1 = np.meshgrid(np.linspace(minX, maxX, plt_points), np.linspace(minY, maxY, plt_points))
        Z1 = cov.decision_function(np.c_[xx1.ravel(), yy1.ravel()])
        Z1 = Z1.reshape(xx1.shape)
        plt.contour(xx1, yy1, Z1, levels=[0], linewidths=2)

        for sname, x, y in zip(snames, data_PCs[:, 0], data_PCs[:, 1]):
            if sname in outliers:
                axs.annotate(sname, (x, y), fontsize=10, ha="center")
    
    return outliers

def nested_cv(data, labels, pipe, param_grid, opt_choice_fun):
    start_time = time.time()
    # Hardcoded settings
    n_jobs=10
    outer_folds = 10
    inner_folds = 5
    scoring = "neg_median_absolute_error"
    # Safety check for data and label shapes
    if data.shape[0] != labels.shape[0]:
        raise ValueError("Shapes of data and labels don't match")
    # Create outer CV object
    cv_outer = KFold(n_splits=outer_folds, shuffle=True)
    # Create dictionary to store results
    outer_results = dict()
    outer_results["cv_results"] = list()
    outer_results["best_estimator"] = list()
    outer_results["best_params"] = list()
    outer_results["best_index"] = list()
    outer_results["outer_MSE"] = list()
    outer_results["outer_MAE"] = list()
    outer_results["outer_r"] = list()
    outer_results["predictions"] = np.empty((3, 0))
    outer_results["run_time"] = 0
    outer_fold_n=0
    for train_ind, test_ind in cv_outer.split(data):
        print("Starting fold " + str(outer_fold_n))
        # Define outer train and test sets
        if type(data) is np.ndarray:
            train_data = data[train_ind, :]
            train_labels = labels[train_ind]
            test_data = data[test_ind, :]
            test_labels = labels[test_ind]
        else:
            train_data = data.iloc[train_ind, :]
            train_labels = labels.iloc[train_ind]
            test_data = data.iloc[test_ind, :]
            test_labels = labels.iloc[test_ind]
        print(data.shape)
        print(train_data.shape)
        print(test_data.shape)
        # Define inner CV
        cv_inner = KFold(n_splits=inner_folds, shuffle=True)
        # optimize hyperparameters by GridSearchCV
#         gridCV = GridSearchCV(pipe, param_grid, scoring="neg_mean_squared_error", n_jobs=n_jobs, 
#                               cv=cv_inner, verbose=False, refit=opt_choice_fun)
        gridCV = GridSearchCV(pipe, param_grid, scoring=scoring, n_jobs=n_jobs, 
                              cv=cv_inner, verbose=False, refit=opt_choice_fun)
        gridCV.fit(train_data, train_labels)
        # predictions
        preds = gridCV.predict(test_data)
        mse = mean_squared_error(test_labels, preds)
        mae = median_absolute_error(test_labels, preds)
        r = np.sqrt(r2_score(test_labels, preds))
        # store results for this fold
        outer_results["predictions"] = np.hstack((outer_results["predictions"], 
                                                  np.vstack(
                                                      (np.tile(outer_fold_n, len(test_labels)), 
                                                       np.squeeze(test_labels), 
                                                       preds))))
        # performance metrics
        outer_results["outer_MSE"].append(mse)
        outer_results["outer_MAE"].append(mae)
        outer_results["outer_r"].append(r)
        outer_results["cv_results"].append(gridCV.cv_results_)
        outer_results["best_estimator"].append(gridCV.best_estimator_)
        outer_results["best_params"].append(gridCV.best_params_)
        outer_results["best_index"].append(gridCV.best_index_)
        outer_fold_n+=1
    outer_results["run_time"] = time.time() - start_time
    return outer_results

def group_nested_cv(data, labels, groups, pipe, param_grid):
    #warnings.filterwarnings("error")
    start_time = time.time()
    # Hardcoded settings
    n_jobs=10
    outer_groups_out = 1
    inner_groups_out = 1
    scoring = "neg_median_absolute_error"
    # Safety check for data and label shapes
    if data.shape[0] != labels.shape[0]:
        raise ValueError("Shapes of data and labels don't match")
    # Create outer CV object
    cv_outer = LeavePGroupsOut(outer_groups_out)
    # Create dictionary to store results
    outer_results = dict()
    outer_results["best_estimator"] = list()
    outer_results["best_params"] = list()
    outer_results["outer_MSE"] = list()
    outer_results["outer_MAE"] = list()
    outer_results["outer_r"] = list()
    outer_results["predictions"] = np.empty((3, 0))
    outer_results["run_time"] = 0
    outer_fold_n=0
    for train_ind, test_ind in cv_outer.split(data, labels, groups):
        print("Starting fold " + str(outer_fold_n))
        # Define outer train and test sets
        train_data = data.iloc[train_ind, :]
        train_labels = labels.iloc[train_ind]
        train_groups = groups[train_ind]

        test_data = data.iloc[test_ind, :]
        test_labels = labels.iloc[test_ind]
        test_groups = groups[test_ind]
        print("Total: {} | Train: {} | Test: {}".format(data.shape, train_data.shape, test_data.shape))
        # Define inner CV
        cv_inner = LeavePGroupsOut(inner_groups_out)
        # optimize hyperparameters by GridSearchCV
        gridCV = GridSearchCV(pipe, param_grid, scoring=scoring, n_jobs=n_jobs, 
                              cv=cv_inner, verbose=False)
        gridCV.fit(train_data, train_labels, groups=train_groups)
        # predictions
        preds = gridCV.predict(test_data)
        mse = mean_squared_error(test_labels, preds)
        mae = median_absolute_error(test_labels, preds)
        r, _ = pearsonr(test_labels, preds)
        # store results for this fold
        outer_results["predictions"] = np.hstack((outer_results["predictions"], 
                                                  np.vstack(
                                                      (np.tile(outer_fold_n, len(test_labels)), 
                                                       np.squeeze(test_labels), 
                                                       preds))))
        # performance metrics
        outer_results["outer_MSE"].append(mse)
        outer_results["outer_MAE"].append(mae)
        outer_results["outer_r"].append(r)
        outer_results["best_estimator"].append(gridCV.best_estimator_)
        outer_results["best_params"].append(gridCV.best_params_)
        outer_fold_n+=1
    outer_results["run_time"] = time.time() - start_time
    return outer_results

# def group_cv_parallel(outer_fold_n, train_data, train_labels, test_data, test_labels, train_groups, pipe, param_grid):
def group_cv_parallel(datapath, metapath):
    
    
    
    # Hardcoded settings
    n_jobs=1
    inner_groups_out = 1
    scoring = "neg_median_absolute_error"
    # Define inner CV
    cv_inner = LeavePGroupsOut(inner_groups_out)
    # Optimize hyperparameters by GridSearchCV
    gridCV = GridSearchCV(pipe, param_grid, scoring=scoring, n_jobs=n_jobs, cv=cv_inner, verbose=False, refit=True)
    gridCV.fit(train_data, train_labels, groups=train_groups)
    # Predictions
    preds = gridCV.predict(test_data)
    mse = mean_squared_error(test_labels, preds)
    mae = median_absolute_error(test_labels, preds)
    r, _ = pearsonr(test_labels, preds)
    return (outer_fold_n, preds, mse, mae, r)

def group_nested_cv_regpath(data, labels, groups, pipe, param_grid):
    #warnings.filterwarnings("error")
    start_time = time.time()
    # Hardcoded settings
    n_jobs=10
    outer_groups_out = 1
    inner_groups_out = 1
#     scoring = "neg_median_absolute_error" # ElasticNetCV doesnt take this argument
    # Safety check for data and label shapes
    if data.shape[0] != labels.shape[0]:
        raise ValueError("Shapes of data and labels don't match")
    # If data is a pandas df convert to numpy
    if type(data) is not np.ndarray:
        data = data.to_numpy()
    # Create outer CV object
    cv_outer = LeavePGroupsOut(outer_groups_out)
    # Create dictionary to store results
    outer_results = dict()
    outer_results["best_estimator"] = list()
    outer_results["best_params"] = list()
    outer_results["outer_MSE"] = list()
    outer_results["outer_MAE"] = list()
    outer_results["outer_r"] = list()
    outer_results["predictions"] = np.empty((3, 0))
    outer_results["run_time"] = 0
    outer_fold_n=0
    # Loop over outer folds
    for train_ind, test_ind in cv_outer.split(data, labels, groups):
        print("Starting fold " + str(outer_fold_n))
        # Define outer train and test sets
        train_data = data[train_ind, :]
        train_labels = labels[train_ind]
        train_groups = groups[train_ind]

        test_data = data[test_ind, :]
        test_labels = labels[test_ind]
        test_groups = groups[test_ind]
        print("Total: {} | Train: {} | Test: {}".format(data.shape, train_data.shape, test_data.shape))
        # Define inner CV
#         cv_inner = LeavePGroupsOut(inner_groups_out) cant do this in ElasticNetCV so will do simple 10 fold
        # Get hyperparameter space
        alphas = param_grid["regressor__alpha"]
        l1_ratios = param_grid["regressor__l1_ratio"]
        # Optimize hyperparameters by ElasticNetCV
        enetCV = ElasticNetCV(l1_ratio=l1_ratios, alphas=alphas, 
                              normalize=True, max_iter=5000, tol=0.0005,
                              cv=10, n_jobs=n_jobs, verbose=False)
        enetCV.fit(train_data, train_labels)
        # predictions
        preds = enetCV.predict(test_data)
        mse = mean_squared_error(test_labels, preds)
        mae = median_absolute_error(test_labels, preds)
#         r = np.sqrt(np.abs(r2_score(test_labels, preds)))
        r, _ = pearsonr(test_labels, preds)
        # store results for this fold
        outer_results["predictions"] = np.hstack((outer_results["predictions"], 
                                                  np.vstack(
                                                      (np.tile(outer_fold_n, len(test_labels)), 
                                                       np.squeeze(test_labels), 
                                                       preds))))
        # performance metrics
        outer_results["outer_MSE"].append(mse)
        outer_results["outer_MAE"].append(mae)
        outer_results["outer_r"].append(r)
        outer_results["best_estimator"].append({"regressor": enetCV})
        outer_results["best_params"].append({"alpha": enetCV.alpha_, "l1_ratio": enetCV.l1_ratio_})
        outer_fold_n+=1
    outer_results["run_time"] = time.time() - start_time
    return outer_results

def make_groups(ngroups, labels):
    nsamp = len(labels)
    labels_sorted_ind = labels.argsort()

    groups = np.zeros(nsamp)
    for i in range(nsamp):
        s_ind = labels_sorted_ind[i]
        groups[s_ind] = i % ngroups
    return groups

def plot_ncv_results(preds, savepath=None):
    fig, axs = plt.subplots(1, 1, figsize=(8, 6))
    sc = axs.scatter(preds["age"], preds["pred"], c=preds["fold"], cmap="tab10")
    plt.colorbar(sc)
    xrange = plt.xlim()
    yrange = plt.ylim()
    axs.plot([0, 100], [0, 100])
    axs.fill_between([0, 100], [-5, 95], [5, 105], alpha=0.1)
    plt.xlim(xrange)
    plt.ylim(yrange)
    if (savepath is not None):
        plt.savefig(savepath + "/preds.png")
    err = preds["age"] - preds["pred"]
    print("RMSE: " + str(np.sqrt(np.mean(err **2))))
    print("MAE: " + str(np.median(np.abs(err))))
    print("r: " + str(pearsonr(preds["age"], preds["pred"])[0]))
    
def get_ncv_results(path, nfold, save=False):
    preds = pd.read_csv(path + "/preds.tsv", sep="\t", index_col=0, names = ["age", "pred", "fold"])
    summary = pd.DataFrame()
    coefs = list()
    for i in range(nfold): 
        gridcv_file = open(path + "/fold_{}/fitted_gridcv".format(i), "rb") 
        gridCV = pickle.load(gridcv_file)
        gridcv_file.close()
        summary.loc[i, "log10(alpha)"] = np.log10(gridCV.best_params_["regressor__alpha"])
        summary.loc[i, "l1_ratio"] = gridCV.best_params_["regressor__l1_ratio"]
        coefs.append(gridCV.best_estimator_["regressor"].coef_)
        summary.loc[i, "non_zero"] = np.sum(coefs[i] != 0)
        sset = preds[preds["fold"] == i]
        summary.loc[i, "RMSE"] = np.sqrt(mean_squared_error(sset["age"], sset["pred"]))
        summary.loc[i, "MAE"] = median_absolute_error(sset["age"], sset["pred"])
        summary.loc[i, "r"], _ = pearsonr(sset["age"], sset["pred"])
    coefs = pd.DataFrame(coefs)
    if save:
        summary.to_csv(path + "/summary.tsv", sep="\t")
        coefs.to_csv(path + "/ncv_coefs.tsv", sep="\t")
    return summary, preds, pd.DataFrame(coefs)

def get_ncv_regpath_results(path, nfold, save=False):
    preds = pd.read_csv(path + "/preds.tsv", sep="\t", index_col=0, names = ["age", "pred", "fold"])
    summary = pd.DataFrame()
    coefs = list()
    for i in range(nfold): 
        pipe_file = open(path + "/fold_{}/fitted_pipe".format(i), "rb") 
        pipe = pickle.load(pipe_file)
        pipe_file.close()
        summary.loc[i, "log10(alpha)"] = np.log10(pipe["regressor"].alpha_)
        summary.loc[i, "l1_ratio"] = pipe["regressor"].l1_ratio_
        coefs.append(pipe["regressor"].coef_)
        summary.loc[i, "non_zero"] = np.sum(coefs[i] != 0)
        sset = preds[preds["fold"] == i]
        summary.loc[i, "RMSE"] = np.sqrt(mean_squared_error(sset["age"], sset["pred"]))
        summary.loc[i, "MAE"] = median_absolute_error(sset["age"], sset["pred"])
        summary.loc[i, "r"], _ = pearsonr(sset["age"], sset["pred"])
    coefs = pd.DataFrame(coefs)
    if save:
        summary.to_csv(path + "/summary.tsv", sep="\t")
        coefs.to_csv(path + "/ncv_coefs.tsv", sep="\t")
    return summary, preds, pd.DataFrame(coefs)
    