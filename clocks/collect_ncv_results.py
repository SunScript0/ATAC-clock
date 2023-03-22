from clock_utils import get_ncv_results, plot_ncv_results
import sys

path = sys.argv[1]

summary, preds, coefs = get_ncv_results(path, 11, True)
plot_ncv_results(preds, path)