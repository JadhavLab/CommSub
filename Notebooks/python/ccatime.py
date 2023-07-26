import numpy as np
from jax import numpy as jnp
import cupy as cp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from scipy.stats import f
from sklearn.model_selection import TimeSeriesSplit, KFold
import time, os
from tqdm import tqdm
tqdm.pandas()
from scipy.stats import f
import multiprocessing, dill
folder = '/Volumes/MATLAB-Drive/Shared/figures/tables/'
name   = 'ZT2powerccatime'
def granger_causality_jax(x, y, lag_order):
    # Compute the lagged versions of x and y
    X = jnp.column_stack([x[t-lag_order:t] for t in range(lag_order, len(x))])
    Y = y[None, lag_order:]
    X = X.T
    Y = Y.T
    # Estimate the autoregressive models
    model_xy = jnp.linalg.lstsq(X, Y, rcond=None)[0]
    model_y  = jnp.linalg.lstsq(X[:, 1:], Y, rcond=None)[0]
    # Compute the residual sum of squares
    rss_xy = jnp.sum((Y - jnp.dot(X, model_xy))**2)
    rss_y  = jnp.sum((Y - jnp.dot(X[:, 1:], model_y))**2)
    # Compute the F-statistic
    f_stat = ((rss_y - rss_xy) / lag_order) / (rss_xy / (len(Y) - 2 * lag_order))
    # Compute the P-value
    p_value = 1 - jnp.sum(np.random.f(1, lag_order, len(Y) - 2 * lag_order) > f_stat) / len(Y)
    return f_stat, p_value
# Define Granger causality function
def granger_causality_cupy(x, y, lag_order, xtest=None, ytest=None):
    x = cp.asarray(x)
    y = cp.asarray(y)
    # Compute the lagged versions of x and y
    X = cp.column_stack([x[t-lag_order:t] for t in range(lag_order, len(x))])
    Y = y[None, lag_order:]
    X = X.T
    Y = Y.T
    # Estimate the autoregressive models
    model_xy = cp.linalg.lstsq(X, Y, rcond=None)[0]
    model_y = cp.linalg.lstsq(X[:, 1:], Y, rcond=None)[0]
    # Compute the residual sum of squares
    if xtest is None:
        yhat_xy = cp.dot(X, model_xy)
        yhat_y  = cp.dot(X[:, 1:], model_y)
    else:
        xtest = cp.asarray(xtest)
        xtest = cp.column_stack([xtest[t-lag_order:t] for t in range(lag_order, len(xtest))])
        xtest = xtest.T
        yhat_xy = cp.dot(xtest, model_xy)
        yhat_y  = cp.dot(xtest[:, 1:], model_y)
    if ytest is None:
        ytest = Y
    else:
        ytest = cp.asarray(ytest)
        ytest = ytest[None, lag_order:]
        ytest = ytest.T
    total   = cp.sum((ytest - cp.mean(ytest))**2)
    rss_xy  = cp.sum((ytest - yhat_xy)**2)
    rss_y   = cp.sum((ytest - yhat_y)**2)
    n = len(ytest)
    sigma_sq = rss_xy / (n - 2 * lag_order)
    LR = (rss_y - rss_xy) / sigma_sq
    # Compute the F-statistic
    f_stat = ((rss_y - rss_xy) / lag_order) / (rss_xy / (len(ytest) - 2 * lag_order))
    # Compute the P-value
    df1 = lag_order
    df2 = len(ytest) - 2 * lag_order
    p_value = 1 - f.cdf(f_stat.get(), df1, df2)  # Use scipy's F-distribution CDF function
    r2_xy = 1 - rss_xy / total # Compute r2 for x->y
    r2_y  = 1 - rss_y / total # Compute r2 for y
    # r2_xyfrac = r2_xy / (r2_xy + r2_y) # Compute fraction of variance explained by x->y
    r2_xy_part = 1 - (rss_y - rss_xy) / total # Compute fraction of variance explained by x->y
    # Compute yhat
    return f_stat.get(), p_value, r2_xy.get(), r2_y.get(), r2_xy_part.get(), LR.get()
# Define directionality test function
def directionality_test(granger_func, x, y, lag_order, xtest=None, ytest=None):
    f_stat_xy, p_value_xy, r2_xy, r2_y, r2_xy_p = granger_func(x, y, lag_order, xtest, ytest)
    f_stat_yx, p_value_yx, r2_yx, r2_x, r2_yx_p = granger_func(y, x, lag_order, xtest, ytest)
    direction = "x->y" if f_stat_xy > f_stat_yx else "y->x"
    significance = "significant" if min(p_value_xy, p_value_yx) < 0.05 else "not significant"
    magnitude = abs(r2_xy - r2_yx)
    magnitude_p = abs(r2_xy_p - r2_yx_p)
    return direction, significance, magnitude, magnitude_p
def is_continuos_var(df, col):
    return df[col].nunique() > 20 
def get_cont_vars(df):
    return [col for col in df.columns if is_continuos_var(df, col)]



# Example usage
x = np.random.randn(100)
y = 0.5 * x + np.random.randn(100)
lag_order = 5
start = time.time()
results  = granger_causality_cupy(x, y, lag_order)
end = time.time()
result2,p2 = granger_causality_jax(x, y, lag_order)
end2 = time.time()
result3 = grangercausalitytests(np.column_stack([x,y]), lag_order, verbose=False)
end3 = time.time()
print(result)
print(result2)
print(result3)
print("cupy")
print("time: ", end-start)
print("jax")
print("time: ", end2-end)
print("statsmodels")
print("time: ", end3-end2)

# Load the data
df = pd.read_csv(os.path.join(folder, f'{name}.csv'))
cont_vars = get_cont_vars(df)
dfc = df[cont_vars]
# zscore dfc columns
dfc = dfc.apply(lambda x: (x - x.mean()) / x.std())
# Define the lag
maxlag = 10
method = 'cupy'
gmeth = granger_causality_jax if method == 'jax' else granger_causality_cupy
# Initialize the KFold object
DF_results = []
n_splits = 5
tscv = TimeSeriesSplit(n_splits=n_splits)

# MULTI-PROCESSING
from multiprocessing import cpu_count
chunk_size = 240_000 // cpu_count()
num_chunks = df.shape[0] // chunk_size
print("num_chunks: ", num_chunks)
# Define a function to process a pair of columns
def process_pair(pair):
    with cp.cuda.Device():
        i, j = pair
        if i == j:
            return None
        print(f"Processing {i}, {j}")
        # Extract the time series for the pair of columns
        X = dfc.iloc[:, [i, j]]
        X = X.dropna().values
        # num_chunks = X.shape[0] // chunk_size
        RESULTS = []
        for chunk_index in tqdm(range(num_chunks), desc='Chunk', total=num_chunks):
            start_index = chunk_index * chunk_size
            end_index = min(start_index + chunk_size, X.shape[0])
            X_chunk = X[start_index:end_index]
            cp._default_memory_pool.free_all_blocks()
            X_chunk = cp.asarray(X_chunk, dtype=cp.float32)
            # Initialize a list to store the results for each split
            results = []
            s, f = 0, 0
            for itest, (train_index, test_index) in enumerate(tscv.split(X_chunk)):
                # Split the data into training and test sets
                X_train, X_test = X_chunk[train_index], X_chunk[test_index]
                if method == 'cupy' or method == 'jax':
                    for lag in range(1, maxlag+1):
                        # try:
                            tup = gmeth(X_train[:, 0], X_train[:, 1], lag+1,
                                        xtest=X_test[:, 0], ytest=X_test[:, 1])
                            results.append(dict(
                                lag=lag,
                                pvalue=tup[1],
                                F=tup[0],
                                r2_xy=tup[2],
                                r2_y=tup[3],
                                r2_xy_p=tup[4],
                                LR=tup[5],
                                itest=itest
                            ))
                            s += 1
                        # except Exception as e:
                        #     print(e)
                        #     results.append(None)
                        #     f += 1
                else:
                    t1=time.time()
                    result = grangercausalitytests(X_train, maxlag=maxlag, verbose=False)
                    t2=time.time()
                    print("time: ", t2-t1)
                    pvalue = result[maxlag][0]['ssr_ftest'][1]
                    F = result[maxlag][0]['ssr_ftest'][0]
                    results.append([list(range(maxlag)), pvalue, F])
                    s += 1
            if f > 0:
                print("percent success rate: ", s/(s+f))
            results = [x for x in results if x is not None]
            df_results = pd.DataFrame(results)
            df_results['column1']     = dfc.columns[i]
            df_results['column2']     = dfc.columns[j]
            df_results['chunk_index'] = chunk_index
            RESULTS.append(df_results)
    RESULTS = pd.concat(RESULTS)
    return RESULTS
# Serialize the function
serialized_func = dill.dumps(process_pair)
# Prepare a list of all pairs of column indices
pairs = [(i, j) for i in range(dfc.shape[1]) for j in range(dfc.shape[1]) if i != j]
# Create a pool of worker processes
from pathos.multiprocessing import ProcessingPool as Pool
with Pool(cpu_count()) as pool:
    results = pool.map(process_pair, pairs)
# Convert results to a dataframe
DF_results = pd.concat(results)
pre = name.replace('ccatime', f'_granger_causality')
DF_results.to_csv(os.path.join(folder, f'{pre}.csv'),
                  index=False)

# Let's add some summary statistics
pre = name.replace('ccatime', f'_granger_causality')
DF_results = pd.read_csv(os.path.join(folder, f'{pre}.csv'))
df = DF_results
df['uid'] = (df['chunk_index'].astype(str) + '_' + df['lag'].astype(str) + '_'
             + df['itest'].astype(str) + '_' + df['column1'] + '_' +
             df['column2'])
row_dict = df.set_index('uid').T.to_dict()
def directionality_test_vectorized(df, row_dict):
    print("Running directionality test")
    uid_reverse = (df['chunk_index'].astype(str) + '_' + df['lag'].astype(str)
                   + '_' + df['itest'].astype(str) + '_' + df['column2'] + '_'
                   + df['column1'])
    print("finding reverse")
    import pdb; pdb.set_trace()
    df_reverse = pd.DataFrame([row_dict[uid] for uid in uid_reverse])
    print("found reverse...")
    f_stat_xy, p_value_xy, r2_xy, r2_y  = (df['F'], df['pvalue'], df['r2_xy'],
                                           df['r2_y'])
    f_stat_yx, p_value_yx, r2_yx, r2_x = (df_reverse['F'],
                                          df_reverse['pvalue'],
                                          df_reverse['r2_xy'],
                                          df_reverse['r2_y'])
    print("calculating directionality")
    direction = np.where(f_stat_xy > f_stat_yx, "x->y", "y->x")
    print("calculating significance")
    significance = np.where(np.minimum(p_value_xy, p_value_yx) < 0.05,
                            "significant", "not significant")
    print("calculating magnitude")
    magnitude = np.abs(r2_xy - r2_yx)
    return pd.DataFrame({'direction': direction, 'significance': significance, 'magnitude': magnitude})
df_results = directionality_test_vectorized(DF_results, row_dict)
assert DF_results.shape[0] == df.shape[0]
DF_results = pd.concat([DF_results, df_results], axis=1)
print(f"saving {pre}_directionality.csv")
DF_results.to_csv(os.path.join(folder, f'{pre}_directionality.csv'), index=False)

# Plots

#  . .     . . .o         |          o          
# -+-+-    | | |.,---.,---|,---.. . ..,---.,---.
# -+-+-    | | |||   ||   ||   || | |||   ||   |
#  ` `     `-'-'``   '`---'`---'`-'-'``   '`---|
#                                          `---'

def capture_windows(df, trigger_column, quantile, w):

    # Remove df columns that are not numeric
    df = df.select_dtypes(include=[np.number])

    # Normalize all columns to be between 0 and 1
    df_normalized = (df - df.min()) / (df.max() - df.min())

    # Compute the threshold
    threshold = df[trigger_column].quantile(quantile)

    # Identify the rows where the trigger column crosses the threshold
    crossing_indices = df.index[df[trigger_column] > threshold]

    # Capture the windows
    windows = []
    for i in crossing_indices:
        if i-w >= 0 and i+w < len(df):
            # Check for gaps greater than 0.2 seconds in the time diff
            window = df_normalized.iloc[i-w:i+w]
            time_diff = window['time'].diff().abs()
            if not any(time_diff > 0.2):
                windows.append(window)

    return windows

wins = capture_windows(df, 'U1', 0.9, 15)


def plot_windows_heatmap(windows):
    # Compute the mean of each column in each window
    window_means = pd.DataFrame([window.mean() for window in windows])

    # Plot the mean of each column over time
    sns.heatmap(window_means.T, vmin=0, center=0.5 vmax=1, cmap='coolwarm')

def plot_windows(windows, quantile, w):
    fig, axs = plt.subplots(len(windows.columns), 1, figsize=(10, 20), constrained_layout=True)
    fig.suptitle(f'Windows for quantile {quantile} and window size {w}', fontsize=16)
    
    for i, col in enumerate(windows.columns):
        data = windows[col]
        mean = data.mean(axis=0)
        # compute 95% confidence interval (1.96 is approximately the 97.5 percentile of a standard normal distribution)
        ci = 1.96 * data.std(axis=0) / np.sqrt(data.shape[0])  
        axs[i].plot(mean, label=f'Mean {col}')
        axs[i].fill_between(range(len(mean)), (mean-ci), (mean+ci), color='b', alpha=.1)
        axs[i].legend(loc='upper right')
    
    plt.show()

plot_windows(captured_windows, quantile, w)

