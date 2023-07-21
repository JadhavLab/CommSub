import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import resample
from tqdm import tqdm

plt.ion()
savefolder='/Volumes/MATLAB-Drive/Shared/figures/python/dimPred/'
if not os.path.exists(savefolder):
    os.makedirs(savefolder)

# Load the dataset
df = pd.read_csv('/Volumes/MATLAB-Drive/Shared/figures/tables/predDim.csv')
u = {key:np.unique(df[key]) for key in df.keys()}
df.loc[:,'highlow'] = df.name.apply(lambda x: 'high' if 'control' not in x else 'low')
# Display the first few rows of the dataset
df.head()
df.loc[:,'patternname'] = df.name.apply(lambda x: x.replace('-control', ''))
df.query('genH != "wpli"', inplace=True)

# -------------------------------------------------------
# Plot the samples of MEA vs. dimension for each dataset
# -------------------------------------------------------

# List of unique datasets
dataset_ids = df['iDataset'].unique()

# Create a subplot for each dataset
fig, axs = plt.subplots(int(len(dataset_ids)/3), 3, figsize=(10, 5 * len(dataset_ids)), sharex=True, sharey=True)

for dataset_id in dataset_ids:
    
    # Filter the data for the current dataset
    data_subset = df[df['iDataset'] == dataset_id]

    j = np.where(data_subset.genH.iloc[0] == np.array(['power', 'coherence', 'wpli']))[0]
    i = np.where(data_subset.animal.iloc[0] == np.array(u["animal"]))[0]

    ax = axs[i,j][0]
    
    # Unique partitions in the current dataset
    partitions = data_subset['iP'].unique()
    
    # Plot each partition as a separate line
    for partition in partitions:
        partition_data = data_subset[data_subset['iP'] == partition]
        ax.plot(partition_data['dims'], partition_data['mea'], linestyle='dotted', c='black', alpha=0.2, markersize=1, lw=0.5)
        ax.set(xlim=(0,10))
    
    # Set the title and labels
    ax.set_title(f'Animal: {data_subset["animal"].iloc[0]}, genH: {data_subset["genH"].iloc[0]}')
    ax.set_xlabel('Dimensions')
    ax.set_ylabel('MEA')

plt.tight_layout()
plt.show()


# ------------------------------o
# Define colors for highlow types
def plot_mea_by_dimensions(df, direction_filter=None):
    """
    Function to plot Mean Error of Approximation (MEA) vs dimensions for different animal and genH types.
    The data can be filtered by a specific directionality.
    
    Parameters:
    - df: DataFrame containing the data
    - direction_filter: Optional; if specified, only data with this directionality will be plotted
    """
    
    # Define colors for highlow types
    highlow_colors = {'high': 'red', 'low': 'blue'}
    
    # If a direction filter is specified, filter the data accordingly
    if direction_filter is not None:
        df = df[df['direction'] == direction_filter]

    # List of unique animals
    unique_animals = df['animal'].unique()
    
    # List of unique genH types
    unique_genH = df['genH'].unique()
    
    # Create a subplot for each combination of animal and genH type
    fig, axs = plt.subplots(len(unique_animals), len(unique_genH), figsize=(10, 18), sharex=True, sharey=True)
    
    for i, animal in enumerate(unique_animals):
        for j, genH_type in enumerate(unique_genH):
            # Filter the data for the current animal and genH type
            data_subset = df[(df['animal'] == animal) & (df['genH'] == genH_type)]
            
            if not data_subset.empty:
                # Unique partitions in the current dataset
                partitions = data_subset['iP'].unique()
                
                # Plot each partition as a separate line
                for partition in partitions:
                    partition_data = data_subset[data_subset['iP'] == partition]
                    
                    # Define color based on highlow type and plot each highlow type as a separate line
                    for highlow in ["high", "low"]:
                        color = highlow_colors[highlow]
                        pD = partition_data[partition_data['highlow'] == highlow]
                        axs[i, j].plot(pD['dims'], pD['mea'], linestyle='dotted', 
                                       color=color, alpha=0.2, markersize=1, lw=0.5)
                        
                axs[i, j].set_xlim(0, 10)  # Set x-axis limits
                
                # Set the title and labels
                axs[i, j].set_title(f'Animal: {animal}, genH: {genH_type}')
                axs[i, j].set_xlabel('Dimensions')
                axs[i, j].set_ylabel('MEA')
    
    fig.suptitle(f'MEA vs. dimensions for different animal and genH types ({direction_filter})')
    plt.tight_layout()
    plt.show()

plot_mea_by_dimensions(df, direction_filter='hpc-hpc')
plot_mea_by_dimensions(df, direction_filter='hpc-pfc')


# ------------------------------
# Plot the samples of MEA vs. fractional dimension for each dataset
# colored by the direction type and the network pattern

def plot_bootstrap_mea_by_dataset(df, pattern_name=None, highlow_filter=None):
    """
    Function to plot bootstrap Mean Error of Approximation (MEA) vs fractional dimensions for each dataset.
    The data can be filtered by a specific highlow value and optionally split by pattern name.
    
    Parameters:
    - df: DataFrame containing the data
    - pattern_name: Optional; if specified, separate subplots will be created for each unique pattern name
    - highlow_filter: Optional; if specified, only data with this highlow value will be plotted
    """
    
    # If a highlow filter is specified, filter the data accordingly
    if highlow_filter is not None:
        df = df[df['highlow'] == highlow_filter]

    # If a pattern name is specified, filter the data accordingly
    if pattern_name is not None:
        df = df[df['name'] == pattern_name]

    # Define colors for different direction types
    direction_colors = {'hpc-hpc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'},
                        'hpc-pfc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'}}

    # Lighten colors for 'hpc-pfc' direction
    for genH in direction_colors['hpc-pfc']:
        direction_colors['hpc-pfc'][genH] = sns.light_palette(direction_colors['hpc-pfc'][genH], n_colors=3)[1]

    # List of unique datasets
    dataset_ids = df['iDataset'].unique()

    # Create a subplot for each dataset and each unique pattern name
    fig, axs = plt.subplots(len(dataset_ids), len(df['name'].unique()), figsize=(10 * len(df['name'].unique()), 5 * len(dataset_ids)))

    # Adjust the margins around the plot grid
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

    # Set the background color to dark gray
    fig.patch.set_facecolor('xkcd:dark grey')

    for i, dataset_id in enumerate(dataset_ids):
        for j, pattern in enumerate(df['name'].unique()):
            # Filter the data for the current dataset and pattern name
            data_subset = df[(df['iDataset'] == dataset_id) & (df['name'] == pattern)]

            # Unique directions in the current dataset
            directions = data_subset['direction'].unique()

            # Plot each direction as a separate line
            for direction in directions:
                direction_data = data_subset[data_subset['direction'] == direction]

                # Initialize arrays to store bootstrapped means and standard deviations
                bootstrap_means = np.zeros((n_iterations, len(direction_data['fracDim'].unique())))
                bootstrap_stds = np.zeros_like(bootstrap_means)

                # Perform bootstrapping
                for k in range(n_iterations):
                    bootstrap_sample = resample(direction_data, replace=True)

                    # Compute the mean and standard deviation for each fractional dimension
                    for l, (fracDim, group_data) in enumerate(bootstrap_sample.groupby('fracDim')):
                        bootstrap_means[k, l] = group_data['mea'].mean()
                        bootstrap_stds[k, l] = group_data['mea'].std()

                # Compute the mean and standard deviation of the bootstrapped means and standard deviations
                mean_bootstrap_means = bootstrap_means.mean(axis=0)
                mean_bootstrap_stds = bootstrap_stds.mean(axis=0)

                # Define color based on genH type and direction
                color = direction_colors[direction][direction_data['genH'].iloc[0]]

                # Plot the mean and standard deviation for each fractional dimension as a shaded region
                # Set the text color of the axes to white
                axs[i, j].xaxis.label.set_color('white')
                axs[i, j].yaxis.label.set_color('white')
                axs[i, j].tick_params(axis='x', colors='white')
                axs[i, j].tick_params(axis='y', colors='white')
                axs[i, j].set_facecolor('xkcd:dark grey')
                axs[i, j].plot(direction_data['fracDim'].unique(), mean_bootstrap_means, color=color)
                axs[i, j].fill_between(direction_data['fracDim'].unique(), 
                                       mean_bootstrap_means - mean_bootstrap_stds, 
                                       mean_bootstrap_means + mean_bootstrap_stds, 
                                       color=color, alpha=0.2)

                # Plot a vertical line at the fractional dimension corresponding to the optimal dimension
                opt_fracDim = direction_data['fracDim'].iloc[direction_data['dims'].tolist().index(direction_data['optDim'].iloc[0])]
                axs[i, j].axvline(opt_fracDim, color=color, linestyle='--', linewidth=2.0)  # increased linewidth for boldness

            # Set the title and labels
            axs[i, j].set_title(f'Animal: {direction_data["animal"].iloc[0]}, genH: {direction_data["genH"].iloc[0]}, pattern: {pattern}')
            axs[i, j].set_xlabel('Fractional Dimension')
            axs[i, j].set_ylabel('Bootstrap MEA')

    plt.tight_layout()
    plt.show()


plot_bootstrap_mea_by_dataset(df, highlow_filter='high')
plot_bootstrap_mea_by_dataset(df, highlow_filter='low')


# ------------------------------

# Here's the corrected function

def plot_bootstrap_mea_by_dataset(df, name_split=False, highlow_filter=None, n_iterations=100,
                                  frac_dim_filter=None):
    """
    Function to plot bootstrap Mean Error of Approximation (MEA) ratio of 'hpc-pfc' over 'hpc-hpc' vs fractional dimensions for each dataset.
    The data can be filtered by a specific highlow value and optionally split by name.

    Parameters:
    - df: DataFrame containing the data
    - name_split: Optional; if True, separate subplots will be created for each unique name
    - highlow_filter: Optional; if specified, only data with this highlow value will be plotted
    - n_iterations: Number of bootstrap iterations
    """

    # If a highlow filter is specified, filter the data accordingly
    if highlow_filter is not None:
        df = df[df['highlow'] == highlow_filter]
    if frac_dim_filter:
        df = df[df['fracDim'] < frac_dim_filter]

    # Define colors for different names
    name_colors = {'theta': 'green', 'ripple': 'blue', 'delta': 'red'}
    for name in df['name'].unique():
        if 'control' in name:
            name_colors[name] = sns.light_palette(name_colors[name.split('-')[0]], n_colors=3)[1]
        else:
            name_colors[name] = name_colors[name]

    # List of unique datasets
    dataset_ids = df['iDataset'].unique()

    # Create a subplot for each dataset and each unique name if name_split is True
    if name_split:
        fig, axs = plt.subplots(len(dataset_ids), len(df['name'].unique()), figsize=(10 * len(df['name'].unique()), 5 * len(dataset_ids)), sharey=True)
    else:
        fig, axs = plt.subplots(len(dataset_ids), 1, figsize=(10, 5 * len(dataset_ids)), sharey=True)

    # Adjust the margins around the plot grid
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

    # Set the background color to dark gray
    fig.patch.set_facecolor('xkcd:dark grey')

    for i, dataset_id in enumerate(dataset_ids):
        if name_split:
            for j, name in enumerate(df['name'].unique()):
                # Filter the data for the current dataset and name
                data_subset = df[(df['iDataset'] == dataset_id) & (df['name'] == name)]
                plot_ratio(data_subset, name_colors, n_iterations, i, j, axs, name)
        else:
            # Filter the data for the current dataset
            data_subset = df[df['iDataset'] == dataset_id]
            plot_ratio(data_subset, name_colors, n_iterations, i, 0, axs)

    plt.tight_layout()
    plt.show()

def plot_ratio(data_subset, name_colors, n_iterations, i, j, axs, name=None):
    # Unique directions in the current dataset
    directions = data_subset['direction'].unique()

    if 'hpc-hpc' in directions and 'hpc-pfc' in directions:
        # Extract the data for the two directions
        direction_data_0 = data_subset[data_subset['direction'] == 'hpc-hpc']
        direction_data_1 = data_subset[data_subset['direction'] == 'hpc-pfc']

        # Initialize arrays to store bootstrapped means
        bootstrap_means_0 = np.zeros((n_iterations, len(direction_data_0['fracDim'].unique())))
        bootstrap_means_1 = np.zeros((n_iterations, len(direction_data_1['fracDim'].unique())))

        # Perform bootstrapping
        for k in range(n_iterations):
            bootstrap_sample_0 = resample(direction_data_0, replace=True)
            bootstrap_sample_1 = resample(direction_data_1, replace=True)

            # Compute the mean for each fractional dimension
            for l, (fracDim, group_data_0) in enumerate(bootstrap_sample_0.groupby('fracDim')):
                group_data_1 = bootstrap_sample_1[bootstrap_sample_1['fracDim'] == fracDim]
                if not group_data_1.empty:
                    bootstrap_means_0[k, l] = group_data_0['mea'].mean()
                    bootstrap_means_1[k, l] = group_data_1['mea'].mean()

        # Compute the mean of the bootstrapped means
        mean_bootstrap_means_0 = bootstrap_means_0.mean(axis=0)
        mean_bootstrap_means_1 = bootstrap_means_1.mean(axis=0)

        # Compute the MEA ratio
        ratio_mea = mean_bootstrap_means_1 / mean_bootstrap_means_0

        # Compute upper and lower bounds of the 95% confidence interval
        lower_bound = np.percentile(ratio_mea, 2.5)
        upper_bound = np.percentile(ratio_mea, 97.5)

        # Define color based on name
        color = name_colors[name]

        # Plot the ratio MEA for each fractional dimension
        axs[i, j].plot(direction_data_0['fracDim'].unique(), ratio_mea, color=color, linewidth=2)
        axs[i, j].fill_between(direction_data_0['fracDim'].unique(), lower_bound, upper_bound, color=color, alpha=0.2)
        axs[i, j].axhline(1, color='white', linestyle='--', linewidth=1)

        # Set the text color of the axes to white
        axs[i, j].xaxis.label.set_color('white')
        axs[i, j].yaxis.label.set_color('white')
        axs[i, j].tick_params(axis='x', colors='white')
        axs[i, j].tick_params(axis='y', colors='white')
        axs[i, j].title.set_color('white')
        axs[i, j].set_facecolor('xkcd:dark grey')

        # Set the title and labels
        if name is not None:
            axs[i, j].set_title(f'Dataset: {i+1}, Name: {name}', color='white')
        axs[i, j].set_xlabel('Fractional Dimension', color='white')
        axs[i, j].set_ylabel('Bootstrap MEA Ratio', color='white')

plot_bootstrap_mea_by_dataset(df, name_split=True, highlow_filter='high',
                              frac_dim_filter=0.5)




# ------------------------------

def plot_bootstrap_mea(df, highlow_filter=None):
    """
    Function to plot bootstrap Mean Error of Approximation (MEA) vs fractional dimensions for different genH types.
    The data can be filtered by a specific highlow value.
    
    Parameters:
    - df: DataFrame containing the data
    - highlow_filter: Optional; if specified, only data with this highlow value will be plotted
    """
    
    # If a highlow filter is specified, filter the data accordingly
    if highlow_filter is not None:
        df = df[df['highlow'] == highlow_filter]

    # List of unique genH types
    unique_genH = df['genH'].unique()

    # Create a subplot for each genH type
    fig, axs = plt.subplots(len(unique_genH), figsize=(10, 18))

    # Adjust the margins around the plot grid
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

    for i, genH_type in tqdm(enumerate(unique_genH), total=len(unique_genH)):
        # Filter the data for the current genH type
        data_subset = df[df['genH'] == genH_type]

        # Unique directions in the current dataset
        directions = data_subset['direction'].unique()

        # Plot each direction as a separate line
        for direction in directions:
            direction_data = data_subset[data_subset['direction'] == direction]

            # Initialize arrays to store bootstrapped means and standard deviations
            bootstrap_means = np.zeros((n_iterations, len(direction_data['fracDim'].unique())))
            bootstrap_stds = np.zeros_like(bootstrap_means)

            # Perform bootstrapping
            for j in range(n_iterations):
                bootstrap_sample = resample(direction_data, replace=True)

                # Compute the mean and standard deviation for each fractional dimension
                for k, (fracDim, group_data) in enumerate(bootstrap_sample.groupby('fracDim')):
                    bootstrap_means[j, k] = group_data['mea'].mean()
                    bootstrap_stds[j, k] = group_data['mea'].std()

            # Compute the mean and standard deviation of the bootstrapped means and standard deviations
            mean_bootstrap_means = bootstrap_means.mean(axis=0)
            mean_bootstrap_stds = bootstrap_stds.mean(axis=0)

            # Define color based on genH type and direction
            color = direction_colors[direction][genH_type]

            # Plot the mean and standard deviation for each fractional dimension as a shaded region
            axs[i].plot(direction_data['fracDim'].unique(), mean_bootstrap_means, color=color)
            axs[i].fill_between(direction_data['fracDim'].unique(), 
                                mean_bootstrap_means - mean_bootstrap_stds, 
                                mean_bootstrap_means + mean_bootstrap_stds, 
                                color=color, alpha=0.2)

            # Plot a vertical line at the fractional dimension corresponding to the optimal dimension
            opt_fracDim = direction_data['fracDim'].iloc[direction_data['dims'].tolist().index(direction_data['optDim'].iloc[0])]
            axs[i].axvline(opt_fracDim, color=color, linestyle='--', linewidth=2.0)  # increased linewidth for boldness

        # Set the title and labels
        axs[i].set_title(f'genH: {genH_type}')
        axs[i].set_xlabel('Fractional Dimension')
        axs[i].set_ylabel('Bootstrap MEA')

    plt.tight_layout()
    plt.show()

plot_bootstrap_mea(df, highlow_filter='high')


# ------------------------------
from sklearn.utils import resample
# Find all unique genH types in the dataset
unique_genH = df['genH'].unique()
# ------------------------------


# Function to calculate bootstrap means of 'mea' for a given highlow category and genH
def calculate_bootstrap_means(df, highlow_category, genH, num_samples=1000, sample_size=500):
    # Initialize list to hold the bootstrap means
    bootstrap_means = []

    # Resample and compute means
    for _ in range(num_samples):
        sample = df[(df['highlow'] == highlow_category) & (df['genH'] == genH)]['mea'].sample(sample_size, replace=True)
        bootstrap_means.append(sample.mean())
    
    return bootstrap_means

# Create a subplot for each unique genH
fig, axs = plt.subplots(len(df['genH'].unique()), figsize=(10, 5 * len(df['genH'].unique())))

# Adjust the margins around the plot grid
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

# Set the background color to dark gray
fig.patch.set_facecolor('xkcd:dark grey')

for i, genH in tqdm(enumerate(df['genH'].unique()), 
                    total=len(df['genH'].unique())):
    # Calculate bootstrap means for high and low categories
    high_means = calculate_bootstrap_means(df, 'high', genH)
    low_means = calculate_bootstrap_means(df, 'low', genH)

    # Plot histograms of bootstrap means
    axs[i].hist(high_means, bins=50, alpha=0.5, label='high')
    axs[i].hist(low_means, bins=50, alpha=0.5, label='low')

    # Set the title and labels
    axs[i].set_title(f'{genH}: Distribution of Bootstrap Sample Means')
    axs[i].set_xlabel('Sample Mean of MEA')
    axs[i].set_ylabel('Frequency')

    # Set the text color of the axes to white
    axs[i].xaxis.label.set_color('white')
    axs[i].yaxis.label.set_color('white')
    axs[i].tick_params(axis='x', colors='white')
    axs[i].tick_params(axis='y', colors='white')
    axs[i].set_facecolor('xkcd:dark grey')

    axs[i].legend()

plt.tight_layout()
plt.show()


def calculate_bootstrap_means(df, highlow_category, genH, num_samples=1000, sample_size=500):
    # Initialize list to hold the bootstrap means
    bootstrap_means = []

    # Resample and compute means
    for _ in range(num_samples):
        sample = df[(df['highlow'] == highlow_category) & (df['genH'] == genH)]['mea']
        if len(sample) > 0:
            sample = sample.sample(sample_size, replace=True)
            bootstrap_means.append(sample.mean())
    
    return bootstrap_means


def plot_bootstrap_means(df, split_by_name=False, num_samples=1000, sample_size=500):
    """
    Function to plot bootstrap means of 'mea' for high and low 'highlow' categories, optionally split by 'name'.
    The data can be split by 'name', resulting in separate subplots for each unique name.

    Parameters:
    - df: DataFrame containing the data
    - split_by_name: Optional; if True, separate subplots will be created for each unique name
    - num_samples: Number of bootstrap samples to draw; default is 1000
    - sample_size: Size of each bootstrap sample; default is 500
    """

    df = df.query('fracDim > fracOptDim')
    
    # List of unique genH and names
    genHs = df['genH'].unique()
    names = df['patternname'].unique() if split_by_name else [None]

    # Create a subplot for each unique genH and name
    fig, axs = plt.subplots(len(genHs), len(names), figsize=(10 * len(names), 5 * len(genHs)))

    # Adjust the margins around the plot grid
    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

    # Set the background color to dark gray
    fig.patch.set_facecolor('xkcd:dark grey')

    for i, genH in tqdm(enumerate(genHs), total=len(genHs), desc='genH'):
        for j, name in tqdm(enumerate(names), total=len(names), desc='pattern'):
            # Filter the data for the current genH and name (if specified)
            data_subset = df[df['genH'] == genH] if name is None else df[(df['genH'] == genH) & (df['patternname'] == name)]

            # Calculate bootstrap means for high and low categories
            high_means = calculate_bootstrap_means(data_subset, 'high', genH, num_samples, sample_size)
            low_means = calculate_bootstrap_means(data_subset, 'low', genH, num_samples, sample_size)

            # Plot histograms of bootstrap means
            if len(names) > 1:  # if split_by_name is True
                ax = axs[i, j]
            else:
                ax = axs[i]

            ax.hist(high_means, bins=50, alpha=0.5, label='high', color='xkcd:red')
            ax.hist(low_means, bins=50, alpha=0.5, label='low', color='xkcd:blue')

            # Set the title and labels
            title = f'{genH}: Distribution of Bootstrap Sample Means'
            if name is not None:
                title += f' for {name}'
            ax.set_title(title)
            ax.set_xlabel('Sample Mean of MEA')
            ax.set_ylabel('Frequency')

            # Set the text color of the axes to white
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.set_facecolor('xkcd:dark grey')

            ax.legend()

    plt.tight_layout()
    plt.show()
    return fig, axs

# Use the function with example parameters

fig, axs = plot_bootstrap_means(df.query('direction == "hpc-hpc"'), split_by_name=True)
minmax = [(ax.get_xlim()[0], ax.get_xlim()[1]) for ax in plt.gcf().axes]
minmax = (min([x[0] for x in minmax]), max([x[1] for x in minmax]))
[ax.set_xlim(minmax) for ax in plt.gcf().axes]
fig.suptitle('HPC-HPC\nBootstrap Means of MEA for High and Low Highlow Categories,\nSplit by Pattern/GenH', fontsize=16, fontweight='bold', color='white')
fig.subplots_adjust(top=0.9)
fig.savefig(os.path.join(savefolder, 'bootstrap_means_hpc-hpc.png'), dpi=300)

fig, axs = plot_bootstrap_means(df.query('direction == "hpc-pfc"'), split_by_name=True)
minmax = [(ax.get_xlim()[0], ax.get_xlim()[1]) for ax in plt.gcf().axes]
minmax = (min([x[0] for x in minmax]), max([x[1] for x in minmax]))
[ax.set_xlim(minmax) for ax in plt.gcf().axes]
fig.suptitle('HPC-PFC\nBootstrap Means of MEA for High and Low Highlow Categories,\nSplit by Pattern/GenH', fontsize=16, fontweight='bold', color='white')
fig.subplots_adjust(top=0.9)
fig.savefig(os.path.join(savefolder, 'bootstrap_means_hpc-pfc.png'), dpi=300)

def plot_3d_bootstrap_means(df, split_by_name=False, num_samples=1000, sample_size=500):
    """
    Function to plot 3D histograms of bootstrap means of 'mea' for high and low 'highlow' categories, optionally split by 'name'.
    The x-axis is the bootstrap mean, the y-axis is the frequency, and the z-axis is the 'fracDim' value.
    Each slice of the 3rd dimension is a histogram measured at a 'fracDim' value.

    Parameters:
    - df: DataFrame containing the data
    - split_by_name: Optional; if True, separate subplots will be created for each unique name
    - num_samples: Number of bootstrap samples to draw; default is 1000
    - sample_size: Size of each bootstrap sample; default is 500
    """
    df = df.query('fracDim > fracOptDim')

    # List of unique genH and names
    genHs = df['genH'].unique()
    names = df['patternname'].unique() if split_by_name else [None]
    fig = plt.figure(figsize=(10, 6))
    rows = len(genHs)
    cols = len(names)

    # Create a 3D plot for each unique genH and name
    for i, genH in tqdm(enumerate(genHs), total=len(genHs), desc='genH'):
        for j, name in tqdm(enumerate(names), total=len(names), desc='pattern'):
            # Filter the data for the current genH and name (if specified)
            data_subset = df[df['genH'] == genH] if name is None \
                    else df[(df['genH'] == genH) & (df['patternname'] == name)]

            # Create 3D histogram
            ax = fig.add_subplot(rows, cols, i * cols + j + 1, projection='3d')
            # Set the title and labels
            title = f'{genH}: 3D Distribution of Bootstrap Sample Means'
            if name is not None:
                title += f' for {name}'
            ax.set_title(title)
            ax.set_xlabel('Sample Mean of MEA')
            ax.set_ylabel('fracDim')
            ax.set_zlabel('Frequency')
            
            # List of unique 'fracDim' values
            fracDims = data_subset['fracDim'].unique()

            # Calculate bootstrap means and plot histogram for each 'fracDim' value
            for k, fracDim in enumerate(fracDims):
                # Filter the data for the current 'fracDim'
                fracDim_data = data_subset[data_subset['fracDim'] == fracDim]

                # Calculate bootstrap means for high and low categories
                high_means = calculate_bootstrap_means(fracDim_data, 'high', genH, num_samples, sample_size)
                low_means = calculate_bootstrap_means(fracDim_data, 'low', genH, num_samples, sample_size)
                if len(high_means) == 0 or len(low_means) == 0:
                    continue

                # Define bins
                bins = np.linspace(min(min(high_means), min(low_means)), max(max(high_means), max(low_means)), 50)
                print("bins", bins)

                # Plot histograms for high and low categories
                hist_high, _ = np.histogram(high_means, bins=bins)
                hist_low, _  = np.histogram(low_means, bins=bins)

                # Add bar collections
                # ax.bar(bins[:-1], hist_high, zs=fracDim, zdir='y', alpha=0.3, color='r', label='high' if k == 0 else "")
                # ax.bar(bins[:-1], hist_low,  zs=fracDim, zdir='y', alpha=0.3, color='b', label='low' if k == 0 else "")
                ax.bar3d(bins[:-1], np.full(len(bins) - 1, fracDim), np.zeros(len(bins) - 1), 0.5, 0.5, hist_high,color='r', alpha=0.3, label='high' if k == 0 else "")
                ax.bar3d(bins[:-1], np.full(len(bins) - 1, fracDim), np.zeros(len(bins) - 1), 0.5, 0.5, hist_low, color='b', alpha=0.3, label='low' if k == 0 else "")

    plt.show()


plot_3d_bootstrap_means(df.query('direction == "hpc-hpc"'), split_by_name=True) 
plot_3d_bootstrap_means(df.query('direction == "hpc-pfc"'), split_by_name=True)
