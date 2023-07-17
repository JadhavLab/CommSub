import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.utils import resample
from tqdm import tqdm

plt.ion()
folder='/Volumes/MATLAB-Drive/Shared/figures/python/dimPred/'
if not os.path.exists(folder):
    os.makedirs(folder)

# Load the dataset
df = pd.read_csv('/Volumes/MATLAB-Drive/Shared/figures/tables/predDim.csv')
u = {key:np.unique(df[key]) for key in df.keys()}
df.loc[:,'highlow'] = df.name.apply(lambda x: 'high' if 'control' not in x else 'low')
# Display the first few rows of the dataset
df.head()

# ------------------------------
# Plot the samples of MEA vs. dimension for each dataset

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

