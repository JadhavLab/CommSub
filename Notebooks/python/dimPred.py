import pandas as pd

# Load the dataset
df = pd.read_csv('/Volumes/MATLAB-Drive/Shared/figures/tables/predDim.csv')

# Display the first few rows of the dataset
df.head()

# ------------------------------

import matplotlib.pyplot as plt
import seaborn as sns
plt.ion()

# List of unique datasets
dataset_ids = df['iDataset'].unique()

# Create a subplot for each dataset
fig, axs = plt.subplots(len(dataset_ids), figsize=(10, 5 * len(dataset_ids)))

for i, dataset_id in enumerate(dataset_ids):
    # Filter the data for the current dataset
    data_subset = df[df['iDataset'] == dataset_id]
    
    # Unique partitions in the current dataset
    partitions = data_subset['iP'].unique()
    
    # Plot each partition as a separate line
    for partition in partitions:
        partition_data = data_subset[data_subset['iP'] == partition]
        axs[i].plot(partition_data['dims'], partition_data['mea'], linestyle='dotted', c='black', alpha=0.5, markersize=
    
    # Set the title and labels
    axs[i].set_title(f'Animal: {data_subset["animal"].iloc[0]}, genH: {data_subset["genH"].iloc[0]}')
    axs[i].set_xlabel('Dimensions')
    axs[i].set_ylabel('MEA')

plt.tight_layout()
plt.show()


# ------------------------------

                    # Define colors for different direction types
direction_colors = {'hpc-hpc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'},
                    'hpc-pfc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'}}

# Lighten colors for 'hpc-pfc' direction
for genH in direction_colors['hpc-pfc']:
    direction_colors['hpc-pfc'][genH] = sns.light_palette(direction_colors['hpc-pfc'][genH], n_colors=3)[1]

# Create a subplot for each dataset
fig, axs = plt.subplots(len(dataset_ids), figsize=(10, 5 * len(dataset_ids)))

# Adjust the margins around the plot grid
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

for i, dataset_id in enumerate(dataset_ids):
    # Filter the data for the current dataset
    data_subset = df[df['iDataset'] == dataset_id]
    
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
        color = direction_colors[direction][direction_data['genH'].iloc[0]]
        
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
    axs[i].set_title(f'Animal: {direction_data["animal"].iloc[0]}, genH: {direction_data["genH"].iloc[0]}')
    axs[i].set_xlabel('Fractional Dimension')
    axs[i].set_ylabel('Bootstrap MEA')

plt.tight_layout()
plt.show()

# ------------------------------

# Create a subplot for each genH type
fig, axs = plt.subplots(len(unique_genH), figsize=(10, 5 * len(unique_genH)))

# Adjust the margins around the plot grid
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

for i, genH_type in enumerate(unique_genH):
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

# ------------------------------
from sklearn.utils import resample

# Find all unique genH types in the dataset
unique_genH = df['genH'].unique()

# Define colors for different direction types
direction_colors = {'hpc-hpc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'},
                    'hpc-pfc': {'coherence': '#0000FF', 'power': '#800080', 'wpli': '#FFA500'}}

# Lighten colors for 'hpc-pfc' direction
for genH in direction_colors['hpc-pfc']:
    direction_colors['hpc-pfc'][genH] = sns.light_palette(direction_colors['hpc-pfc'][genH], n_colors=3)[1]

# Adjust the number of bootstrap iterations
n_iterations = 100

# Create a subplot for each genH type
fig, axs = plt.subplots(len(unique_genH), figsize=(10, 5 * len(unique_genH)))

# Adjust the margins around the plot grid
plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9)

for i, genH_type in enumerate(unique_genH):
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

# ------------------------------


