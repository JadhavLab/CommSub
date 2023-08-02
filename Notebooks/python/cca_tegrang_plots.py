import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

animal   = ''
spectype = 'power'
folder   = '/Volumes/MATLAB-Drive/Shared/figures/tables/'

# Load the data
data = pd.read_csv(os.path.join(folder,
    f'{animal}{spectype}_granger_causality_directionality.csv')
)
# Check the first few rows of the dataframe
data.head()

def average_lower_diagonal(df):
    # Create a copy of the dataframe to avoid modifying the original
    df_copy = df.copy()
    # Get the values below the diagonal (including the diagonal)
    lower_diagonal = np.tril(df_copy.values)
    # Get the values above the diagonal (excluding the diagonal), and take their negatives
    upper_diagonal = -np.triu(df_copy.values, k=1)
    # Calculate the average of the two matrices, element-wise
    averaged_values = (lower_diagonal + upper_diagonal) / 2
    # Create a new dataframe with the averaged values
    df_averaged = pd.DataFrame(averaged_values, index=df_copy.index, columns=df_copy.columns)
    return df_averaged

# We need to apply the get_significant_interactions function to the data
# But first, we need to fix an issue with the function where it references a variable 'direction' that is not defined
# The line should be 'raise ValueError(f"Invalid direction: {group.direction.iloc[0]}")'
# Let's correct this and define a new version of the function
#
def get_significant_interactions(df, max_lag_order):
    # Filter the data to only include rows with the maximum lag order
    df_max_lag = df[df['lag'] == max_lag_order]
    # Group by 'column1', 'column2' and 'direction'
    grouped_max_lag = df_max_lag.groupby(['column1', 'column2', 'direction'])
    # Combine p-values and calculate average magnitude
    def combine_and_average(group):
        combined_pstat, combined_pvalue = scipy.stats.combine_pvalues(group['pvalue'], method='fisher')
        avg_magnitude = group['magnitude'].mean()
        if group.direction.iloc[0] == "x->y":
            pass
        elif group.direction.iloc[0] == "y->x":
            avg_magnitude = -avg_magnitude
            combined_pstat = -combined_pstat
        else:
            raise ValueError(f"Invalid direction: {group.direction.iloc[0]}")
        return pd.Series([combined_pstat, combined_pvalue, avg_magnitude, combined_pvalue < 0.05], index=['combined_pstat', 'combined_pvalue', 'avg_magnitude', 'is_significant'])
    # Apply the function to each group
    combined_max_lag = grouped_max_lag.apply(combine_and_average)
    # Filter out the non-significant combinations
    significant_max_lag = combined_max_lag[combined_max_lag['is_significant']]
    # Create a pivot table for the heatmap
    heatmap_data = significant_max_lag.pivot_table(index='column1', columns='column2', values='avg_magnitude')
    # Calculate signed logarithm of average magnitude
    significant_max_lag['signed_log_magnitude'] = np.sign(significant_max_lag['avg_magnitude']) * np.log10(np.abs(significant_max_lag['avg_magnitude']))
    # Create a pivot table for the heatmap
    heatmap_data_signed_log = significant_max_lag.pivot_table(index='column1', columns='column2', values='signed_log_magnitude')
    return significant_max_lag, heatmap_data, heatmap_data_signed_log



# Apply the function to the data for each lag
lag_max = data['lag'].max()
significant_max_lag, heatmap_data, heatmap_data_signed_log = {}, {}, {}
for lag in range(1, lag_max+1):
    significant_max_lag[lag], heatmap_data[lag], heatmap_data_signed_log[lag] = get_significant_interactions(data, max_lag_order=lag)

# Apply the average_lower_diagonal function to the heatmap data for each lag
averaged_heatmap_data = {lag: average_lower_diagonal(df) for lag, df in heatmap_data.items()}
averaged_heatmap_data_signed_log = {lag: average_lower_diagonal(df) for lag, df in heatmap_data_signed_log.items()}

# Let's visualize the first lag for each
sns.heatmap(averaged_heatmap_data[1], cmap="coolwarm", center=0)
plt.show()

sns.heatmap(averaged_heatmap_data_signed_log[1], cmap="coolwarm", center=0)
plt.show()


# Let's visualize the first lag for each
sns.heatmap(averaged_heatmap_data[1], cmap="coolwarm", center=0)
plt.show()

sns.heatmap(averaged_heatmap_data_signed_log[1], cmap="coolwarm", center=0)
plt.show()

# Re-create signed log versions from averaged raw dataframes
averaged_heatmap_data_signed_log = {lag: np.sign(df) * np.log10(np.abs(df)) for lag, df in averaged_heatmap_data.items()}

# Let's visualize the first lag for each again
sns.heatmap(averaged_heatmap_data[1], cmap="coolwarm", center=0)
plt.show()

sns.heatmap(averaged_heatmap_data_signed_log[1], cmap="coolwarm", center=0)
plt.show()


# Determine global min and max
global_min = min(df.min().min() for df in averaged_heatmap_data_signed_log.values())
global_max = max(df.max().max() for df in averaged_heatmap_data_signed_log.values())

def create_heatmap(df, lag, vmin, vmax, cmap="coolwarm"):
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(df, cmap=cmap, center=0, vmin=vmin, vmax=vmax)
    plt.title(f"Lag order: {lag}")
    plt.savefig(f"heatmap_lag_{lag}.png")
    plt.close()

# Generate and save heatmaps
for lag, df in averaged_heatmap_data_signed_log.items():
    create_heatmap(df, lag, global_min, global_max)

# Create GIF
images = []
for lag in range(1, lag_max + 1):
    images.append(imageio.imread(f"heatmap_lag_{lag}.png"))
# Append the first image to make the GIF loop
images.append(images[0])
imageio.mimsave('heatmap.gif', images, duration=1)  # duration is the time between frames in seconds

# --- Transfer entropy heatmap ---

def get_transfer_entropy_heatmap(df, max_lag_order):
    # Filter the data to only include rows with the maximum lag order
    df_max_lag = df[df['lag'] == max_lag_order]
    # Group by 'column1', 'column2' and 'direction'
    grouped_max_lag = df_max_lag.groupby(['column1', 'column2', 'direction'])
    # Calculate average transfer entropy
    def average_transfer_entropy(group):
        avg_transfer_entropy = group['transfer_entropy'].mean()
        if group.direction.iloc[0] == "y->x":
            avg_transfer_entropy = -avg_transfer_entropy
        return avg_transfer_entropy
    # Apply the function to each group
    averaged_max_lag = grouped_max_lag.apply(average_transfer_entropy)
    # Create a pivot table for the heatmap
    heatmap_data = averaged_max_lag.reset_index().pivot_table(index='column1', columns='column2', values=0)
    return heatmap_data

# Apply the function to the data for the first lag
transfer_entropy_heatmap_data = get_transfer_entropy_heatmap(data, max_lag_order=1)

# Visualize the heatmap
sns.heatmap(transfer_entropy_heatmap_data, cmap="coolwarm", center=0)
plt.show()

# --- Transfer entropy heatmap with signed log ---

# Create a signed log version of the transfer entropy heatmap data
transfer_entropy_heatmap_data_signed_log = np.sign(transfer_entropy_heatmap_data) * np.log10(np.abs(transfer_entropy_heatmap_data))

# Visualize the heatmap
sns.heatmap(transfer_entropy_heatmap_data_signed_log, cmap="coolwarm", center=0)
plt.show()


# --- Transfer entropy heatmap with signed log for each lag ---

def create_gif(df_dict, vmin=None, vmax=None, cmap="coolwarm", filename="heatmap.gif"):
    # Determine global min and max if not provided
    if vmin is None:
        vmin = min(df.min().min() for df in df_dict.values())
    if vmax is None:
        vmax = max(df.max().max() for df in df_dict.values())
        
    # Create and save each heatmap as an image
    for lag, df in df_dict.items():
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(df, cmap=cmap, center=0, vmin=vmin, vmax=vmax)
        plt.title(f"Lag order: {lag}")
        plt.savefig(f"heatmap_lag_{lag}.png")
        plt.close()
    
    # Create GIF from images
    images = [imageio.imread(f"heatmap_lag_{lag}.png") for lag in df_dict.keys()]
    # Append the first image to make the GIF loop
    images.append(images[0])
    imageio.mimsave(filename, images, duration=1)  # duration is the time between frames in seconds

# Create signed log versions of the heatmaps
transfer_entropy_heatmap_data_signed_log = {lag: np.sign(df) * np.log10(np.abs(df)) 
                                            for lag, df in transfer_entropy_heatmap_data.items()}

# Create GIF
create_gif(transfer_entropy_heatmap_data_signed_log, filename='transfer_entropy_heatmap.gif')

# ------------ IMPORTANT: TE FROM UV to Spectral components  -----------------
# Filter the data to include only the rows with 'column1' as one of the U and V
# components
uv_data = data[data['column1'].str.startswith(('U', 'V'))]

# Group by 'column1' and 'column2' and calculate the mean transfer entropy
grouped_uv_data = uv_data.groupby(['column1', 'column2'])['transfer_entropy'].mean()

# Reset the index of the DataFrame
grouped_uv_data = grouped_uv_data.reset_index()

# Create a new column 'spectral_type' based on the 'column2' column
grouped_uv_data['spectral_type'] = grouped_uv_data['column2'].str.extract('(Cavg|S1|S2|wpli)')

# Create a color map for the different spectral types
color_map = {'Cavg': 'blue', 'S1': 'green', 'S2': 'red', 'wpli': 'purple'}

# Create a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(data=grouped_uv_data, x='column1', y='transfer_entropy', hue='spectral_type', palette=color_map)
plt.title('Mean Transfer Entropy Between U and V Components and Other Columns')
plt.xlabel('Column1 (U and V Components)')
plt.ylabel('Mean Transfer Entropy')
plt.legend(title='Spectral Type')
plt.show()


# -------------- SPLIT BY LAG AND BY THETA DELTA RIPPLE ---------------

# Define a function to plot the bar plot for a specific lag order
def plot_transfer_entropy_barplot_dodge(lag_order):
    # Filter the data to include only the rows with 'column1' as one of the U and V components and the specified lag order
    uv_data_lag = data[(data['column1'].str.startswith(('U', 'V'))) & (data['lag'] == lag_order)]
    # Group by 'column1', 'column2', and 'lag' and calculate the mean transfer entropy
    grouped_uv_data_lag = uv_data_lag.groupby(['column1', 'column2', 'lag'])['transfer_entropy'].mean()
    # Reset the index of the DataFrame
    grouped_uv_data_lag = grouped_uv_data_lag.reset_index()
    # Create new columns 'spectral_type' and 'frequency_band' based on the 'column2' column
    grouped_uv_data_lag['spectral_type'] = grouped_uv_data_lag['column2'].str.extract('(Cavg|S1|S2|wpli)')
    grouped_uv_data_lag['frequency_band'] = grouped_uv_data_lag['column2'].str.extract('(theta|delta|ripple)')
    # Create a bar plot with subplot rows for lags and subplot columns for frequency bands, with bars dodged by spectral type
    g = sns.FacetGrid(grouped_uv_data_lag, row='lag', col='frequency_band', sharex=False)
    g.map_dataframe(sns.barplot, x='column1', y='transfer_entropy', hue='spectral_type', palette=color_map, dodge=True)
    plt.subplots_adjust(top=0.92)
    g.fig.suptitle('Mean Transfer Entropy Between U and V Components and Other Columns (Split by Frequency Band)', fontsize=16)
    g.set_axis_labels('Column1 (U and V Components)', 'Mean Transfer Entropy')
    g.add_legend(title='Spectral Type')
    plt.show()

# Plot the bar plot for the first 5 lag orders
for lag in range(1, 6):
    plot_transfer_entropy_barplot_dodge(lag)

# -------------- GRANGER VERION OF SPLIT ------------------------------
# Define a function to plot the bar plot for a specific lag order
def plot_F_barplot_dodge(lag_order):
    # Filter the data to include only the rows with 'column1' as one of the U and V components and the specified lag order
    uv_data_lag = data[(data['column1'].str.startswith(('U', 'V'))) & (data['lag'] == lag_order)]
    # Group by 'column1', 'column2', and 'lag' and calculate the mean F value
    grouped_uv_data_lag = uv_data_lag.groupby(['column1', 'column2', 'lag'])['F'].mean()
    # Reset the index of the DataFrame
    grouped_uv_data_lag = grouped_uv_data_lag.reset_index()
    # Create new columns 'spectral_type' and 'frequency_band' based on the 'column2' column
    grouped_uv_data_lag['spectral_type'] = grouped_uv_data_lag['column2'].str.extract('(Cavg|S1|S2|wpli)')
    grouped_uv_data_lag['frequency_band'] = grouped_uv_data_lag['column2'].str.extract('(theta|delta|ripple)')
    # Create a bar plot with subplot rows for lags and subplot columns for frequency bands, with bars dodged by spectral type
    g = sns.FacetGrid(grouped_uv_data_lag, row='lag', col='frequency_band', sharex=False)
    g.map_dataframe(sns.barplot, x='column1', y='F', hue='spectral_type', palette=color_map, dodge=True)
    plt.subplots_adjust(top=0.9)
    g.fig.suptitle('Mean F Value Between U and V Components and Other Columns (Split by Frequency Band)')
    g.set_axis_labels('Column1 (U and V Components)', 'Mean F Value')
    g.add_legend(title='Spectral Type')
    plt.show()

# Plot the bar plot for the first 5 lag orders
for lag in range(1, 6):
    plot_F_barplot_dodge(lag)

