import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.utils import resample
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm
import itertools

# Set the flag for shading the confidence intervals to False
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
shade_confidence_intervals   = False
stratify_spectral_components = True  # Change this to False if you don't want to stratify and shade the spectral components
ci = True

# Define the trajbounds for each column of the subplot grid
column_trajbounds = [0, 1]

# Update the components for each row of the subplot grid
row_components = [["U1", "U2", "U3"], ["V1", "V2", "V3"], ["Cavgtheta",
                                                           "S1theta",
                                                           "S2theta",
                                                           "wpli_avgtheta"],
                  ["Cavgdelta", "S1delta", "S2delta", "wpli_avgdelta"],
                  ["Cavgripple", "S1ripple", "S2ripple", "wpli_avgripple"]]


component_colors = {
    "U1": "darkred", "U2": "red", "U3": "lightcoral",
    "V1": "darkblue", "V2": "blue", "V3": "lightblue",
    "Cavgtheta": "black", "Cavgdelta": "black", "Cavgripple": "black",
    "S1theta": "darkred", "S1delta": "darkred", "S1ripple": "darkred",
    "S2theta": "darkblue", "S2delta": "darkblue", "S2ripple": "darkblue",
    "wpli_avgtheta": "black", "wpli_avgdelta": "black", "wpli_avgripple": "black"
}

component_fill_bases = {
    "U1": 0, "U2": 1, "U3": 2,
    "V1": 0, "V2": 1, "V3": 2,
    "Cavgtheta": 0, "Cavgdelta": 0, "Cavgripple": 0,
    "S1theta": 1, "S1delta": 1, "S1ripple": 1,
    "S2theta": 1, "S2delta": 1, "S2ripple": 1,
    "wpli_avgtheta": 0, "wpli_avgdelta": 0, "wpli_avgripple": 0
}
# Scale distances
if ci:
    between_scale = 2
else:
    between_scale = 5/8
component_fill_bases = {k: v * between_scale for k, v in component_fill_bases.items()}

scale = 2

line_styles = {
    "U1": "-", "U2": "-", "U3": "-",
    "V1": "-", "V2": "-", "V3": "-",
    "Cavgtheta": "-", "Cavgdelta": "-", "Cavgripple": "-",
    "S1theta": "-", "S1delta": "-", "S1ripple": "-",
    "S2theta": "-", "S2delta": "-", "S2ripple": "-",
    "wpli_avgtheta": "--", "wpli_avgdelta": "--", "wpli_avgripple": "--"
}
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
# - - - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - - 
print("Loading data...")
folder = '/Volumes/MATLAB-Drive/Shared/figures/tables/'
# name   = 'ZT2powerccatime'
name   = 'powerccatime'
append = '_50bin'
datafile = os.path.join(folder, f'{name}.*')
datafile = glob.glob(datafile)[0]
if os.path.splitext(datafile)[1] == '.csv':
    df = pd.read_csv(datafile)
elif os.path.splitext(datafile)[1] == '.parquet':
    df = pd.read_parquet(datafile)
else:
    raise ValueError(f'Unknown file extension: {os.path.splitext(datafile)[1]}')
print("Done.")
def label_epochs(df, time_col='time', time_threshold=4):
    # Calculate time difference
    df.loc[:,'time_diff'] = df[time_col].diff()
    # Initialize the epoch column with 0
    df.loc[:,'epoch'] = 0
    # If time difference is more than threshold, increment epoch
    df.loc[df['time_diff'] > time_threshold, 'epoch'] = 1
    # Calculate the cumulative sum of the epoch column
    df.loc[:, 'epoch'] = df['epoch'].cumsum()
    # Drop the time_diff column as we don't need it anymore
    df.drop(columns=['time_diff'], inplace=True)
    return df
# Test the function on your dataframe
df = df.sort_values(['time','animal']).groupby('animal').apply(label_epochs).reset_index()
if not name.endswith("_epoch"):
    name = name + "_epoch"

# Number of bins and bootstrap samples
n_bins = 45
n_bootstrap_samples = 500

# Update the columns to bootstrap
columns_to_bootstrap = ["U1", "U2", "U3", "V1", "V2", "V3", "Cavgtheta", 
                        "Cavgdelta", "Cavgripple", "S1theta", "S1delta", 
                        "S2theta", "wpli_avgtheta", "S1ripple", "S2ripple",
                        "S2delta", "wpli_avgripple", "wpli_avgdelta"]

# Reinitialize the DataFrame to hold the results
bootstrap_means_combined = []
# Bin the lindist data
df["lindist_bin"] = pd.cut(df["lindist"], bins=n_bins)

# Rerun the bootstrap for each bin and each column, for each trajectory bound
print("Running bootstrap with {} bins and {} bootstrap samples...".format(n_bins, n_bootstrap_samples))
for trajbound in tqdm([0, 1], desc="trajbound", total=2):
    # Filter the data based on trajbound
    data_trajbound = df[df["trajbound"] == trajbound]
    
    # Loop over the unique combinations of column, trajbound, and lindist_bin
    iters = list(itertools.product(columns_to_bootstrap, 
                                   data_trajbound["lindist_bin"].unique().categories,
                                   data_trajbound["epoch"].unique()))

    for (column, bin_label, epoch) in tqdm(iters, desc="column, lindist_bin", total=len(iters)):
        # Get the data for this column, trajbound, and bin_label
        data = data_trajbound.loc[(data_trajbound["lindist_bin"] == bin_label) &
                                  (data_trajbound["epoch"] == epoch),
                                  ["animal",column]]
        
        # Find the minimum number of points available for each animal
        min_points_per_animal = data.groupby("animal").size().min()
        
        # Generate bootstrap samples
        for iboot in range(n_bootstrap_samples):
            # Initialize an empty list to hold the data for this bootstrap sample
            bootstrap_sample = []
            
            # Sample the minimum number of points from each animal's data
            for animal, animal_data in data.groupby("animal"):
                sample = (animal_data.sample(min_points_per_animal, replace=True))
                
                # Compute the mean of the bootstrap sample
                bootstrap_mean = sample.mean(numeric_only=True).astype(float)
                
                # Add the result to the DataFrame
                bootstrap_means_combined.append({
                    "iboot": iboot,
                    "lindist_bin": bin_label,
                    "column": column,
                    "bootstrap_mean": bootstrap_mean,
                    "trajbound": trajbound,
                    "animal": animal
                })


# Convert the list of dictionaries to a DataFrame
bootstrap_means_combined = pd.DataFrame(bootstrap_means_combined)
# Create a new lindist_bin_ind column
bootstrap_means_combined["lindist_bin_ind"] = bootstrap_means_combined["lindist_bin"].apply(lambda x: x.right)
bootstrap_means_combined.loc[:,'bootstrap_mean'] = bootstrap_means_combined.bootstrap_mean.astype(float)
bootstrap_means_combined.to_parquet(os.path.join(folder, f'{name}_bootstrap{append}.parquet'), index=False)

# Smooth
bootstrap_means_combined.sort_values(by=["animal", "column", "trajbound", "iboot", "lindist_bin_ind"], inplace=True)
bootstrap_means_combined["bootstrap_mean_smooth"] = bootstrap_means_combined.groupby(["animal","column", "trajbound", "iboot"])["bootstrap_mean"].transform(lambda x: x.rolling(7, 1).mean())
bootstrap_means_combined["bootstrap_mean_smooth"] = bootstrap_means_combined.groupby(["animal","column", "trajbound", "iboot"])["bootstrap_mean_smooth"].transform(lambda x: x.interpolate())

# ----------------------------------------------------
# Normalize the bootstrap_mean values
# ----------------------------------------------------
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
scaler = MinMaxScaler()
iters = itertools.product(columns_to_bootstrap, bootstrap_means_combined.animal.unique())
for column, animal in tqdm(iters, desc="feature engineering", total=len(columns_to_bootstrap) * bootstrap_means_combined.animal.nunique()):
    # Get the data for this column and animal
    data = \
            bootstrap_means_combined[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal)]["bootstrap_mean"].values.reshape(-1, 1)
    data_smooth  = \
            bootstrap_means_combined[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal)]["bootstrap_mean_smooth"].values.reshape(-1, 1)
    # Scale the data
    scaled_data        = scaler.fit_transform(data)
    scaled_data_smooth = scaler.fit_transform(data_smooth)
    # Update the DataFrame
    bootstrap_means_combined.loc[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal), "bootstrap_mean"] = scaled_data
    bootstrap_means_combined.loc[(bootstrap_means_combined["column"] == column) & \
            (bootstrap_means_combined["animal"] == animal), "bootstrap_mean_smooth"] = scaled_data_smooth
# bootstrap_means_combined.head()
bootstrap_means_combined.to_parquet(
        os.path.join(folder, f'{name}_bootstrap_normalized{append}.parquet'), 
                                    index=False)
# ----------------------------------------------------
# Read Parquet
# ----------------------------------------------------
bootstrap_means_combined = pd.read_parquet(os.path.join(folder, f'{name}_bootstrap_normalized{append}.parquet'))

# Add the lindist_bin_mid column back to the DataFrame
print("Adding lindist_bin_mid column...")
bootstrap_means_combined["lindist_bin_mid"] = \
    bootstrap_means_combined["lindist_bin"].apply(lambda x: x.mid)


# -------------------------Checking for duplicates------------------------- ##

component_fill_bases = {
    "U1": 0, "U2": 1, "U3": 2,
    "V1": 0, "V2": 1, "V3": 2,
    "Cavgtheta": 0, "Cavgdelta": 0, "Cavgripple": 0,
    "S1theta": 1, "S1delta": 1, "S1ripple": 1,
    "S2theta": 1, "S2delta": 1, "S2ripple": 1,
    "wpli_avgtheta": 0, "wpli_avgdelta": 0, "wpli_avgripple": 0
}
# Scale distances
scale = 1
skipripple = True
if ci:
    between_scale = 5/8 
else:
    between_scale = 5/8
component_fill_bases = {k: v * between_scale for k, v in component_fill_bases.items()}

animal_bootstrap_means_combined = bootstrap_means_combined.groupby(["iboot", "column", "trajbound", "lindist_bin_mid"]).mean().reset_index()
rows = 4 if skipripple else 5

for epoch in tqdm(animal_bootstrap_means_combined["epoch"].unique(), desc="epoch"):

    D = animal_bootstrap_means_combined[animal_bootstrap_means_combined["epoch"] == epoch]

    fig, axes = plt.subplots(nrows=rows, ncols=2, figsize=(15, 25))
    smooth = True
    field = "bootstrap_mean_smooth" if smooth else "bootstrap_mean"
    # Set the overall title
    fig.suptitle('Epoch=' + epoch, fontsize=20)
    for i, components in enumerate(row_components[:rows]):
        for j, trajbound in enumerate(column_trajbounds):
            # Get the data for this subplot
            data = D[D["column"].isin(components) & (D["trajbound"] == trajbound)]
            print("select components", components, "and unique components", data["column"].unique())
            print("select trajbound", trajbound, "and unique trajbounds", data["trajbound"].unique())
            # Create the subplot
            for component in components:
                if skipripple and "ripple" in component:
                    continue
                component_data = data[data["column"] == component].sort_values(by="lindist_bin_mid")
                # Add the base value to the bootstrap_mean if stratification is active
                # if stratify_spectral_components and component.startswith(("S1", "S2","U", "V")):
                component_data[field] += component_fill_bases[component]
                component_data[field] *= scale
                # Plot the curve
                sns.lineplot(x="lindist_bin_mid", y=field, data=component_data, color=component_colors[component], 
                             errorbar="se", 
                             ax=axes[i, j], linestyle=line_styles[component])
                axes[i,j].axhline(y=scale*(component_fill_bases[component]+0.5), color=component_colors[component], linestyle="dotted")
                if ci:
                    ci_99 = component_data.groupby("lindist_bin_mid").quantile(0.99)[field]
                    ci_1 = component_data.groupby("lindist_bin_mid").quantile(0.01)[field]
                    axes[i,j].plot(ci_99.index, ci_99, color=component_colors[component], linestyle="dashed")
                    axes[i,j].plot(ci_1.index, ci_1, color=component_colors[component], linestyle="dashed")
                if i == 2 or i == 3:
                    axes[i,j].set_ylim([0.3, 1.4])
                # Fill under the curve based on the component type and the flag
                # if component in ["U1", "U2", "U3", "V1", "V2", "V3"] or (stratify_spectral_components and component.startswith(("S1", "S2"))):
                #     axes[i, j].fill_between(component_data["lindist_bin_mid"], 
                #                             component_fill_bases[component], 
                #                             component_data[field], 
                #                             color=component_colors[component], 
                #                             alpha=0.3)
                # Shade the confidence interval around the curve if the flag is set
                # if shade_confidence_intervals:
                #     component_data_bootstrap_samples = component_data[field].tolist()
                #     confidence_interval = np.percentile(component_data_bootstrap_samples, [2.5, 97.5])
                #     axes[i, j].fill_between(component_data["lindist_bin_mid"], confidence_interval[0], confidence_interval[1], color=component_colors[component], alpha=0.1)
            # Set the title and labels
            axes[i, j].set_title(f'Trajbound = {trajbound}')
            axes[i, j].set_xlabel("Lindist Bin Midpoint")
            axes[i, j].set_ylabel("Normalized Bootstrap Mean")
            # Rotate the x-axis labels for readability
            axes[i, j].tick_params(axis='x', rotation=45)
            # Add a legend
            axes[i, j].legend()

# Improve the layout
plt.tight_layout()
# Add space for the overall title
# fig.subplots_adjust(top=0.92)
plt.show()

figfolder = '/Volumes/MATLAB-Drive/Shared/figures/lindist_bootstrap/'
if not os.path.exists(figfolder):
    os.makedirs(figfolder)
plt.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals.png', dpi=300)
plt.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals.svg', dpi=300)
plt.savefig(figfolder + f'lindist_bootstrap{append}_{field}_balancedanim_collapseanimals.pdf', dpi=300)


# ------------------------------------------------------
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def plot_by_epoch(df, column, epochs, field="bootstrap_mean_smooth"):
    # Get the unique epochs
    unique_epochs = sorted(df[epochs].unique())
    
    # Create colormap based on the number of epochs
    colors = cm.rainbow(np.linspace(0, 1, len(unique_epochs)))
    
    # Create a dictionary to map epoch to color
    epoch_color_map = dict(zip(unique_epochs, colors))
    
    # Create a subplot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Group data by epoch and plot each group
    for epoch, group in df.groupby(epochs):
        # Subset to the column of interest
        group = group[group["column"] == column]
        
        # Sort values for consistent plotting
        group = group.sort_values(by="lindist_bin_mid")
        
        # Plot the line for this epoch
        sns.lineplot(x="lindist_bin_mid", y=field, data=group, ax=ax, color=epoch_color_map[epoch], label=f'Epoch {epoch}')
    
    # Set title and labels
    ax.set_title(f'{column} over epochs')
    ax.set_xlabel('lindist_bin_mid')
    ax.set_ylabel(field)
    
    # Show the plot
    plt.show()

# Call the function
plot_by_epoch(animal_bootstrap_means_combined, "S1theta", "epoch")


