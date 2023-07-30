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

# Number of bins and bootstrap samples
n_bins = 50
n_bootstrap_samples = 1000

# Update the columns to bootstrap
columns_to_bootstrap = ["U1", "U2", "U3", "V1", "V2", "V3", "Cavgtheta", 
                        "Cavgdelta", "Cavgripple", "S1theta", "S1delta", 
                        "S2theta", "wpli_avgtheta", "S1ripple", "S2ripple",
                        "S2delta", "wpli_avgripple", "wpli_avgdelta"]

# Reinitialize the DataFrame to hold the results
bootstrap_means_combined = []

# Rerun the bootstrap for each bin and each column, for each trajectory bound
for trajbound in tqdm([0, 1], desc="trajbound", total=2):
    # Filter the data based on trajbound
    data_trajbound = df[df["trajbound"] == trajbound]
    # Bin the lindist data
    data_trajbound["lindist_bin"] = pd.cut(data_trajbound["lindist"], bins=n_bins)
    
    # Loop over the unique combinations of column, trajbound, and lindist_bin
    iters = list(itertools.product(columns_to_bootstrap, [0, 1], data_trajbound["lindist_bin"].unique().categories))
    for column, trajbound, bin_label in tqdm(iters, desc="column, trajbound, lindist_bin", total=len(iters)):
        # Get the data for this column, trajbound, and bin_label
        data = data_trajbound.loc[(data_trajbound["lindist_bin"] == bin_label), ["animal",column]]
        
        # Find the minimum number of points available for each animal
        min_points_per_animal = data.groupby("animal").size().min()
        
        # Generate bootstrap samples
        for iboot in range(n_bootstrap_samples):
            # Initialize an empty list to hold the data for this bootstrap sample
            bootstrap_sample = []
            
            # Sample the minimum number of points from each animal's data
            for animal, animal_data in data.groupby("animal"):
                bootstrap_sample.append(animal_data.sample(min_points_per_animal, replace=True))
            
            # Concatenate the samples to form the bootstrap sample
            bootstrap_sample = pd.concat(bootstrap_sample)
            
            # Compute the mean of the bootstrap sample
            bootstrap_mean = bootstrap_sample.mean(numeric_only=True)
            
            # Add the result to the DataFrame
            bootstrap_means_combined.append({
                "iboot": iboot,
                "lindist_bin": bin_label,
                "column": column,
                "bootstrap_mean": bootstrap_mean,
                "trajbound": trajbound
            })


# Convert the list of dictionaries to a DataFrame
bootstrap_means_combined = pd.DataFrame(bootstrap_means_combined)
# Create a new lindist_bin_ind column
bootstrap_means_combined["lindist_bin_ind"] = bootstrap_means_combined["lindist_bin"].apply(lambda x: x.right)
bootstrap_means_combined.to_csv(os.path.join(folder, f'{name}_bootstrap{append}.csv'), index=False)

# Smooth
bootstrap_means_combined.sort_values(by=["column", "trajbound", "iboot", "lindist_bin_ind"], inplace=True)
bootstrap_means_combined["bootstrap_mean"] = bootstrap_means_combined.groupby(["column", "trajbound", "iboot"])["bootstrap_mean"].transform(lambda x: x.rolling(7, 1).mean())
bootstrap_means_combined["bootstrap_mean"] = bootstrap_means_combined.groupby(["column", "trajbound", "iboot"])["bootstrap_mean"].transform(lambda x: x.interpolate())


# ----------------------------------------------------
# Normalize the bootstrap_mean values
# ----------------------------------------------------
scaler = MinMaxScaler()
for column in tqdm(columns_to_bootstrap, desc="feature engineering", total=len(columns_to_bootstrap)):
    # Get the data for this column
    data = bootstrap_means_combined[bootstrap_means_combined["column"] == column]["bootstrap_mean"].values.reshape(-1, 1)
    # Scale the data
    scaled_data = scaler.fit_transform(data)
    # Update the DataFrame
    bootstrap_means_combined.loc[bootstrap_means_combined["column"] == column, "bootstrap_mean"] = scaled_data

bootstrap_means_combined.head()
bootstrap_means_combined.to_csv(os.path.join(folder, f'{name}_bootstrap_normalized{append}.csv'), index=False)
# ----------------------------------------------------

# Define the trajbounds for each column of the subplot grid
column_trajbounds = [0, 1]

# Update the components for each row of the subplot grid
row_components = [["U1", "U2", "U3"], ["V1", "V2", "V3"], ["Cavgtheta",
                                                           "S1theta",
                                                           "S2theta",
                                                           "wpli_avgtheta"],
                  ["Cavgdelta", "S1delta", "S2delta", "wpli_avgdelta"],
                  ["Cavgripple", "S1ripple", "S2ripple", "wpli_avgripple"]]


# Add the lindist_bin_mid column back to the DataFrame
bootstrap_means_combined["lindist_bin_mid"] = \
    bootstrap_means_combined["lindist_bin"].apply(lambda x: x.mid)

# # Create a 5x2 grid of subplots
# fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(15, 25))
#
# # Set the overall title
# fig.suptitle('Normalized Bootstrap Means for Different Components and Trajectories', fontsize=20)
#
# for i, components in enumerate(row_components):
#     for j, trajbound in enumerate(column_trajbounds):
#         # Get the data for this subplot
#         data = bootstrap_means_combined[
#             bootstrap_means_combined["column"].isin(components) &
#             (bootstrap_means_combined["trajbound"] == trajbound)
#         ]
#         # Create the subplot
#         sns.lineplot(x="lindist_bin_mid", y="bootstrap_mean", hue="column",
#                      data=data, ax=axes[i, j])
#         # Set the title and labels
#         axes[i, j].set_title(f'Trajbound = {trajbound}')
#         axes[i, j].set_xlabel("Lindist Bin Midpoint")
#         axes[i, j].set_ylabel("Normalized Bootstrap Mean")
#         # Rotate the x-axis labels for readability
#         axes[i, j].tick_params(axis='x', rotation=45)
#
# # Improve the layout
# plt.tight_layout()
# # Add space for the overall title
# fig.subplots_adjust(top=0.92)
plt.show()

# ----------------------------------------------------

# Set the flag for shading the confidence intervals to False
shade_confidence_intervals   = False
stratify_spectral_components = True  # Change this to False if you don't want to stratify and shade the spectral components

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

line_styles = {
    "U1": "-", "U2": "-", "U3": "-",
    "V1": "-", "V2": "-", "V3": "-",
    "Cavgtheta": "-", "Cavgdelta": "-", "Cavgripple": "-",
    "S1theta": "-", "S1delta": "-", "S1ripple": "-",
    "S2theta": "-", "S2delta": "-", "S2ripple": "-",
    "wpli_avgtheta": "--", "wpli_avgdelta": "--", "wpli_avgripple": "--"
}

# Create a 5x2 grid of subplots
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(15, 25))

# Set the overall title
fig.suptitle('Normalized Bootstrap Means for Different Components and Trajectories', fontsize=20)

for i, components in enumerate(row_components):
    for j, trajbound in enumerate(column_trajbounds):
        # Get the data for this subplot
        data = bootstrap_means_combined[
            bootstrap_means_combined["column"].isin(components) &
            (bootstrap_means_combined["trajbound"] == trajbound)
        ]
        # Create the subplot
        for component in components:
            component_data = data[data["column"] == component].sort_values(by="lindist_bin_mid")
            # Add the base value to the bootstrap_mean if stratification is active
            if stratify_spectral_components and component.startswith(("S1", "S2","U", "V")):
                component_data["bootstrap_mean"] += component_fill_bases[component]
            # Plot the curve
            axes[i, j].plot(component_data["lindist_bin_mid"], component_data["bootstrap_mean"], color=component_colors[component], label=component,
                            linestyle=line_styles[component])
            # Fill under the curve based on the component type and the flag
            if component in ["U1", "U2", "U3", "V1", "V2", "V3"] or (stratify_spectral_components and component.startswith(("S1", "S2"))):
                axes[i, j].fill_between(component_data["lindist_bin_mid"], component_fill_bases[component], component_data["bootstrap_mean"], color=component_colors[component], alpha=0.3)
            # Shade the confidence interval around the curve if the flag is set
            if shade_confidence_intervals:
                component_data_bootstrap_samples = component_data["bootstrap_mean"].tolist()
                confidence_interval = np.percentile(component_data_bootstrap_samples, [2.5, 97.5])
                axes[i, j].fill_between(component_data["lindist_bin_mid"], confidence_interval[0], confidence_interval[1], color=component_colors[component], alpha=0.1)
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
plt.savefig(figfolder + f'lindist_bootstrap{append}_balancedanim.png', dpi=300)
plt.savefig(figfolder + f'lindist_bootstrap{append}_balancedanim.svg', dpi=300)
plt.savefig(figfolder + f'lindist_bootstrap{append}_balancedanim.pdf', dpi=300)
