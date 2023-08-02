import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

# Load the provided CSV file
df = pd.read_csv('/Volumes/MATLAB-Drive/Shared/figures/tables/eventuv.csv')
# Display the first few rows of the dataframe
df.head()
# Function to calculate the distance from a point to the line y = x
def distance_to_line(u, v):
    return np.abs(u - v) / np.sqrt(2)
# Function to calculate the Euclidean distance from a point to the origin
def distance_to_origin(u, v):
    return np.sqrt(u**2 + v**2)
# Function to scale a series to the range [0, 1]
def scale(series):
    scaler = MinMaxScaler()
    return scaler.fit_transform(series.values.reshape(-1, 1)).ravel()

# Remove rows with NaN in 'event_u_values' or 'event_v_values'
df_clean = df.dropna(subset=['event_u_values', 'event_v_values'])
# Calculate the distance from each point to the line y = x
df_clean['distance_to_line'] = distance_to_line(df_clean['event_u_values'], df_clean['event_v_values'])
# Calculate the Euclidean distance from each point to the origin
df_clean['distance_to_origin'] = distance_to_origin(df_clean['event_u_values'], df_clean['event_v_values'])
# Group the DataFrame and apply the scaling function to the 'distance_to_line' within each group
df_clean['scaled_distance_to_line']   = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal'])['distance_to_line'].transform(scale)
df_clean['scaled_distance_to_origin'] = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal'])['distance_to_origin'].transform(scale)
# Calculate 'on_commsub' as one minus the scaled distance to line
df_clean['on_commsub'] = 1 - df_clean['scaled_distance_to_line']
# Multiply the two distances to get the desired measure
df_clean['on_commsub_mag'] = df_clean['on_commsub'] * df_clean['scaled_distance_to_origin']
df_summary = df_clean.groupby(['genH', 'patterns', 'uv_components', 'animal', 'events']).mean().reset_index()

# Display the first few rows of the dataframe
df_clean.head()
uv_components = df_clean['uv_components'].unique()
patterns = df_clean['patterns'].unique()
genH_values = df_clean['genH'].unique()

# ------------- PLOT: Hist plot of on/off commsub -------------

# Create a subplot grid with one row for each uv_component and one column for each pattern
def plot_hist(df_clean, thing = 'on_commsub_mag'):
    fig, axs = plt.subplots(len(uv_components), len(patterns), figsize=(4*len(patterns), 4*len(uv_components)), sharey=True)
    # Adjust for case when there is only one row or column
    if len(uv_components) == 1:
        axs = [axs]
    if len(patterns) == 1:
        axs = [[ax] for ax in axs]
    
    q = df_clean[thing].quantile(0.99)
    q_half = df_clean[thing].mean()
    # For each uv_component
    for i, uv_component in enumerate(uv_components):
        # For each pattern
        for j, pattern in enumerate(patterns):
            ax = axs[i][j]
            # For each genH value
            for genH in genH_values:
                # Filter the DataFrame for the current uv_component, pattern, and genH
                df_sub = df_clean[(df_clean['uv_components'] == uv_component) & (df_clean['patterns'] == pattern) & (df_clean['genH'] == genH)]
                # Plot a histogram of the measure values
                sns.histplot(df_sub[thing], ax=ax, kde=False, label=genH, binrange=(0, 10))
            if i == 0:
                ax.set_title(f'Pattern {pattern}')
            if j == 0:
                ax.set_ylabel(f'UV Component {int(uv_component)}')
            ax.legend(title='genH')
            ax.set_xlim(0, q)
            ax.axvline(q_half, color='red')
    plt.tight_layout()
    plt.suptitle(f'{thing}')
    plt.show()

plot_hist(df_clean, thing = 'on_commsub_mag')
plot_hist(df_clean, thing = 'on_commsub')

# ------------- PLOT: Scatter plot of event_u_values vs event_v_values -------------


