import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from scipy.stats import combine_pvalues
import scipy

# Reload the data

animal = 'ZT2'
typ    = 'power' # coherence
folder = '/Volumes/MATLAB-Drive/Shared/figures/tables/'
df = pd.read_csv(f'/{folder}/{animal}{typ}_granger_causality_directionality.csv')
# If two copies of same column, drop one
df = df.loc[:, ~df.columns.duplicated()]

def get_significant_interactions(df, max_lag_order=df['lag'].max(),
                                 drop_columns=['time', 'vel']):
    # Filter the data to only include rows with the maximum lag order
    df_max_lag = df[df['lag'] == max_lag_order]
    # Group by 'column1', 'column2' and 'direction'
    grouped_max_lag = df_max_lag.groupby(['column1', 'column2', 'direction'])
    # Combine p-values and calculate average magnitude
    def combine_and_average(group):
        combined_pvalue, _ = combine_pvalues(group['pvalue'], method='fisher')
        avg_magnitude = group['magnitude'].mean()
        return pd.Series([combined_pvalue, avg_magnitude, combined_pvalue < 0.05], index=['combined_pvalue', 'avg_magnitude', 'is_significant'])
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
    # ..........................................................................
    # Drop the 'time' and 'vel' columns
    # ..........................................................................
    def flexdrop(df, values, axis):
        columnsindices = df.columns if axis == 1 else df.index
        for value in values:
            if value in columnsindices:
                df = df.drop(value, axis=axis)
        return df
    heatmap_data_signed_log = flexdrop(heatmap_data_signed_log, drop_columns, 1)
    heatmap_data = flexdrop(heatmap_data, drop_columns, 1)
    heatmap_data_signed_log = flexdrop(heatmap_data_signed_log, drop_columns, 0)
    heatmap_data = flexdrop(heatmap_data, drop_columns, 0)
    return significant_max_lag, heatmap_data, heatmap_data_signed_log

sml, hd, hds = dict(), dict(), dict()
for lag in range(1, df['lag'].max()):
    significant_max_lag, heatmap_data, heatmap_data_signed_log = get_significant_interactions(df, max_lag_order=lag)
    print(f'lag: {lag}')
    print(f'number of significant interactions: {len(significant_max_lag)}')
    print(f'number of columns: {len(heatmap_data_signed_log.columns)}')
    print(f'number of rows: {len(heatmap_data_signed_log.index)}')
    sml[lag] = significant_max_lag
    hd[lag] = heatmap_data
    hds[lag] = heatmap_data_signed_log

# --------------------------------------------------------------------------
# Granger causality heatmap
# --------------------------------------------------------------------------
# Set up the matplotlib figure
def plot_heatmap(heatmap_data_signed_log):
    f, ax = plt.subplots(figsize=(12, 9))

    # Draw the heatmap
    cmap = sns.diverging_palette(250, 10, s=80, l=55, n=9, as_cmap=True)
    sns.heatmap(heatmap_data_signed_log, cmap=cmap, center=0, annot=True, fmt=".2f",
                square=True, linewidths=.5, cbar_kws={"shrink": .5})

    # Set the title
    ax.set_title('Signed Logarithm of Average Magnitude of Influence for Maximum Lag Order (excluding \'time\' column)')

    # Show the plot
    plt.show()

plot_heatmap(hds[1])
plot_heatmap(hds[2])
plot_heatmap(hds[3])
plot_heatmap(hds[4])
plot_heatmap(hds[5])
plot_heatmap(hds[6])

# --------------------------------------------------------------------------
# Granger causality heatmap with zones
# --------------------------------------------------------------------------

# Define the order of the zones
zones = ['Cavg', 'S1', 'S2', 'U', 'V', 'behavior']

# Define the mapping of column names to zones
zone_mapping = {'Cavg': ['Cavgtheta', 'Cavgdelta', 'Cavgripple'],
                'S1': ['S1theta', 'S1delta', 'S1ripple'],
                'S2': ['S2theta', 'S2delta', 'S2ripple'],
                'U': ['U1', 'U2', 'U3'],
                'V': ['V1', 'V2', 'V3'],
                'behavior': ['vel', 'accel', 'lindist']}

# Define the order of the columns based on the zones
column_order = [column for zone in zones for column in sorted(heatmap_data_signed_log.columns) if column in zone_mapping[zone]]

# Reorder the columns
heatmap_data_signed_log = heatmap_data_signed_log[column_order]

# Create a list to hold the lines that separate the zones
zone_lines = []

# Calculate the position of the lines
for i in range(1, len(zones)):
    zone_lines.append(heatmap_data_signed_log.columns.get_loc(zone_mapping[zones[i]][0]) - 0.5)

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(14, 11))

# Draw the heatmap
sns.heatmap(heatmap_data_signed_log, cmap=cmap, center=0, annot=True, fmt=".2f",
            square=True, linewidths=.5, cbar_kws={"shrink": .5}, linecolor='black')

# Draw the lines that separate the zones
for line in zone_lines:
    ax.axhline(line + 0.5, color='black', linewidth=3)
    ax.axvline(line + 0.5, color='black', linewidth=3)

# Set the title
ax.set_title('Signed Logarithm of Average Magnitude of Influence for Maximum Lag Order (excluding \'time\' column)')

# Show the plot
plt.show()

# --------------------------------------------------------------------------

def plot_siglag_barplot(sml, lag):
    # Check the significant_max_lag dataframe for significant interactions for bullet 2 and bullet 4
    significant_max_lag = sml[lag].reset_index()
    cca_to_spectral_col1_original = significant_max_lag[
            (significant_max_lag['column1'].str.contains('U|V') & significant_max_lag['column2'].str.contains('S1|S2|Cavg'))]
    cca_to_spectral_col2_original = significant_max_lag[
            (significant_max_lag['column2'].str.contains('U|V') & significant_max_lag['column1'].str.contains('S1|S2|Cavg'))]

    # Display the results
    cca_to_spectral_col1_original, cca_to_spectral_col2_original
    # # Print the dataframes to check their contents
    # print(spectral_to_cca_col2)
    # print(cca_to_spectral_col2_original)

    # Create a 2x2 subplot
    fig, axs = plt.subplots(2, 2, figsize=(14, 16))

    # Create a barplot for spectral things to CCA (column1 to column2)
    sns.barplot(data=spectral_to_cca_col1, x='column1', y='avg_magnitude', hue='column2', ax=axs[0, 0], ci=None)
    axs[0, 0].set_title('Spectral to CCA (column1 to column2)')
    axs[0, 0].set_xlabel('Spectral Things (S1, S2, Cavg)')
    axs[0, 0].set_ylabel('Average Magnitude of Influence')
    axs[0, 0].tick_params(axis='x', rotation=45)

    # Create a barplot for CCA to spectral things (column1 to column2)
    if not cca_to_spectral_col1_original.empty:
        sns.barplot(data=cca_to_spectral_col1_original, x='column1', y='avg_magnitude', hue='column2', ax=axs[0, 1], ci=None)
        axs[0, 1].set_title('CCA to Spectral (column1 to column2)')
        axs[0, 1].set_xlabel('CCA Communication Subspace (U, V)')
        axs[0, 1].set_ylabel('Average Magnitude of Influence')
        axs[0, 1].tick_params(axis='x', rotation=45)
    else:
        axs[0, 1].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Create a barplot for spectral things to CCA (column2 to column1)
    if not spectral_to_cca_col2.empty:
        sns.barplot(data=spectral_to_cca_col2, x='column2', y='avg_magnitude', hue='column1', ax=axs[1, 0], ci=None)
        axs[1, 0].set_title('Spectral to CCA (column2 to column1)')
        axs[1, 0].set_xlabel('Spectral Things (S1, S2, Cavg)')
        axs[1, 0].set_ylabel('Average Magnitude of Influence')
        axs[1, 0].tick_params(axis='x', rotation=45)
    else:
        axs[1, 0].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Create a barplot for CCA to spectral things (column2 to column1)
    if not cca_to_spectral_col2_original.empty:
        sns.barplot(data=cca_to_spectral_col2_original, x='column2', y='avg_magnitude', hue='column1', ax=axs[1, 1], ci=None)
        axs[1, 1].set_title('CCA to Spectral (column2 to column1)')
        axs[1, 1].set_xlabel('CCA Communication Subspace (U, V)')
        axs[1, 1].set_ylabel('Average Magnitude of Influence')
        axs[1, 1].tick_params(axis='x', rotation=45)
    else:
        axs[1, 1].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Show the plots
    plt.tight_layout()
    plt.show()

plot_siglag_barplot(sml, 1)

# --------------------------------------------------------------------------


# Filter significant_max_lag for cases where spectral things influence CCA communication subspace components and vice versa
spectral_to_cca = significant_max_lag[(significant_max_lag['column1'].str.contains('S1|S2|Cavg') & significant_max_lag['column2'].str.contains('U|V')) |
                                      (significant_max_lag['column1'].str.contains('U|V') & significant_max_lag['column2'].str.contains('S1|S2|Cavg'))]

# Display the results
spectral_to_cca

# --------------------------------------------------------------------------

def plot_interactions_for_lag(df, lag):
    # Filter the dataframe for the given lag
    df_lag = df[df['lag'] == lag]

    # Calculate the combined p-value and average magnitude of influence for each pair of columns
    grouped_df_lag = df_lag.groupby(['column1', 'column2']).apply(lambda group: pd.Series({
        'combined_pvalue': scipy.stats.combine_pvalues(group['pvalue'], method='fisher')[1],
        'avg_magnitude': group['F'].mean()
    })).reset_index()

    # Filter for significant results
    significant_df_lag = grouped_df_lag[grouped_df_lag['combined_pvalue'] < 0.05]

    # Separate the data into four categories
    spectral_to_cca_col1 = significant_df_lag[significant_df_lag['column1'].str.contains('S1|S2|Cavg') & significant_df_lag['column2'].str.contains('U|V')]
    cca_to_spectral_col1 = significant_df_lag[significant_df_lag['column1'].str.contains('U|V') & significant_df_lag['column2'].str.contains('S1|S2|Cavg')]
    spectral_to_cca_col2 = significant_df_lag[significant_df_lag['column2'].str.contains('S1|S2|Cavg') & significant_df_lag['column1'].str.contains('U|V')]
    cca_to_spectral_col2 = significant_df_lag[significant_df_lag['column2'].str.contains('U|V') & significant_df_lag['column1'].str.contains('S1|S2|Cavg')]

    # Create a 2x2 subplot
    fig, axs = plt.subplots(2, 2, figsize=(14, 16))

    # Create a barplot for spectral things to CCA (column1 to column2)
    if not spectral_to_cca_col1.empty:
        sns.barplot(data=spectral_to_cca_col1, x='column1', y='avg_magnitude', hue='column2', ax=axs[0, 0], ci=None)
        axs[0, 0].set_title('Spectral to CCA (column1 to column2)')
        axs[0, 0].set_xlabel('Spectral Things (S1, S2, Cavg)')
        axs[0, 0].set_ylabel('Average Magnitude of Influence')
        axs[0, 0].tick_params(axis='x', rotation=45)
    else:
        axs[0, 0].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Create a barplot for CCA to spectral things (column1 to column2)
    if not cca_to_spectral_col1.empty:
        sns.barplot(data=cca_to_spectral_col1, x='column1', y='avg_magnitude', hue='column2', ax=axs[0, 1], ci=None)
        axs[0, 1].set_title('CCA to Spectral (column1 to column2)')
        axs[0, 1].set_xlabel('CCA Communication Subspace (U, V)')
        axs[0, 1].set_ylabel('Average Magnitude of Influence')
        axs[0, 1].tick_params(axis='x', rotation=45)
    else:
        axs[0, 1].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Create a barplot for spectral things to CCA (column2 to column1)
    if not spectral_to_cca_col2.empty:
        sns.barplot(data=spectral_to_cca_col2, x='column2', y='avg_magnitude', hue='column1', ax=axs[1, 0], ci=None)
        axs[1, 0].set_title('Spectral to CCA (column2 to column1)')
        axs[1, 0].set_xlabel('Spectral Things (S1, S2, Cavg)')
        axs[1, 0].set_ylabel('Average Magnitude of Influence')
        axs[1, 0].tick_params(axis='x', rotation=45)
    else:
        axs[1, 0].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    # Create a barplot for CCA to spectral things (column2 to column1)
    if not cca_to_spectral_col2.empty:
        sns.barplot(data=cca_to_spectral_col2, x='column2', y='avg_magnitude', hue='column1', ax=axs[1, 1], ci=None)
        axs[1, 1].set_title('CCA to Spectral (column2 to column1)')
        axs[1, 1].set_xlabel('CCA Communication Subspace (U, V)')
        axs[1, 1].set_ylabel('Average Magnitude of Influence')
        axs[1, 1].tick_params(axis='x', rotation=45)
    else:
        axs[1, 1].text(0.5, 0.5, 'No significant interactions', ha='center', va='center')

    plt.tight_layout()
    plt.show()

plot_interactions_for_lag(df, 1)
plot_interactions_for_lag(df, 2)

