function DataHigh_prep(raw_data)
% Load your raw data, assumed to be a 3D matrix [neurons x time x trials]


% Define a time bin size
binWidth = 15e-3;  % bin width in seconds

% Convert your raw data to the DataHigh struct format
D = struct();
for trial = 1:size(raw_data, 3)
    D(trial).data = raw_data(:, :, trial);
    D(trial).type = 'traj';  % assuming your data represents trajectories
end

% Use DataHigh to perform dimensionality reduction and visualize high-dimensional population activity
DataHigh(D, 'DimReduce');


