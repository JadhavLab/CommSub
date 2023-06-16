function locationDistribution = calculateLocationDistribution(X_pos, Y_pos, grid_size, components)
% this function calculate the distrubution of all pattern priniciple
% components across the track given a grid size 
%
% INPUT:
% X_pos/Y_pos: the position of the animal at different timePoints
% grid_size: size of the grid to partition the arena
% components: the 1*T matrix of a component's strength throughout time
%
% OUTPUT:
% locationDistribution: the distribution of component strength throughout time

%% divide the arena into grid 
X_coordinate = min(X_pos):grid_size:max(X_pos);
Y_coordinate = min(Y_pos):grid_size:max(Y_pos);
locationDistribution = nan(numel(X_coordinate)-1, numel(Y_coordinate)-1);
%% add the component stength within each position grid
for i = 1:numel(X_pos)
    % find which location the animal is at
    for j = 1:numel(X_coordinate)-1
        if X_pos(i) > X_coordinate(j) && X_pos(i) <= X_coordinate(j+1)
            for k = 1:numel(Y_coordinate)-1
                if Y_pos(i) > Y_coordinate(k) && Y_pos(i) <= Y_coordinate(k+1)
                    if isnan(locationDistribution(j,k))
                        locationDistribution(j,k) = 0;
                    end
                    % add the sum of the components of all patterns???
                    locationDistribution(j,k) = locationDistribution(j,k) + components(i);
                end
            end
        end
    end
end
end

