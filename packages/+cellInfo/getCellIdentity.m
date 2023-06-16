function all_unique_cells = getCellIdentity(animal, cell_index)

% this function returns the list of cell properties that aligns with later
% numerical indices (1:numCells) used in marking source/target neurons

cellinfo = ndb.load(animal,'cellinfo');
indices = ndb.indicesMatrixForm(cellinfo);

%%
total_unique_indices = unique(cell_index,'rows');
num_cell = size(total_unique_indices,1);
%%
unique_indices = [];
places_of_unique_indices = [];
% m = containers.Map;
% for i = 1:size(indices,1)
for i = 1:size(cell_index,1)
    currIndex = indices(i, 3:4);
    
    unique = false;
    if size(places_of_unique_indices) == 0
        unique = true;
    end
    
    for j = 1:size(places_of_unique_indices,2)
%         disp(size(places_of_unique_indices,2)  )
        if currIndex(1) == unique_indices(j,1) && currIndex(2) == unique_indices(j,2)
            unique = false;
            break;
        end
        unique = true;
    end
    
    if unique
        unique_indices = [unique_indices; currIndex];
        places_of_unique_indices = [places_of_unique_indices, i];
    end
    
end
%% put all the unique cells into an array of cells

all_unique_cells = cell(1,num_cell);
for i = 1:num_cell
    currcell = ndb.get(cellinfo, indices(places_of_unique_indices(i),:));
    try
        all_unique_cells{i} = currcell;
    catch
        keyboard
    end
end

end

