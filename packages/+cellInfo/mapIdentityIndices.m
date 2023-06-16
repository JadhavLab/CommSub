function regional_indexed_cells = mapIdentityIndices(cells, areaPerNeuron)

% maps the indices to the respective cell identities, as a field of the
% cellinfo struct


nCells = numel(cells);
pfcIndex = 0;
hpcIndex = 0;

for i = 1:nCells
    currCell = cells{i};
    currArea = areaPerNeuron(i);
    
    currCell.area = currArea;
    
    if currArea == "PFC"
        pfcIndex = pfcIndex+1;
        currCell.regionalIndex = pfcIndex;
    else
        hpcIndex = hpcIndex+1;
        currCell.regionalIndex = hpcIndex;
    end
    cells{i} = currCell;
end

regional_indexed_cells = cells;
for i =1:nCells
    disp(regional_indexed_cells{i}.regionalIndex)
end
end

