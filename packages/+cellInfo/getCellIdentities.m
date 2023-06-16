function celllookup = getCellIdentities(animal, cell_index, areaPerNeuron)

cellinfo = ndb.load(animal,'cellinfo');
inds = ndb.indicesMatrixForm(cellinfo);


cellinfo = ndb.toNd(cellinfo);

numCells = size(cell_index, 1);

regionalIndex = zeros(1,numCells);
regionalCsi = zeros(1,numCells);
regionalSpikeWidth = zeros(1,numCells);
regionalMeanrate = zeros(1,numCells);
hpcIndex = 0;
pfcIndex = 0;
%%
for i = 1:numCells
    all_csi = [];
    all_spikewidth = [];
    all_meanrate = [];
    if areaPerNeuron(i) == "CA1"
        hpcIndex = hpcIndex+1;
        regionalIndex(i) = hpcIndex;
    elseif areaPerNeuron(i) == "PFC"
        pfcIndex = pfcIndex+1;
        regionalIndex(i) = pfcIndex;
    else
        keyboard
    end
    
    for j = 1:size(cellinfo, 2)
        if ~ismember([j, cell_index(i,:)], inds(:,2:end), 'rows')
            continue
        end
        
        info = cellinfo(1,j,cell_index(i,1), cell_index(i,2));
        
        if ~isempty(info.csi)
            all_csi = [all_csi, info.csi];
        end
        if ~isempty(info.meanrate)
            all_meanrate = [all_meanrate, info.meanrate];
        end
        
        if ~isempty(info.spikewidth)
            all_spikewidth = [all_spikewidth, info.spikewidth];
        end
        
    end
    
    regionalCsi(i) = median(all_csi);
    regionalMeanrate(i) = nanmedian(all_meanrate);
    regionalSpikeWidth(i) = nanmedian(all_spikewidth);
    
end

%% categorize the cells
cellType = string.empty;
for i = 1:numCells
    if regionalMeanrate(i) >= 10 && regionalSpikeWidth(i)<0.4
        cellType(i) = "int";
    else
        cellType(i) = "pyr";
    end
end

%% make table
celllookup = table(areaPerNeuron', cell_index, regionalIndex', regionalMeanrate', regionalSpikeWidth', cellType');
celllookup.Properties.VariableNames = {'region' 'cell_index' 'index_by_region' 'Meanrate' 'Spikewidth' 'cellType'};
end