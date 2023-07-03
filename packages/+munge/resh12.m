function Data = resh12(Data)
%ANIMPART Redo a partitioned data structure nPartitions x nDataset x ... x ...
%   so that top two dimensions are combined into one dimension
    tempData = squeeze(Data);
    newSize      = size(tempData);
    newSize      = [prod(newSize(1:2)), newSize(3:end)];
    Data = reshape(tempData, newSize);
end
