function T = splitLFPfactors(T, dim, nType)

nFrequencyByType =  size(T, dim);
new_dims = [nFrequencyByType/nType, nType];

T_size = size(T);
new_T_size = [T_size(1:dim-1) new_dims T_size(dim+1:end)];

T = reshape(T, new_T_size);
