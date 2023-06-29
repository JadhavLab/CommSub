function seqnmf_save(file, data)

for field = fieldnames(data)'
    if isequal(field{1},'data')
        continue
    end
    field, 
    try
        assign(field{1},data.(field{1})); 
        save([file filesep field{1}], field{1}); 
    catch ME
        keyboard
        warning('Skipping %s\n', field{1}); 
    end
end

if isfield(data,'data')
    disp('Saving data.mat')
    data = data.data;
    save([file filesep 'data'], 'data');
end

