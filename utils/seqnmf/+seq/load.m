function data = seqnmf_load(file)

data = struct();
for fieldfile = dir([file filesep '*.mat'])'
    [~,field] = fileparts(fieldfile.name);
    field, 
    try
        tmp = load([string(fieldfile.folder) + filesep + string(fieldfile.name)]);
        data.(field) = tmp.(field);
    catch ME
        warning('Skipping %s\n', field); 
    end
end

