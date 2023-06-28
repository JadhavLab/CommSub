% Converts mlx to m files if mod times are newer or file doesn't exist

cd(fullfile(codedefine, 'Notebooks'));
m_folder = fullfile(codedefine, 'Notebooks', 'mfiles');
if ~exist(m_folder,'dir')
    mkdir(m_folder);
end
mlx_files = dir('*.mlx');
for i = 1:length(mlx_files)
    mlx_file = mlx_files(i).name;
    m_file = [mlx_file(1:end-4) '.m'];
    if ~exist(m_file,'file') || ...
            datenum(mlx_files(i).date) > datenum(dir(m_file).date)
        export(mlx_file, fullfile(m_folder, m_file));
    end
end
