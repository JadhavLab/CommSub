% Converts m to mlx files if mod times are newer or file doesn't exist

cd(fullfile(codedefine, 'Notebooks'));
m_folder = fullfile(codedefine, 'Notebooks', 'mfiles');
if ~exist(m_folder,'dir')
    mkdir(m_folder);
end
m_files = dir(fullfile(m_folder, '*.m'));
for i = 1:length(m_files)
    m_file = m_files(i);
    m_file_path   = fullfile(m_folder, m_file.name);
    mlx_file_path = fullfile(codedefine, 'Notebooks', m_file.name(1:end-2) + "mlx");
    if ~exist(mlx_file_path,'file') || m_file.datenum > dir(mlx_file_path).datenum
        matlab.internal.liveeditor.openAndConvert(m_file_path, mlx_file_path);
    end
end
