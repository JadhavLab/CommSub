function directory = animalToPath(animal, varargin)
% Adds an animals folder information into the path for this matlab session.
%
% animal, char or cell array
%   set of animals
% exclude_eeg_folder, logical
%   whether to exclude EEG folders

exclude_eeg_folder = false;
varargin = optlistassign(who,varargin{:});
assert(isempty(varargin), 'Unidentified optional arg');

if ischar(animal)
  animal = {animal};
end

for a = animal
  info = animaldef(a{1});
  directory = info{2};

  directoryset = genpath(directory);
  if exclude_eeg_folder
      folds = split(directory,':');
      hasEEG = contains(folds,'EEG');
      folds = string(folds(~hasEEG));
      directoryset = join(folds,':');
  end

  % Add all folders in the project folder into the path
  addpath(directoryset);

  % Remove anything related to a .git repository in this folder
  rmpath( genpath( fullfile(directory,'.git') ) );
end

