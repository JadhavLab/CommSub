function folder = codedefine(varargin)

persistent announced
if ~exist('announced', 'var') || isempty(announced)
    announced = false;
end

if ispc
    if ~announced
        disp('Defining code folder for Ziyi machine');
        announced = true;
    end
    folder = 'C:\\Users\BrainMaker\Matlab Drive\Shared\';
elseif ismac
    if ~announced
        disp('Defining code folder for Ryan''s machine');
        announced = true;
    end
    folder = '~/Data/Matlab-DRIVE/Shared/';
elseif isunix
    if ~announced
        disp('Defining code folder for Ryan''s linux machine');
        announced = true;
    end
    folder = '/Volumes/MATLAB-Drive/Shared/';
end

folder = string(folder);
if nargin > 0
    folder = fullfile(folder, varargin{:});
end

end
