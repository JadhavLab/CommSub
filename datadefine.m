function folder = datadefine()
% Automates pointing to the data  folders location on our computers

if ispc
    disp('Defining data folder for Ziyi machine');
    folder = " C:\Users\BrainMaker\commsubspace\SingleDayExpt";
elseif ismac
    disp('Defining data folder for Ryan''s machine');
    %folder = '~/Data/commsubspace/';
    folder(1) = "/Volumes/sharespace-commsub/data";
    folder(2) = "~/Data/commsubspace";
elseif isunix
    disp('Defining data folder for Ryan''s linux machine');
    folder = '~/Data/Raw/SingleDayExpt';
end
