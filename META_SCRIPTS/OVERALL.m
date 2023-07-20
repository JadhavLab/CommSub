% HOUSES A FLOW TO RECREATE ANALYSES/FIGS

                                                
% o  /    ,---.          |                        
%   /     |---|,---.,---.|    ,   .,---.,---.,---.
%  /      |   ||   |,---||    |   |`---.|---'`---.
% /  o    `   '`   '`---^`---'`---|`---'`---'`---'
%                             `---'               

RunAll % Calls TheScript using various Option struct options
% RunAll_spectra
% RunAll_coh

                                       
% o  /    ,---.o                         
%   /     |__. .,---..   .,---.,---.,---.
%  /      |    ||   ||   ||    |---'`---.
% /  o    `    ``---|`---'`    `---'`---'
%               `---'                    
% Documentation of the data generated per animal
% (... more to come ... )

% -------------------------------------------------------
% Mostly knockoff figures from the inspired paper source 
% (some novel)
% -------------------------------------------------------
SemedoPaperFigures

% -----------------------------
% Grammer of graphics summaries
% -----------------------------
GrammPlots

% --------------
% Python figures
% --------------
py = pyenv("Version", "~/miniconda3/bin/python")
disp(pyversion); % Check python version
script1 = ...
py.importlib.import_module(codedefine("Notebooks","python","dimPred.py"))
script2 = ...
py.importlib.import_module(codedefine("Notebooks","python","prediction.py"))
% ----------------------------------------------
