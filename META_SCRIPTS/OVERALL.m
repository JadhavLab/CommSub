py = pyenv("Version", "~/miniconda3/bin/python")
disp(pyversion); % Check python version

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
SemedoPaperFigures_m

% -----------------------------
% Grammer of graphics summaries
% -----------------------------
GrammPlots

% -----------------------------
% Subspace angle
% -----------------------------
plots.subspace.angle.Run;
% Python plots
py.importlib.import_module(codedefine("Notebooks","python","create_clustergram.py"));
py.importlib.import_module(codedefine("Notebooks","python","create_graph.py"));

% --------------
% Python figures
% --------------
script1 = ...
py.importlib.import_module(codedefine("Notebooks","python","dimPred.py"));
script2 = ...
py.importlib.import_module(codedefine("Notebooks","python","prediction.py"));
script3 = ...
py.importlib.import_module(codedefine("Notebooks","python","ccatime.py"));
% ----------------------------------------------

