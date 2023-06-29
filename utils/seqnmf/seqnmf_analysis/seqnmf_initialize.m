

% Parameters
% ----------
animal_list = {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'};
%fields      = {'S1','S2','wpli', 'phi'};
fields      = {'S1','S2','wpli'};
timescale = 10;

% Aliases
% -------
lowhigh_defintions = [0.2, 0.8]; % percentile definitions for high and low respectively
CA1                = 1; PFC   = 2; CA1PFC = 3;
theta              = 1; delta = 2;
high               = 2; low   = 1;
frequency          = 2; time  = 1;

% Windowing
% ---------
bands      = [6 12; 0.5 4]; % 1- theta, 2 - delta
samprate   = 0.1;
win        = [10, 10]; % 10 seconds before and after
plot_range = [0, 40];

% Level of analysis
% -----------------
if ~exist('nmfs', 'var')
    nmfs = {'overall'}; % modes : animal_overall (animal-wise overall patterns), overall (patterns for entire set of animals at once)
end

% SeqNMF and fieldpack args
% -------------------------
fieldpack_kws = {...
    'zscore', fields(1:2),...
    'minmax', fields(3),...
    'floorceil', [-1,2]...
    };
seqnmf_arg_names = ["", "Parts-based", "Events-based"];
seqnmf_kws  = {{'lambda', 0.0009, 'lambdaOrthoW', 0, 'lambdaOrthoH', 0, 'maxiter', 1.2e3, 'K', 4},...
               {'lambda', 0.0009, 'lambdaOrthoW', 0, 'lambdaOrthoH', 1, 'maxiter', 1.2e3, 'K', 4},...
               {'lambda', 0.0009, 'lambdaOrthoW', 1, 'lambdaOrthoH', 0, 'maxiter', 1.2e3, 'K', 4}};

% Setup folder
% ------------
if ismember('phi', fields)
    folder = 'seqnmf_wPhi';
else
    folder = 'seqnmf';
end
folder = [folder '_' num2str(timescale)];

% Function
% --------
seqnmf_factorize = @seqNMF_gpu;
