function [opt, optKws, optDict] = initialize(varargin)

% INTIALIZE initializes parameters for seqnmf analysis
debug = false;
% Parse optional arguments
% ------------------------
ip = inputParser;

% Data related
ip.addParameter('storeDatBeforeSubsamp', true);
ip.addParameter('subsample', 0.05);
ip.addParameter('animal_list', {'JS12','JS13','JS14','JS15','JS21','ER1','KL8'});
% Seqnmf related
ip.addParameter('maxiter', []);
ip.addParameter('epoch_type', 'all', @(x) ischar(x) || isstring(x));
ip.addParameter('timescale',  [], @isnumeric);
ip.addParameter('K',  []);
ip.addParameter('L',  nan);
ip.addParameter('lambda',  []); % 0.0009 decent
ip.addParameter('lambdaOrthoH',  []);
ip.addParameter('lambdaOrthoW',  []);
ip.addParameter('orthoName', "", @isstring);

% Field specific 
ip.addParameter('fieldstr',   [], @(x) ischar(x) || isstring(x));

% Specicic to this function
ip.addParameter('skipList', {}); %Variables to skip in the assignnennt process at theend of this function
ip.addParameter('assignInCaller', true) %Variables to skip in the assignnennt process at theend of this function
ip.addParameter('opt',[])
ip.parse(varargin{:})
opt = ip.Results;

% Parameters
% ----------
opt.onehotencoding = false;

switch opt.fieldstr
case 'wSPLFP'
    opt.fields = {'S1','S2','wpli', 'ca1SPpfcLFP_C', 'pfcSPca1LFP_C'};
case 'wPhi'
    opt.fields = {'S1','S2','wpli', 'phi'};
case 'PhiC'
    opt.fields = {'phi','C','ca1SPpfcLFP_C', 'pfcSPca1LFP_C'};
case 'wPhi'
    opt.fields = {'S1','S2','wpli', 'phi'};
case 'speccoh'
    opt.fields = {'S1','S2','wpli'};
case 'spectra'
    opt.fields = {'S1','S2'};
case 'traj'
    opt.fields = {'S1', 'S2', 'wpli', 'trajdist', 'aglindist', 'velocity', 'reward', 'trajbound'};
case 'trajNoWPLI'
    opt.fields = {'S1', 'S2', 'trajdist', 'velocity', 'reward'};
otherwise
    opt.fields = {'S1','S2','wpli','velocity'};
end

switch lower(opt.orthoName)
case "parts-based"
    opt.lambdaOrthoW = opt.lambda * 1.0;
    opt.lambda       = opt.lambda * 0.0;
case "events-based"
    opt.lambdaOrthoH = opt.lambda * 1.0;
    opt.lambda       = opt.lambda * 0.0;
end

% SeqNMF and fieldpack args and assign in caller
% ----------------------------------------------
fieldpack_kws = {...
    'zscore', opt.fields(:),...
...'minmax', fields(3),...
    'floorceil', [-1,2.5]...
    };
kws  = {'lambda', opt.lambda, ...
        'lambdaOrthoW', opt.lambdaOrthoW, ...
        'lambdaOrthoH', opt.lambdaOrthoH, ...
        'maxiter', opt.maxiter, ...
        'K', opt.K};
opt.kws = kws;
opt.fieldpack_kws = fieldpack_kws;

seqnmf_kws = containers.Map('KeyType','char','ValueType','any');
for k = 1:2:numel(kws)
    seqnmf_kws(kws{k}) = kws{k+1};
end

if opt.assignInCaller
    for name = who()'
        if isequal(name{1}, 'opt') || ismember(name, opt.skipList)
            continue
        end
        assignin('caller', name{1}, eval(name{1}));
    end
end


% Setup folder
% ------------
%if fieldstr
%    folder = ['seqnmf_' fieldstr];
%    disp(['seqnmf with ' fieldstr]);
%end
%folder = [folder '_' epoch_type];
%folder = [folder '_' num2str(timescale)];
%folder = [folder '_lambda=' sprintf('%1.0e',lambda)];
%folder

% Function
% --------
seqnmf_factorize = @seqNMF_gpu; % use faster gpu accelerated version

addpath(genpath('~/Code/Src_Matlab/ry_Utility'))
for animal = opt.animal_list
    animalToPath(animal{1});
end
 

%% Aliases
%% -------
%lowhigh_defintions = [0.2, 0.8]; % percentile definitions for high and low respectively
%CA1                = 1; PFC   = 2; CA1PFC = 3;
%theta              = 1; delta = 2;
%high               = 2; low   = 1;
%frequency          = 2; time  = 1;
%
%% Windowing
%% ---------
%bands      = [6 12; 0.5 4]; % 1- theta, 2 - delta
%samprate   = 0.1;
%win        = [10, 10]; % 10 seconds before and after
%plot_range = [0, 40];
