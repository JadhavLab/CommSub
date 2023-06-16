function [sol, model, u, v] = UVtensorDecomp(t, b_, optD, rank_setting, varargin)
% Uses structured data fusion to find a tensor decomposition of 
% the time-varying tensor commsubspace, optionally constraining the
% cellular dimensions to the U, V vectors of the SVD of b_.

ip = inputParser;
ip.addParameter('nonneg', ["time","trial"]); % Impose non-negativity?
ip.addParameter('constrainUV', false); % Constrain to b_ UV svd?j
ip.parse(varargin{:})
Opt = ip.Results;

if isstruct(t)
    t = ful(t);
end
if isempty(rank_setting)
    rank_setting = min(size(b_));
end

% Obtain uv components of regression
[u,~,v] = svd(b_);

% Build a structured data fusion model
model = struct;

% Initialize variables for decomp
model.variables.time = randn(size(t,4), rank_setting);
model.variables.trial = randn(size(t,3), rank_setting);
model.variables.area_x = u(:,1:rank_setting);
model.variables.area_y = v(:,1:rank_setting);

% Make fusion aware that these are factors
model.factors.T   = {'time'};
model.factors.t  = {'trial'};
model.factors.X = {'area_x'};
model.factors.Y = {'area_y'};

% Apply constraints?
if Opt.constrainUV
    for factor_name = ["X","Y"]
        model.factors.(factor_name) = {model.factors.(factor_name){:}, @(x) struct_const(x,[])};
    end
end
if ~isempty(Opt.nonneg) && ~isequal(Opt.nonneg,0) && ~isequal(Opt.nonneg, false)
    factors = string(fieldnames(model.factors))';
    if ~islogical(Opt.nonneg)
        if isnumeric(Opt.nonneg)
            factors = factors(Opt.nonneg);
        else
            factors = factors(factors == Opt.nonneg);
        end
    end
    for factor_name = factors
        model.factors.(factor_name) = {model.factors.(factor_name){:}, @struct_nonneg};
    end
end

% What to factorize
model.factorizations.tensor.data = t;
model.factorizations.tensor.cpd  = {'X','Y','t','T'};
%if ~Opt.constrainUV
model.factorizations.matrix.data = b_;
model.factorizations.matrix.cpd  = {'X','Y'};
%end

% Check the model's correctness
sdf_check(model,'print');

% Find the solution
tic;
sol = sdf_nls(model, 'Display', 5, 'TolX', 1e-5);
time = toc;
disp("Model took " + time + " seconds");
