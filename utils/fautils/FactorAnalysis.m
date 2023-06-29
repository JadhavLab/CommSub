function [L, psi, logLike] = FactorAnalysis(S, q, varargin)
%
% [L, psi] = FactorAnalysis(S, q) applies Factor Analysis to the sample
% covariance matrix S, using a latent state with q dimensions.
%
%   p: data dimensionality
%   q: latent dimensionality
%   N: number of data points
%
% INPUTS: 
%
% S - sample covariance matrix (p x p)
% q - number of factors
%
% OUTPUTS:
%
% L       - factor loadings (p x q)
% psi     - diagonal of uniqueness matrix (p x 1)
% logLike - log likelihood at final EM iteration
%
% OPTIONAL ARGUMENTS (NAME-VALUE PAIRS):
%
% 'Method' - 'FA' (default) or 'PPCA'
%
% Code adapted from fastfa.m by Byron Yu.
%
% @ 2018 Joao Semedo -- joao.d.semedo@gmail.com


C_TOL = 1e-8;	      %	Stopping criterion for EM
C_MAX_ITER = 1e8;     %	Maximum number of EM iterations
C_MIN_FRAC_VAR = .01; % Fraction of overall data variance for each
%						observed dimension to set as the private variance
%						floor. This is used to combat Heywood cases, where
%						ML parameter learning returns one or more zero
%						private variances
%						(See Martin & McDonald, Psychometrika, 1975).

method = 'FA';
r = [];
A = [];
kws = {};
for i = 1:2:numel(varargin)
	switch upper(varargin{i})
		
	case 'METHOD'
		method = varargin{i+1};
        
    case 'RANK'
		r = varargin{i+1};
    case 'ORIG'
		A = varargin{i+1};
    case 'GPU'
        kws = {'gpuArray'};
	end
end

% Excludes neurons with 0 variability from the analysis.
s = diag(S);
%idxs = find( abs(s) < sqrt( eps(class(s)) ) );
idxs = abs(s) < sqrt( eps(class(gather(s))) );
if any(idxs)
	auxP = size(S, 1);
	S(idxs,:) = [];
	S(:,idxs) = [];
	auxIdxs = (1:auxP)';
	auxIdxs(idxs) = [];
end

p = size(S, 1);
if isempty(r)
    r = rank(S);
end
if r == p
	scale = exp( 2*sum( log(diag( chol(S) )) )/p );
else
    if ~isempty(A)
        [~,sigma,~] = svd(A,'econ');
        eigenvalues = diag(sigma).^2;
    else
        eigenvalues = eig(S);
    end
	d = sort( eigenvalues, 'descend' );
    % potential speedup in accepted answer here: https://www.mathworks.com/matlabcentral/answers/262276-how-can-i-increase-the
	scale = geomean( d(1:r) );
end
L = randn(p, q, kws{:})*sqrt(scale/q);
psi = diag(S);

varFloor = C_MIN_FRAC_VAR*psi;

I = eye(q, kws{:});
c = -p/2*log(2*pi);
logLike = 0;

for i = 1:C_MAX_ITER
	
	% ==================================================================
	%   E-step
	% ==================================================================
	
	% C = L*L' + Psi is the estimated sample covariance matrix. inv(C)
	% is computed using the Woodbury identity to speed up computation
	% time, inverting a q-by-q matrix instead of directly inverting C,
	% which is p-by-p (q << p).
	invPsi = diag( 1./psi );
	invPsiTimesL = invPsi*L;
	invC = invPsi - invPsiTimesL/( I + L'*invPsiTimesL )*L'*invPsi;
	
	% EZ = (X - M)*inv(C)*L;
	% EZZ = 1/n*(n*I - n*L'*inv(C)*L + EZ'*EZ)
	%     = I - L'*inv(C)*L + L'*inv(C)*S*inv(C)*L
	V = invC*L;
	StimesV = S*V;
	EZZ = I - V'*L + V'*StimesV;
	
	% ==================================================================
	%   Compute log-likelihood
	% ==================================================================
	
	prevLogLike = logLike;
    try
        ldm = sum( log(diag( chol(invC) )) );
    catch ME
        %save('errordata_doNotDeleteMePlease.mat','-v7.3');
        %keyboard
        %  Ryan Y addeed  this: if the matrix is not positive semi definite,
        %  nearestSPD  will find the closest semidef matrix.
        warning('Cholskey not positive semidef. Finding nearest.')
        invC = nearestSPD(invC);
        ldm = sum( log(diag( chol(invC) )) );
    end
	logLike = c + ldm - .5*sum( sum(invC.*S) );
	
	if i <= 2
		baseLogLike = logLike;
	elseif (logLike-baseLogLike) < (1+C_TOL)*(prevLogLike-baseLogLike)
		break
	end
	
	% ==================================================================
	%   M-step
	% ==================================================================
	
	% L = ((X - M)'*EZ)*inv(EZZ)
	L = StimesV/EZZ;
	
	% Psi = diag( S - 1/n*((X - M)'*EZ)*L' )
	psi = diag(S) - sum( StimesV .* L, 2 );
	switch upper(method)
	case 'PPCA'
		psi = mean(psi) * ones(p, 1, kws{:});
	case 'FA'
		psi = max(varFloor, psi);
	end
	
end

if any(idxs)
	
	auxL = L;
	L = zeros(auxP, q, kws{:});
	L(auxIdxs,:) = auxL;
	
	auxPsi = psi;
	psi = zeros(auxP, 1, kws{:});
	psi(auxIdxs) = auxPsi;
	
end

end
