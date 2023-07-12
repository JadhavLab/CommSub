function [B_rrr, B_slice] = getReducedB_(B, V, nSource, nTarget, optDim)
% this function returns the reduced B_ from RRR results?
% to be used to remove the preditive dimensions

assert(numel(nTarget)==1)

B_slice = B(2:nSource+1, (optDim-1)*nTarget+1:optDim*nTarget);
B_rrr = B_slice/V(:,1:optDim)';

end
