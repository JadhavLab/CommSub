function [performance, loss, Yhat] = getUncorrelatedPerformance(B_,X_source, X_target,...
                                                    optDim, numDimsUsedForPrediction, cvLoss)

% this function removes the predictive dimensions by their rank in
% principle components, and returns the predictive performance subsequently
% after removing each prominent dimension

% B_slice is the prediction matrix with the optDim rank PLUS the intercept

M = B_' * cov(X_source');

[U,D,V] = svd(M);

[~, nCol] = size(D);

numSingular = 0;
for i = 1:nCol
    if sum(D(:,i))==0
        numSingular = numSingular+1;
    end
end
 
Q = V(:,nCol-numSingular+1:end);

unCorrelated = M*Q;

% check if the matrix is almost zero/less than eplison

Xhat = X_source' * Q;

optDimLambdaReducedRankRegress = ModelSelect...
	(cvLoss, numDimsUsedForPrediction);

B = RidgeRegress(X_target', Xhat, optDimLambdaReducedRankRegress);
[loss, Yhat] = RegressPredict(X_target', Xhat, B);

performance = 1 - loss;

end




