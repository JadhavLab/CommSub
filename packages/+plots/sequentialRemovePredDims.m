function [performance, full_model] = sequentialRemovePredDims(X_source, X_target,...
    B_, optDim, cvLoss, numDimsUsedForPrediction, fromMostPredicitive, varargin)


% sequentially removing the most predicitive dimensions
% by taking away the largest eigen-vectores from the V matrix

ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.addParameter('do_plot', true, @islogical);
ip.addParameter('normalized', true, @islogical);
ip.addParameter('fromMostPredicitive', true, @islogical);

ip.parse(varargin{:});
opt = ip.Results;

%%

M = B_' * cov(X_source');

[U,~,V] = svd(M);

performance = zeros(1,optDim);

full_model = plots.plotPredictiveDimensions(numDimsUsedForPrediction,...
    cvLoss, "optDim", optDim, "mode", "rr",  "averaged", false, "do_plot", false);

for k= 1:optDim
    if ~fromMostPredicitive
        try
            V_removed = V(:, (optDim - (k-1)):end); % change
        catch
            keyboard
        end
    else
        V_removed = V(:, k:end);
    end
    
    V_removed = V(:, k:end);
    Xhat = X_source' * V_removed;
    
    optDimLambdaReducedRankRegress = ModelSelect...
        (cvLoss, numDimsUsedForPrediction);
    
    B = RidgeRegress(X_target', Xhat, optDimLambdaReducedRankRegress);
    [loss, ~] = RegressPredict(X_target', Xhat, B);
    
    performance(k) = 1 - loss;
    
    if opt.normalized
        performance(k)= performance(k)/full_model;
    end
end


% [B_slices] = takeSlicesOfB(B, optDim, withIntercept);
% performance = zeros(1,optDim);
%
% for i = 1:optDim
%     [performance(i),~] = plots.removePredictiveDimensions(B_, X_source',...
%                            X_target, numDimsUsedForPrediction, optDim, B_slices{i});
% end
%
% if opt.do_plot
%     xlabel('# predictive dimensions')
%     ylabel('Performance')
%
%     % plot where the optimal dimension falls
%     if ~isempty(opt.optDim)
%         lineObject=line([optDim,optDim],[0 full_model]);
%         lineObject.LineStyle = ':'; % Make line dotted
%         lineObject.LineWidth = 2;  % Thicken the line
%         lineObject.Color = opt.color; % Color it
%     end
% end

end

