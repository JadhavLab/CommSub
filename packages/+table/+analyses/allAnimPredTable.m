function [combinedPatternsTable, singlePredTable] = allAnimPred(Patterns, Option)
    % Get the size of the original array                                                                                                                                                                                                                                        
    const = option.constants();
    regions = [const.HPC, const.PFC];

    sz = size(Patterns);
    % Calculate the product of the dimensions to collapse                                                                                                                                                                                                                       
    collapsed_size = prod(sz(1:2));                                                                                                                                                                                                                                            
    % Reshape the array
    ReshapedPatterns = squeeze(reshape(Patterns, [collapsed_size, sz(end-2:end)]));
    % Preallocate growing vectors
    combinedPatternsTable   = table();
    singlePredTable         = table();
    patternPerformanceTable = table();
    for p = progress(1:size(ReshapedPatterns, 1), 'Title', 'Animal*GenH*Partition')
        for i = 1:Option(1).nPatternAndControl
            for d = regions
                % Existing code here
                % Get current nTarget and nSource
                currSource = ReshapedPatterns(p, d, i).X_source';
                currTarget = ReshapedPatterns(p, d, i).X_target';
                currB = ReshapedPatterns(p, d, i).rankRegress.B;
                nSource = repmat(size(currSource, 2), size(currTarget, 2), 1);
                nTarget = repmat(size(currTarget, 2), size(currTarget, 2), 1);
                iTarget = 1:nTarget;
                iSource = 1:nSource;
                animal = repmat(ReshapedPatterns(p, d, i).animal, nTarget(1), 1);
                genH   = repmat(ReshapedPatterns(p, d, i).genH_name, nTarget(1), 1);
                direction = repmat(ReshapedPatterns(p, d, i).directionality, nTarget(1), 1);
                name = repmat(ReshapedPatterns(p, d, i).name, nTarget(1), 1);
                [perf, mu, dev] = plots.calculatePredictionPerformance(currSource, currTarget, currB);
                mu = repmat(mu, nTarget(1), 1);
                dev = repmat(dev, nTarget(1), 1);
                iTarget = iTarget';
                perf = perf';
                newtab = table(animal, genH, name, direction, nSource, nTarget, iTarget, perf, mu, dev);
                combinedPatternsTable = [combinedPatternsTable; newtab];
                for iSource = 1:nSource(1)
                    if ~isfield(ReshapedPatterns(p,d,i), 'rankRegress') || ...
                        ~isfield(ReshapedPatterns(p,d,i).rankRegress, 'singlesource_B')
                        continue
                    end
                    curr_singleB = ReshapedPatterns(p,d,i).rankRegress.singlesource_B;
                    if isempty(curr_singleB)
                        continue
                    end
                    iSource;
                    curr_singleB = curr_singleB{iSource};
                    curr_singlesource = currSource(:,iSource);
                    [perf, mu, dev] = plots.calculatePredictionPerformance(curr_singlesource, currTarget, curr_singleB);
                    mu = repmat(mu, nTarget(1), 1);
                    dev = repmat(dev, nTarget(1), 1);
                    iSource = repmat(iSource, nTarget(1), 1);
                    iTarget = iTarget(:);
                    perf = perf(:);
                    newtab = table(animal, genH, direction, nSource, nTarget, iSource, iTarget, perf, mu, dev);
                    singlePredTable = [singlePredTable; newtab];
                end
            end
        end
    end
    if ~exist(figuredefine("tables"), 'dir'), mkdir(figuredefine("tables")); end
    writetable(combinedPatternsTable, fullfile(figuredefine("tables"), 'fig2_prediction.csv'));
    writetable(singlePredTable, fullfile(figuredefine("tables"), 'fig2_singlePrediction.csv'));

    % Just to not break legacy downstream code
    combinedPatternsTable.Var8 = combinedPatternsTable.perf;
    combinedPatternsTable.Var7 = combinedPatternsTable.iTarget;
    singlePredTable.Var9 = singlePredTable.perf;
    singlePredTable.Var8 = singlePredTable.iTarget;
