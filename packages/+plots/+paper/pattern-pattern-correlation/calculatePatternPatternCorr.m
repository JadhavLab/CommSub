%Patterns = Patterns_AllAnimals;
figFolder = fullfile(figuredefine, 'calcPatternPatternCorr');
if ~exist(figFolder)
    mkdir(figFolder);
end
skws = {'export_path',figFolder};

Opt.fischerZ = true;
Opt.fischerZRaw = false;
Opt.upperval = 0.99;
Opt.lowerval = 0.01;
Opt.truemin = false;
Opt.controlAvg = true; % use only control period for rmean?
Opt.normalize = true;

nPatterns = numel(patternNames);
nGen = numel(X_hpc);
fish = @(x)1/2 * log((1+x)./(1-x));

% Determine pattern average
% -------------------------
hpc_mean = spikeRateMatrix(areaPerNeuron == "CA1",:)';
pfc_mean = spikeRateMatrix(areaPerNeuron == "PFC",:)';
if ~Opt.controlAvg
    Rmean = corrcoef([hpc_mean,pfc_mean]);
else
    clear Rmean
    for iGen = 1:numel(Option.generateH)
        for iPat = 4:6
            rtmp{iPat-3} = [X_hpc{1}{iPat}', X_pfc{1}{iPat}'];
        end
        Rmean{iGen} = cat(1, rtmp{:});
        Rmean{iGen} = corrcoef(Rmean{iGen});
    end
end

% Shortcut for cross versus same interactions
% -------------------------------------------
same_hpc = {1:size(hpc_mean, 2), 1:size(hpc_mean, 2)};
same_pfc = {size(hpc_mean,2) + 1:size(pfc_mean, 2), size(hpc_mean,2) + 1:size(pfc_mean, 2)};
different = {same_hpc{1}, same_pfc{1}};

% Shortcut for cross versus same interactions
fig("Pattern Correlation Matrices")
AX = tiledlayout(nGen,nPatterns);
Rset={};
genH = shortcut.generateH(string(Option.generateH));
for iGen = 1:nGen
    for iPat = 1:nPatterns
        nexttile;
        X = [X_hpc{iGen}{iPat}', X_pfc{iGen}{iPat}'];
        X = zscore(X,[],1);
        R = corrcoef(X);
        if Opt.fischerZRaw
            R = fish(R);
        end
        if Opt.controlAvg
            mu = Rmean{iGen};
        else
            mu = Rmean;
        end
        if Opt.normalize
            if Opt.fischerZ
                R  = fish(R) - fish(mu);
            else
                R = R - mu;
            end
        end
        hold off;
        imagesc(R)
        cmocean('balance')
        set(gca,'clim',[-1,1])
        colorbar
        Rset{iGen, iPat} = R;
        Pat{iGen,iPat} =  iPat * ones(size(R));
        GenH{iGen,iPat} = iGen * ones(size(R));
        title("GenH = " + genH(iGen) +newline + " Pat = " + shortcut.patternSymbols(patternNames(iPat)));
    end
end
% Normalize colors
ravel = @(x) x(:);
notInf = @(x) x(~isinf(x));
maxVal = cellfun(@(x) quantile(ravel(notInf(x)), Opt.upperval), get(findobj(gcf,'type','Image'),'Cdata'));
minVal =  cellfun(@(x) quantile(ravel(notInf(x)), Opt.lowerval), get(findobj(gcf,'type','Image'),'Cdata'));
maxVal = max(maxVal);
if Opt.truemin
    minVal = min(minVal);
else
    minVal = min(minVal);
    minVal = -max(abs(minVal), abs(maxVal));
    assert(isscalar(minVal));
end
set(findobj(gcf,'type','axes'),'clim',[minVal,maxVal])


% Plot mean
fig("Mean R, all times")
if Opt.controlAvg
    for iGen = 1:numel(Option.generateH)
        nexttile
        imagesc(Rmean{iGen});
        cmocean('balance');
        set(gca,'clim',[minVal,maxVal]);
        colorbar;
    end
else
    imagesc(Rmean);
    cmocean('balance');
    set(gca,'clim',[minVal,maxVal]);
    colorbar;
end

%%------------------------
%% total interaction terms
%%------------------------
C = @(x) categorical(x);
stitch =  @(x) cat(1, x{:});
homo   = @(y) stitch(cellfun(@(x) [ravel(x(same_hpc{:})); ravel(x(same_pfc{:}))],...
    y, 'UniformOutput', false));
hetero = @(y) stitch(cellfun(@(x) ravel(x(different{:})),...
    y, 'UniformOutput', false));
kws = {'VariableNames', ["genH","pattern", "corr"]};
rt_homo   = table(homo(GenH), homo(Pat), homo(Rset),kws{:});
rt_homo.samediff = repmat("within brain area", height(rt_homo), 1);
rt_hetero =  table(hetero(GenH), hetero(Pat), hetero(Rset),kws{:});
rt_hetero.samediff = repmat("between brain area", height(rt_hetero), 1);
rt = [rt_homo;rt_hetero];
rt.patternLevel = repmat("",height(rt),1);
rt(rt.pattern > 3,:).patternLevel = repmat("low",sum(rt.pattern>3),1);
rt(rt.pattern <= 3,:).patternLevel = repmat("high",sum(rt.pattern<=3),1);
rt.patternAbstract = mod(rt.pattern-1,3)+1;
genH = shortcut.generateH(string(Option.generateH));
rt.patternType = genH(rt.genH(:))';
rt.patternName = patternNames(rt.pattern)';
rt.patternNameAbs = patternNames(rt.patternAbstract)';

figc("CDF Distributions of homo/hetero interactions" +newline + "high and low")
g = gramm('x', rt.corr, 'color', C(rt.patternLevel), 'linestyle', C(rt.samediff), 'lightness', C(rt.samediff))
g.stat_density('function','cdf')
g.set_title(""+newline + newline +"Summary: Patterns vs control potentially time of increased interaction")
g.draw()
g.export('file_name', 'cdf-overall',skws{:})

figc("CDF Distributions of homo/hetero interactions, split genH x pat" +newline + "high and low")
g = gramm('x', rt.corr, 'color', C(rt.patternLevel), 'linestyle', C(rt.samediff), 'lightness', C(rt.samediff))
g.facet_grid(C(rt.patternType), C(rt.patternNameAbs))
g.stat_density('function','cdf','npoints',1000)
g.set_title(""+newline + newline +"Summary: Patterns vs control potentially time of increased interaction")
g.draw()
g.export('file_name', 'cdf-split by pattern',skws{:})

figc("Distributions of homo/hetero interactions" +newline + "high and low")
g = gramm('x', rt.corr, 'color', C(rt.patternLevel), 'linestyle', C(rt.samediff), 'lightness', C(rt.samediff))
g.stat_density()
g.set_title(""+newline + newline +"Summary: Patterns vs control potentially time of increased interaction")
g.axe_property('xlim',[-0.3,0.3]);
g.draw()
g.export('file_name', 'overall',skws{:})

figc("Distributions of homo/hetero interactions, split genH x pat" +newline + "high and low")
g = gramm('x', rt.corr, 'color', C(rt.patternLevel), 'linestyle', C(rt.samediff), 'lightness', C(rt.samediff))
g.facet_grid(C(rt.patternType), C(rt.patternNameAbs))
g.stat_density('function','pdf','npoints',1000)
g.set_title(""+newline + newline +"Summary: Patterns vs control potentially time of increased interaction")
g.axe_property('xlim',[-0.3,0.3]);
g.draw()
g.export('file_name', 'split by pattern',skws{:})

fig("log Distributions of homo/hetero interactions, split genH x pat" +newline + "high and low")
g = gramm('x', rt.corr, 'color', C(rt.patternLevel), 'linestyle', C(rt.samediff), 'lightness', C(rt.samediff))
g.facet_grid(C(rt.patternType), C(rt.patternNameAbs))
g.stat_bin('normalization','countdensity')
g.set_title(""+newline + newline +"Summary: Patterns vs control potentially time of increased interaction")
g.draw()


