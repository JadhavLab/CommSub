function T = pairwiseCorr(Spk, source_area, target_area)
    % Get indices for the source and target areas
    source = Spk.areaPerNeuron == source_area;
    target = Spk.areaPerNeuron == target_area;

    % Get the spike count matrices for the source and target areas
    overall_source = Spk.spikeCountMatrix(source,:);
    overall_target = Spk.spikeCountMatrix(target,:);

    if source_area == "CA1"
        source_area = "hpc";
    elseif source_area == "PFC"
        source_area = "pfc";
    end
    if target_area == "CA1"
        target_area = "hpc";
    elseif target_area == "PFC"
        target_area = "pfc";
    end

    % Get the cell arrays of the source and target sampled during different network states
    X_source = Spk.(source_area).X;
    X_target = Spk.(target_area).X;

    % Calculate the overall correlation matrix
    % [~,allzerosSource] = clean.throwOutAllZeros(overall_source');
    overall_corr = corrcoef([overall_source', overall_target']);
    source_indices = 1:size(overall_source,1);
    target_indices = size(overall_source,1)+1:size(overall_source,1)+size(overall_target,1);
    overall_corr = overall_corr(source_indices, target_indices);
    sa = source_area;
    ta = target_area;

    % Initialize the table
    VariableNames = {'source_index', char("target_index"), 'correlation_value', 'overall_corr_value', 'diff', 'source_area','target_area','pattern', 'diff_norm', 'diff_norm2', 'diff_norm3','other_value'};
    T = cell(1, numel(X_source));

    for p = progress(1:numel(X_source), 'Title', 'Calculating pairwise correlations')
        compare = mod(p+3, 6) + 1;
        C = corrcoef([X_source{p}', X_target{p}']);
        c = C(source_indices, target_indices);
        control = corrcoef([X_source{compare}', X_target{compare}']);
        c_control = control(source_indices, target_indices);
        % keyboard
        % figure;tiledlayout(2,1);nexttile;imagesc(c);colorbar;title('c');nexttile;imagesc(overall_corr);colorbar;title('overall_corr');axis(findobj('type','axes'),'square')
        % Calculate the differences and add to the table
        diff       = c - overall_corr;
        diff_norm  = (c- overall_corr)./overall_corr;
        diff_norm2 = (c- overall_corr)./(overall_corr + c);
        diff_norm3 = c./overall_corr;
        [source_idx, target_idx] = ind2sub(size(c), 1:numel(c));
        source_area = repmat(sa, numel(source_idx), 1);
        target_area = repmat(ta, numel(target_idx), 1);
        pattern     = repmat(p,  numel(source_idx), 1);
        t = table(source_idx', target_idx', c(:), overall_corr(:), diff(:), source_area, target_area, pattern, diff_norm(:), diff_norm2(:), diff_norm3(:), control_value(:)...
        'VariableNames', VariableNames);
        T{p} = t; 
    end
    T = vertcat(T{:});

end

