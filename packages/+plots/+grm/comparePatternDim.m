function plotCornerHist(T, genH_name)
%PLOTCORNERHIST Plot a corner histogram of the predictive dimensions of

    % Define subsets
    hpcsubset = T(T.directionality == "hpc-hpc" & T.genH_name == genH_name,:);
    pfcsubset = T(T.directionality == "hpc-pfc" & T.genH_name == genH_name,:);

    % Define variables for plotting
    x = hpcsubset.rrDim;
    y = pfcsubset.rrDim;

    % Calculate the maximum absolute value of x and y
    M = max([abs(x); abs(y)]);

    % Set up figure
    clf
    f = figure();
    corner_kws = {'edges', -M:(1/20 * M):M, 'aspect', 1, 'location', [M/2], 'fill', 'transparent', 'normalization', 'countdensity'};
    set(f, 'Position', get(0, 'ScreenSize'));

    % Create gramm plot
    g = gramm('x', x, 'y', y, 'subset', hpcsubset.control == "control");
    g.facet_grid(categorical(hpcsubset.patternAbstractSymbol), []);
    g.geom_jitter('alpha', 0.01, 'width', 0.35, 'height', 0.35);
    g.stat_cornerhist(corner_kws{:});
    g.set_point_options('base_size', 5);
    g.set_text_options('label_scaling', 1.5, 'base_size', 10);
    g.set_names('x', "HPC" + newline + "predictive dimensions", 'y', "PFC" + newline + "predictive dimensions", 'Color', 'Pattern', 'row', '', 'Lightness', 'Treatment/Control');
    g.set_color_options('chroma', 0);
    g.set_text_options('interpreter', 'latex', 'base_size', 10);

    % Update gramm plot
    g.update('subset', hpcsubset.control ~= "control", 'color', categorical(hpcsubset.patternAbstract));
    g.set_color_options();
    g.stat_cornerhist(corner_kws{:});
    g.geom_abline('style', 'k:');
    g.geom_jitter('alpha', 0.30, 'width', 0.35, 'height', 0.35);

    % Draw the plot
    g.draw();

    % Save the plot
    g.export("file_name", figuredefine("gramm", "patterndim", "compare-dir-" + henH_name + "_cornerhist"), "file_type", "pdf");
    g.export("file_name", figuredefine("gramm", "patterndim", "compare-dir-" + henH_name + "_cornerhist"), "file_type", "png");

end

