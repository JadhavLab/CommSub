figFolder = fullfile(figuredefine, 'dimensionRemoval');
if ~exist(figFolder)
    mkdir(figFolder);
end
skws = {'export_path',figFolder};

[groups, genHs] = findgroups(rt.genHshort);
uGroups = unique(groups);
for G = uGroups'

    RT = rt(groups == G,:);
    figc("dimensionRemoval perPattern" + genHs(G));
    g = gramm(...
        'x',RT.dimensionRemoved, ...
        'marker', categorical(RT.targetArea),...
        'y', RT.performance, ...
        'color', categorical(RT.sameDirectionLabel),'subset',~RT.sameDirection);
    g = g.facet_grid([],categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar');
    g = g.stat_summary('geom','area');
    g = g.set_color_options('chroma',0,'lightness',30);
    g.set_text_options("interpreter",'latex');
    g = g.set_names('x','dims removed',...
        'y','peformance',...
        'column','Interaction',...
        'row','Pattern', ...
        'color','Brain area removed');
    snapnow;
    g=g.update('subset',RT.sameDirection);
    %g = g.facet_grid([],categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar');
    g = g.stat_summary('geom','area');
    g = g.set_color_options(); % Restore default color
    g.draw();
    warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
    sgtitle("")
    g.expoRT('file_name', get(gcf,'Name') + ".svg",skws{:})
end
