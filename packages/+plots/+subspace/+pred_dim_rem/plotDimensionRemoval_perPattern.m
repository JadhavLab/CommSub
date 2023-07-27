figFolder = fullfile(figuredefine, 'dimensionRemoval');
if ~exist(figFolder)
    mkdir(figFolder);
end
skws = {'export_path',figFolder};

[groups, genHs] = findgroups(rt.genH);
uGroups = unique(groups);
for g = uGroups'

    figc("dimensionRemoval perPatternRegion " + genHs(g));
    RT = rt(groups == g, :);
    g = gramm(...
        'x',RT.dimensionRemoved, ...
        'y', RT.performance, ...
        'color', categorical(RT.removePattern), ...
        'lightness', categorical(RT.sameDirectionLabel),...
        'linestyle', categorical(RT.sameDirectionLabel));
    g = g.facet_grid(categorical(RT.basePatternLabel), categorical(RT.targetArea));
    g = g.stat_summary('geom','line');
    g = g.stat_summary('geom','point');
    g = g.stat_summary('geom','errorbar');
    g = g.stat_summary('geom','area');
    g.set_text_options("interpreter",'latex');
    g = g.set_names('x','dims removed',...
        'y','peformance',...
        'column','Interaction',...
        'row','Pattern', ...
        'color','Pattern dims removed',...
        'lightness',"Remove Same/Different"+newline+"Pred. Target",...
        'linestyle',"Remove Same/Different"+ newline+"Pred. Target");
    g.draw();
    warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
    sgtitle("")
    g.export('file_name', get(gcf,'Name') + ".svg",skws{:})
    g.export('file_name', get(gcf,'Name') + ".png",skws{:})

end
