function dimensionRemoved(rt, genH)
% DIMENSIONREMOVED 

rtnew  = rt(rt.method == genH,:);
assert(~isempty(rtnew), "No data for method " + genH);

figure;
g = gramm(...
    'x', rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color',     categorical(rtnew.removePattern), ...
    'lightness', categorical(rtnew.sameDirectionLabel),...
    'linestyle', categorical(rtnew.sameDirectionLabel),...
    'subset', ~contains(rtnew.removePattern, "control") & ...
              ~contains(rtnew.basePattern,   "control"));
g = g.facet_grid(categorical(rtnew.basePatternLabel), ...
                 categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','type','sem','dodge',0.5);
g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g = g.set_names('x','dims removed',...
                'y','peformance',...
                'column','Interaction',...
                'row','Pattern', ...
                'color','Pattern dims removed',...
                'lightness',"Remove Same/Different"+ newline+"Pred. Target",...
                'linestyle',"Remove Same/Different"+ newline+"Pred. Target");
g.axe_property('XLim',[0 6]);
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval",genH + ": remove dims from same+diff - zoom"),'file_type',["svg","png"]);
poststeps(g)

%%  REMOVE PATTERN FROM SAME TARGET AREA ONLY
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',rtnew.sameDirection);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5);
% g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of similar target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval",genH + ": remove dims from same target area only"),'file_type',["pdf","svg"]);
poststeps(g)

% REMOVE PATTERN FROM SAME TARGET AREA ONLY, zoom into [0,6]
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',rtnew.sameDirection);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type','sem');
% g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of similar target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
set(g.facet_axes_handles,'xlim',[0 6])
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + ": remove dims from same target area only - zoom"),'file_type',["pdf","svg"]);
poststeps(g)

%% REMOVE PATTERN FROM DIFFERENT TARGET AREA
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',~rtnew.sameDirection);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g = g.set_color_options('lightness',200);
g.set_text_options("interpreter",'latex');
g = g.set_names('x', 'dims removed',...
                'y', 'peformance',...
                'column', 'Interaction',...
                'row', 'Pattern', ...
                'color', 'Pattern dims removed');
g=g.set_title('Removing different only');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + ": remove dims from different target area only"),'file_type',["pdf","svg"]);
poststeps(g)

%% REMOVE PATTERN FROM DIFFERENT TARGET AREA, zoom into [0,6]
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.removePattern), ...
    'subset',~rtnew.sameDirection);
g = g.facet_grid(categorical(rtnew.basePatternLabel), categorical(rt.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar','dodge',0.5,'type','sem');
% g = g.stat_summary('geom','area');
g.set_text_options("interpreter",'latex');
g=g.set_title("Removing dimensions" + newline + "(of different target area only)");
g = g.set_names('x','dims removed',...
    'y','peformance',...
    'column','Interaction',...
    'row','Pattern', ...
    'color','Pattern dims removed');
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export('file_name',figuredefine("dimensionRemoval", genH + ": remove dims from different target area only - zoom"),'file_type',["pdf","png"]);
poststeps(g)

%% Area summary
figure;
g = gramm(...
    'x',rtnew.dimensionRemoved, ...
    'marker', rtnew.targetArea,...
    'y', rtnew.performance, ...
    'color', categorical(rtnew.sameDirectionLabel),'subset',~rt.sameDirection);
g = g.facet_grid([],categorical(rtnew.targetArea));
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
g=g.update('subset',rtnew.sameDirection);
%g = g.facet_grid([],categorical(rtnew.targetArea));
g = g.stat_summary('geom','line');
g = g.stat_summary('geom','point');
g = g.stat_summary('geom','errorbar');
g = g.stat_summary('geom','area');
g = g.set_color_options(); % Restore default color
set(g.facet_axes_handles, 'xlim', [0 6])
g.draw();
warning off; set(g.facet_axes_handles, 'yscale','log'); warning on;
steps()

g.export("file_name", figuredefine("dimensionRemoval/", genH + ": area-summary"), "file_type", ["pdf","png"]);
poststeps(g)

    function steps()
        sgtitle(genH);
        set(gcf, 'Position',  get(0, 'Screensize'));
    end
    function poststeps(g)
        set(gcf, 'Position',  get(0, 'Screensize'));
        g.export('file_name',... 
        figuredefine("dimensionRemoval","all"),'file_type',["pdf","svg"]);
    end

end
