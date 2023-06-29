SO = SimAnalTableOverall;
SC = SimAnalTableComponent;
SO.pattern = SO.pattern1 + " - " + SO.pattern2;
SC.pattern = SC.pattern1 + " - " + SC.pattern2;
SO.control = SO.pattern.count('control');
SC.control = SC.pattern.count('control');
SO.patternSimple = SO.pattern;
SO.patternSimple(SO.control > 0) = "control";
SC.patternSimple = SC.pattern;
SC.patternSimple(SC.control > 0) = "control";

SC(SC.corrMax == 0,:) = [];
SO(SO.corrMaxMean == 0,:) = [];
SO(isnan(SO.corrMean),:) = [];

%{
% Overall
%}
% figure;
% g = gramm( 'x', categorical(SO.patternSimple), 'y', SO.corrMean, 'color', categorical(SO.patternSimple))
% g = g.facet_wrap(categorical(SO.field), categorical(SO.directionality))
% g= g.stat_summary('geom','bar', 'width', 2, 'type','95percentile')
% g= g.stat_summary('geom','black_errorbar', 'type','95percentile')
% g.draw()
% 
% figure;
% g = gramm( 'x', categorical(SO.patternSimple), 'y', SO.corrMean, 'color', categorical(SO.patternSimple))
% g = g.facet_wrap(categorical(SO.field), categorical(SO.directionality))
% g= g.stat_summary('geom','bar', 'width', 2)
% g= g.stat_summary('geom','black_errorbar')
% g.draw()

%{
% Comparing the top two components of each pattern combo
%}
dim = 2;
for field = ["U","V","UV"]
    figure;
    filt = SC.component <= dim & SC.field == field;
    g(1,1) = gramm( 'x', categorical(SC.patternSimple), 'y', SC.corrMax, 'color', categorical(SC.patternSimple), 'subset', filt);
    g(1,1) = g(1,1).facet_wrap(categorical(SC.directionality));
    g(1,1)= g(1,1).stat_summary('geom','bar', 'width', 2, 'type','bootci');
    g(1,1)= g(1,1).stat_summary('geom','black_errorbar',  'type','bootci');
    g(1,1).set_title('Top two pattern max similarity');
    g(1,1) = g(1,1).axe_property('ylim',[0.3,0.7],'xticklabelrotation',45);
    filt = SC.component <= 2& SC.field == field;
    g(2,1) = gramm( 'x', categorical(SC.patternSimple), 'y', SC.corrMin, 'color', categorical(SC.patternSimple), 'subset', filt);
    g(2,1) = g(2,1).facet_wrap(categorical(SC.directionality));
    g(2,1)= g(2,1).stat_summary('geom','bar', 'width', 2, 'type','bootci');
    g(2,1)= g(2,1).stat_summary('geom','black_errorbar',  'type','bootci');
    g(2,1) = g(2,1).axe_property('ylim',[-0.6,-0.3],'xticklabelrotation',45);
    g(2,1).set_title('Top two pattern max dissimilarity')
    g.set_names('x','Pattern-comparison');
    g(1,1).set_names('y',"Maximum Similarity," + newline + "Top 2 components");
    g(2,1).set_names('y',"Maximum Dissimilarity, " + newline + "Top 2 components");
    g.set_title("Field = " + field)
    g.draw()
end

%{
% Comparing the top two components of each pattern combo
%}
dim = 2;
violin_kws={'npoints',500, 'bandwidth',0.03};
for field = ["U","V","UV"]
    figure;
    filt = SC.component <= dim & SC.field == field;
    g(1,1) = gramm( 'x', categorical(SC.patternSimple), 'y', SC.corrMax, 'color', categorical(SC.patternSimple), 'subset', filt);
    g(1,1) = g(1,1).facet_wrap(categorical(SC.directionality));
    g(1,1)= g(1,1).stat_violin('width',1.5, violin_kws{:});
    g(1,1).set_title('Top two pattern max similarity');
    g(1,1) = g(1,1).axe_property('ylim',[0,1.2],'xticklabelrotation',45);
    filt = SC.component <= dim & SC.field==field;
    g(2,1) = gramm( 'x', categorical(SC.patternSimple), 'y', SC.corrMin, 'color', categorical(SC.patternSimple), 'subset', filt);
    g(2,1) = g(2,1).facet_wrap(categorical(SC.directionality));
    g(2,1)= g(2,1).stat_violin('width',1.5, violin_kws{:});
    g(2,1) = g(2,1).axe_property('ylim',[-1.2,-0],'xticklabelrotation',45);
    g(2,1).set_title('Top two pattern max dissimilarity')
    g.set_names('x','Pattern-comparison');
    g(1,1).set_names('y',"Maximum Similarity," + newline + "Top " + dim + " components");
    g(2,1).set_names('y',"Maximum Dissimilarity, " + newline + "Top " + dim + " components");
    g.set_title("Field = " + field)
    g.draw()
end
