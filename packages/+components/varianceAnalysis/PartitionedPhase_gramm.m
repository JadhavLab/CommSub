
% Table columns
% animal epoch phase pattern component tperf val

savedir = '/Volumes/Colliculus/commsubspace/component';
mkdir(savedir);
filt_hpc = combined_hpc.pattern <= 3 & ~isnan(combined_hpc.val) & combined_hpc.tperf >= 0.4;
filt_pfc = combined_pfc.pattern <= 3 & ~isnan(combined_pfc.val) & combined_pfc.tperf >= 0.4;
combined_hpc = combined_hpc(filt_hpc,:);
combined_pfc = combined_pfc(filt_pfc,:);
combined_hpc.patternSym = patternSym(mod(combined_hpc.pattern-1,3)+1)';
combined_pfc.patternSym = patternSym(mod(combined_pfc.pattern-1,3)+1)';
phaseNames = ["reward", "error", "inBoundChoice","outBoundChoice"];
C=@categorical;
order =  unique(C(combined_hpc.patternSym));
order = order(end:-1:1);

% Absolute mean tperf
fig("Abs mean, tperf versus value"); clf
clear G;
G(1,2) = gramm('x', combined_hpc.tperf, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component);
G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
G(1,1) = gramm('x', combined_pfc.tperf, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component);
G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
arrayfun(@(g) g.stat_summary('type','fitnormalci','geom',{'bar','black_errorbar'}, 'bin_in', 6,'dodge',0, 'width', 0.9),G(:));
arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
G.set_color_options('lightness_range',[15 80]);
G.set_names('x','SS performance Bin','y','abs($\mu$)','column', 'Component');
G.set_order_options('color',order);
G.set_text_options('interpreter', 'latex');
G.set_title('Abs. mean over performance')
G.draw() 

% Absolute mean epoch
fig("Abs mean, epoch versus value"); clf
clear G;
G(1,2) = gramm('x', combined_hpc.epoch, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component);
G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
G(1,1) = gramm('x', combined_pfc.epoch, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component);
G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
arrayfun(@(g) g.stat_summary('type','fitnormalci','geom',{'bar','black_errorbar'}, 'dodge',0, 'width', 0.9),G(:));
arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
G.set_color_options('lightness_range',[15 80]);
G.set_names('x','Epoch','y','abs($\mu$)');
G.set_order_options('color',order);
G.set_text_options('interpreter', 'latex');
G.set_title('Abs. mean over Epochs')
G.draw()

%                                      
% .    ,          o                    
% |    |,---.,---..,---.,---.,---.,---.
%  \  / ,---||    |,---||   ||    |---'
%   `'  `---^`    ``---^`   '`---'`---'
                                     

varfun = @(y) [var(y); bootci(100, @var, y)];

% var tperf
fig("Var, tperf versus value"); clf
clear G;
G(1,2) = gramm('x', combined_hpc.tperf, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component);
G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
G(1,1) = gramm('x', combined_pfc.tperf, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component);
G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
arrayfun(@(g) g.stat_summary('type',varfun,'geom',{'bar','black_errorbar'}, 'bin_in', 6,'dodge',0, 'width', 0.9,),G(:));
%arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
G.set_color_options('lightness_range',[15 80]);
G.set_names('x','SS performance Bin','y','\sigma^2','column', 'Component');
G.set_order_options('color',order);
G.set_text_options('interpreter', 'latex');
G.set_title('Var over performance')
G.draw() 

% var epoch
fig("Var, Epoch versus value"); clg
clear G;
G(1,2) = gramm('x', combined_hpc.epoch, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component);
G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
G(1,1) = gramm('x', combined_pfc.epoch, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component);
G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
arrayfun(@(g) g.stat_summary('type',varfun,'geom',{'bar','black_errorbar'}, 'dodge',0, 'width', 0.9),G(:));
%arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
G.set_color_options('lightness_range',[15 80]);
G.set_names('x','Epoch','y','sigma^2');
G.set_order_options('color',order);
G.set_text_options('interpreter', 'latex');
G.set_title('Var over Epochs')
G.draw()


for phase = progress(1:numel(phaseNames))
    % Absolute mean tperf
    fig(phaseNames(phase) + " abs mean, tperf versus value"); clf
    clear G;
    G(1,2) = gramm('x', combined_hpc.tperf, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component, 'subset', combined_hpc.phase == phase);
    G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
    G(1,1) = gramm('x', combined_pfc.tperf, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component, 'subset', combined_pfc.phase == phase);
    G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
    arrayfun(@(g) g.stat_summary('type','fitnormalci','geom',{'bar','black_errorbar'}, 'bin_in', 6,'dodge',0, 'width', 0.9),G(:));
    arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
    G.set_color_options('lightness_range',[15 80]);
    G.set_names('x','SS performance Bin','y','abs($\mu$)','column', 'Component');
    G.set_order_options('color',order);
    G.set_text_options('interpreter', 'latex');
    G.set_title(phaseNames(phase) + ' abs. mean over performance')
    G(1).set_title('hpc-pfc')
    G(2).set_title('hpc-hpc')
    G.draw() 
end
util.fig.saveByName(savedir, 'close', true, 'screenSize', [0.5, 0.5]);

for phase = progress(1:numel(phaseNames))
    % Absolute mean epoch
    fig(phaseNames(phase) + " abs mean, Epoch versus value"); clf
    clear G;
    G(1,2) = gramm('x', combined_hpc.epoch, 'y', abs(combined_hpc.val), 'color', C(combined_hpc.patternSym), 'lightness', combined_hpc.component, 'subset', combined_hpc.phase == phase);
    G(1,2).facet_grid(combined_hpc.patternSym, combined_hpc.component);
    G(1,1) = gramm('x', combined_pfc.epoch, 'y', abs(combined_pfc.val), 'color', C(combined_pfc.patternSym), 'lightness', combined_pfc.component, 'subset', combined_pfc.phase == phase);
    G(1,1).facet_grid(combined_pfc.pattern, combined_pfc.component);
    arrayfun(@(g) g.stat_summary('type','fitnormalci','geom',{'bar','black_errorbar'}, 'dodge',0, 'width', 0.9),G(:));
    arrayfun(@(g) g.axe_property('ylim',[0,20]), G(:), 'UniformOutput', false);
    G.set_color_options('lightness_range',[15 80]);
    G.set_names('x','Epoch','y','abs($\mu$)');
    G.set_order_options('color',order);
    G.set_text_options('interpreter', 'latex');
    G.set_title(phaseNames(phase) + ' abs. mean over Epochs')
    G(1).set_title('hpc-pfc')
    G(2).set_title('hpc-hpc')
    G.draw()
end
util.fig.saveByName(savedir, 'close', true, 'screenSize', [0.5, 0.5]);
