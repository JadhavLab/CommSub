
%% Get averages of bands
%% ---------------------
% Obtain band averages
filter = sieve.bands(data); % Obtain filtration ranges for all meaningful bands
bandavg = @(x) mean(x,2);   % Function for averaging
clear wpli s1 s2 C
for band = sieve.bandset(); 
    wpli.(band) = bandavg(data.('wpli')(:,filter.(band))); 
    %C.(band) = bandavg(data.('C')(data.subsample_indices,filter.(band))); 
    s1.(band) = bandavg(data.('S1')(:,filter.(band))); 
    s2.(band) = bandavg(data.('S2')(:,filter.(band))); 
end


%                                                   
%. . .,---.|    |    ,---.o          |              
%| | ||---'|    |    |---'.,---.,---.|    ,---.,---.
%| | ||    |    |    |  \ ||   ||   ||    |---'`---.
%`-'-'`    `---'`    `   ``|---'|---'`---'`---'`---'
%                          |    |                   
% ------------


% WPLI_ripple related to HPC and PFC, but most to PFC
% ----------------------------------------------------
bins = 8;
subset = s1.ripple < quantile(s1.ripple, 0.99) & wpli.ripple < quantile(s1.ripple, 0.99);
figure;
clear g
g(1,1) = gramm('x', wpli.ripple, 'y', s1.ripple, 'subset', [])
g(1,1).set_names('x',{'WPLI','150-200hz'}, 'y', {'S_{hpc}','150-200hz'})
g.stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,1) = gramm('x', s1.ripple, 'y',wpli.ripple)
g(2,1).set_names('x',{'S_{hpc}','150-200hz'}, 'y', {'WPLI','150-200hz'})
g(2,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2) = gramm('x', s2.ripple, 'y',wpli.ripple)
g(1,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2).set_names('x',{'S_{pfc}','150-200hz'}, 'y', {'WPLI','150-200hz'})
g(2,2) = gramm('x', wpli.ripple, 'y',s2.ripple)
g(2,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,2).set_names('x',{'WPLI','150-200hz'}, 'y', {'S_{pfc}','150-200hz'})
g(1,3) = gramm('x', s2.ripple, 'y',s1.ripple)
g(1,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,3).set_names('x',{'S_{pfc}','150-200hz'}, 'y', {'S_{hpc}','150-200hz'})
g(2,3) = gramm('x', s1.ripple, 'y',s2.ripple)
g(2,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,3).set_names('x',{'S_{hpc}','150-200hz'}, 'y', {'S_{pfc}','150-200hz'})
for x = 1:size(g,1)
    for y = 1:size(g,2)
        g(x,y).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
        g(x,y).set_text_options('interpreter','tex')
    end
end
g.draw()

% Theta coherence
% --------
bins = 8;
subset = s1.theta < quantile(s1.theta, 0.99) & wpli.theta < quantile(s1.theta, 0.99);
figure;
clear g
g(1,1) = gramm('x', wpli.theta, 'y', s1.theta, 'subset', [])
g(1,1).set_names('x',{'WPLI','6-12hz'}, 'y', {'S_{hpc}','6-12hz'})
g.stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,1) = gramm('x', s1.theta, 'y',wpli.theta)
g(2,1).set_names('x',{'S_{hpc}','6-12hz'}, 'y', {'WPLI','6-12hz'})
g(2,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2) = gramm('x', s2.theta, 'y',wpli.theta)
g(1,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2).set_names('x',{'S_{pfc}','6-12hz'}, 'y', {'WPLI','6-12hz'})
g(2,2) = gramm('x', wpli.theta, 'y',s2.theta)
g(2,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,2).set_names('x',{'WPLI','6-12hz'}, 'y', {'S_{pfc}','6-12hz'})
g(1,3) = gramm('x', s2.theta, 'y',s1.theta)
g(1,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,3).set_names('x',{'S_{pfc}','6-12hz'}, 'y', {'S_{hpc}','6-12hz'})
g(2,3) = gramm('x', s1.theta, 'y',s2.theta)
g(2,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,3).set_names('x',{'S_{hpc}','6-12hz'}, 'y', {'S_{pfc}','6-12hz'})
for x = 1:size(g,1)
    for y = 1:size(g,2)
        g(x,y).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
        g(x,y).set_text_options('interpreter','tex')
    end
end
g.draw()

% Delta
% --------
bins = 8;
subset = s1.delta < quantile(s1.delta, 0.99) & wpli.delta < quantile(s1.delta, 0.99);
figure;
clear g
g(1,1) = gramm('x', wpli.delta, 'y', s1.delta, 'subset', [])
g(1,1).set_names('x',{'WPLI','0.5-4hz'}, 'y', {'S_{hpc}','0.5-4hz'})
g.stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,1) = gramm('x', s1.delta, 'y',wpli.delta)
g(2,1).set_names('x',{'S_{hpc}','0.5-4hz'}, 'y', {'WPLI','0.5-4hz'})
g(2,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2) = gramm('x', s2.delta, 'y',wpli.delta)
g(1,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,2).set_names('x',{'S_{pfc}','0.5-4hz'}, 'y', {'WPLI','0.5-4hz'})
g(2,2) = gramm('x', wpli.delta, 'y',s2.delta)
g(2,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,2).set_names('x',{'WPLI','0.5-4hz'}, 'y', {'S_{pfc}','0.5-4hz'})
g(1,3) = gramm('x', s2.delta, 'y',s1.delta)
g(1,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,3).set_names('x',{'S_{pfc}','0.5-4hz'}, 'y', {'S_{hpc}','0.5-4hz'})
g(2,3) = gramm('x', s1.delta, 'y',s2.delta)
g(2,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,3).set_names('x',{'S_{hpc}','0.5-4hz'}, 'y', {'S_{pfc}','0.5-4hz'})
for x = 1:size(g,1)
    for y = 1:size(g,2)
        g(x,y).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
        g(x,y).set_text_options('interpreter','tex')
    end
end
g.draw()

% PROBLEM: No rewarded outbound or unrewarded inbound ... check this! If true, then wtf.
figure;
clear g
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
coordinate = [0 0];
rewList = ["t_{-2}","t_{-1}","t_{0}","t_{1}","t_{2}"];
bandList = ["delta", "theta", "ripple"];
oRew = ["Null","Error","Correct"];
tmp  = double(data.orig_reward)+2;
tmp(isnan(tmp)) = 1;
tmp = int32(round(tmp));
oRew = reshape(oRew(tmp(:)),size(tmp));
for band = bandList
    coordinate(1,1) = coordinate(1,1) + 1;
    coordinate(1,2) = 0;
for reward = rewList
    coordinate(1,2) = coordinate(1,2) + 1;
    C=num2cell(coordinate);
    %C
    g(C{:}) = gramm('x',data.orig_trajdist(:), 'y', wpli.(band)(:), 'lightness', categorical(oRew(:,coordinate(1,2))));
    g(C{:}).geom_hline('yintercept',0)
    g(C{:}).geom_vline('xintercept',0)
    g(C{:}).set_stat_options('alpha',0.001);
    g(C{:}).set_color_options('lightness_range', [0,95], 'chroma_range', [20, 80], 'hue_range', [0 + (coordinate(1)-1)*360/numel(bandList), 0 + coordinate(1)*360/numel(bandList)])
    g(C{:}).stat_summary('bin_in',50,'interp_in',100,'interp','linear','geom',{'solid_area'}, 'type','ci','setylim',true)
    g(C{:}).set_names('x','Trajdist', 'y', {'WPLI', sieve.bandstring(band)}, 'color', 'Correct/Error')
    g(C{:}).set_order_options('lightness', ["Error","Correct"]);
    g(C{:}).set_text_options('interpreter','tex')
    %g(C{:}).set_color_options('interpreter','tex')
    g(C{:}).no_legend()
end
end
%g().axe_property('ylim',[-0.025 0.05]);
g.draw()

% PROBLEM: No rewarded outbound or unrewarded inbound ... check this! If true, then wtf.
figure;
clear g
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
coordinate = [0 0];
rewList = ["t_{-2}","t_{-1}","t_{0}","t_{1}","t_{2}"];
bandList = ["beta", "ripplegamma", "lowgamma", "highgamma", "epsilon"];
oRew = ["Null","Error","Correct"];
tmp  = double(data.orig_reward)+2;
tmp(isnan(tmp)) = 1;
tmp = int32(round(tmp));
oRew = reshape(oRew(tmp(:)),size(tmp));
for band = bandList
    coordinate(1,1) = coordinate(1,1) + 1;
    coordinate(1,2) = 0;
for reward = rewList
    coordinate(1,2) = coordinate(1,2) + 1;
    C=num2cell(coordinate);
    %C
    g(C{:}) = gramm('x',data.orig_trajdist(:), 'y', wpli.(band)(:), 'lightness', categorical(oRew(:,coordinate(1,2))));
    g(C{:}).geom_hline('yintercept',0)
    g(C{:}).geom_vline('xintercept',0)
    g(C{:}).set_stat_options('alpha',0.001);
    g(C{:}).set_color_options('lightness_range', [0,95], 'chroma_range', [20, 80], 'hue_range', [0 + (coordinate(1)-1)*360/numel(bandList), 0 + coordinate(1)*360/numel(bandList)])
    g(C{:}).stat_summary('bin_in',50,'interp_in',100,'interp','linear','geom',{'solid_area'}, 'type','ci','setylim',true)
    g(C{:}).set_names('x','Trajdist', 'y', {band, 'WPLI', sieve.bandstring(band)}, 'color', 'Correct/Error')
    g(C{:}).set_order_options('lightness', ["Error","Correct"]);
    g(C{:}).set_text_options('interpreter','tex')
    %g(C{:}).set_color_options('interpreter','tex')
    g(C{:}).no_legend()
end
end
%g().axe_property('ylim',[-0.025 0.05]);
g.draw()

% Ripple comparison, but first have to actually HAVE samples!
figure;
field = 'ripple'
I = ["inbound","outbound"];
I = I(double(data.orig_trajdist>0) + 1);
g = gramm('x',categorical(I), 'y', wpli.(field), 'color', data.reward(:,4))
g.set_stat_options('alpha',0.001);
g.geom_hline('yintercept',0)
g.stat_summary('bin_in',2,'geom',{'bar', 'black_errorbar'}, 'type','ci');
g.draw()
outbound= data.orig_trajdist<=0;
correct = data.reward(:,4) > 0;
XT = xticks();
[p, obsDiff, effectSize] = permutationTest(wpli.(field)(outbound & correct), wpli.(field)(outbound & ~correct), 1000, 'plotresult',0);
sigstar(num2cell(XT(2:3)), p)
[p, obsDiff, effectSize] = permutationTest(wpli.(field)(~outbound & correct), wpli.(field)(~outbound & ~correct), 1000, 'plotresult', 0);
sigstar(num2cell(XT(4:5)), p)



%                                                                      
%,---.     |             |                   |    |              o     
%|---.,---.|--- ,---.    |---.,   .,---.,---.|--- |---.,---.,---..,---.
%|   ||---'|    ,---|    |   ||   ||   ||   ||    |   ||---'`---.|`---.
%`---'`---'`---'`---^    `   '`---||---'`---'`---'`   '`---'`---'``---'
%                             `---'|                                   
% --------
% Show theta better predicts prefrontal gamma than beta
clear g
zd = zscore(s1.delta, 1, 1);
zt = zscore(s1.theta, 1, 1);
zb = zscore(s1.beta,  1, 1);
zg1 =zscore(s1.gamma, 1, 1);
zg2 =zscore(s2.gamma, 1, 1);
bins = 5;
figure;
D = [];
% Beta/Theta, (Value of theta or beta), value of gamma_hpc, value of gamma_pfc
D = [repmat(0, size(s1.beta)), zb, zg1, zg2; ...
    repmat(1, size(s1.beta)), zt, zg1, zg2];
lab = ["Beta","Theta"];
lab = lab(double(D(:,1))+ones(size(D,1),1))';
g = gramm('x',D(:,2), 'y', D(:,end-1), 'color', categorical(lab), 'subset', D(:,2) < 6);
g.stat_summary('bin_in',40, 'interp_in',60)
g.draw()
clear g
zd = zscore(s1.delta, 1, 1);
zt = zscore(s1.theta, 1, 1);
zb = zscore(s1.beta,  1, 1);
zg1 =zscore(s1.gamma, 1, 1);
zg2 =zscore(s2.gamma, 1, 1);
bins = 5;
figure;
D = [];
% Beta/Theta, (Value of theta or beta), value of gamma_hpc, value of gamma_pfc
D = [repmat(0, size(s1.beta)), zb, zg1, zg2; ...
    repmat(1, size(s1.beta)), zt, zg1, zg2];
lab = ["Beta","Theta"];
lab = lab(double(D(:,1))+ones(size(D,1),1))';
g = gramm('x',D(:,2), 'y', D(:,end), 'color', categorical(lab), 'subset', D(:,2) < 6);
g.stat_summary('bin_in',40, 'interp_in',60)
g.draw()
clear g

clear g
zd = zscore(s1.delta,1);
zt = zscore(s1.theta,1);
zb = zscore(s1.beta,1);
zg1 =zscore(s1.gamma,1); 
zg2 =zscore(s2.gamma,1); 
bins = 5;
figure;
D = [];
for cshift = 0
     D = [D; [repmat(0, size(s1.beta)), zb+zb, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zb, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
end

% ZB > 0
g(1,1) = gramm('x',zb, 'y', D(:,3), 'color', D(:,1), 'subset', zb>0)
g(1,1).facet_grid([],D(:,4), 'scale', 'free_y')
g(1,1).set_names('row','', 'column', '')
g(1,1).geom_hline('yintercept',0)
g(1,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,1).set_names('row','', 'column', '', 'x', {'HPC \beta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(1,1).set_text_options('interpreter','tex')
%D = [];
%for cshift = 0
%     D = [D; [repmat(0, size(s1.beta)), zt+zt, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zt, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
%end
g(1,2) = gramm('x',zt, 'y', D(:,3), 'color', D(:,1), 'subset', zb>0)
g(1,2).facet_grid([],D(:,4), 'scale', 'free_y')
g(1,2).geom_hline('yintercept',0)
g(1,2).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(1,2).set_text_options('interpreter','tex')
g(1,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
% ZB < 0
g(2,1) = gramm('x',zb, 'y', D(:,3), 'color', D(:,1), 'subset', zb<0)
g(2,1).facet_grid([],D(:,4), 'scale', 'free_y')
g(2,1).set_names('row','', 'column', '')
g(2,1).geom_hline('yintercept',0)
g(2,1).set_names('row','', 'column', '', 'x', {'HPC \beta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(2,1).set_text_options('interpreter','tex')
g(2,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
%D = [];
%for cshift = 0
%     D = [D; [repmat(0, size(s1.beta)), zt+zt, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zt, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
%end
g(2,2) = gramm('x',zt, 'y', D(:,3), 'color', D(:,1), 'subset', zb<0)
g(2,2).facet_grid([], D(:,4), 'scale', 'free_y')
g(2,2).geom_hline('yintercept',0)
g(2,2).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(2,2).set_text_options('interpreter','tex')
g(2,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,1).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
g(2,2).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
g.set_title('')
g.draw()

%zd = zscore(mean(data.S1,2),1);
zt = zscore(s1.theta,1);
zb = zscore(s1.beta,1);
zg1 =zscore(s1.gamma,1); 
zg2 =zscore(s2.gamma,1); 
figure;
D = [];
for cshift = 0
     D = [D; [repmat(0, size(s1.beta)), zb+zb, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zb, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
end
% ZB > 0
g(1,1) = gramm('x',zb, 'y', D(:,3), 'color', D(:,1), 'subset', zb>0)
g(1,1).facet_grid([],D(:,4), 'scale', 'free_y')
g(1,1).set_names('row','', 'column', '')
g(1,1).geom_hline('yintercept',0)
g(1,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,1).set_names('row','', 'column', '', 'x', {'HPC \beta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(1,1).set_text_options('interpreter','tex')
%D = [];
%for cshift = 0
%     D = [D; [repmat(0, size(s1.beta)), zt+zt, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zt, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
%end
g(1,2) = gramm('x',zt, 'y', D(:,3), 'color', D(:,1), 'subset', zb>0)
g(1,2).facet_grid([],D(:,4), 'scale', 'free_y')
g(1,2).geom_hline('yintercept',0)
g(1,2).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(1,2).set_text_options('interpreter','tex')
g(1,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(1,3) = gramm('x',zd, 'y', D(:,3), 'color', D(:,1), 'subset', zb>-0.25 & zb<0.25)
g(1,3).facet_grid([],D(:,4), 'scale', 'free_y')
g(1,3).geom_hline('yintercept',0)
g(1,3).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(1,3).set_text_options('interpreter','tex')
g(1,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
% ZB < 0
g(2,1) = gramm('x',zb, 'y', D(:,3), 'color', D(:,1), 'subset', zb<0)
g(2,1).facet_grid([],D(:,4), 'scale', 'free_y')
g(2,1).set_names('row','', 'column', '')
g(2,1).geom_hline('yintercept',0)
g(2,1).set_names('row','', 'column', '', 'x', {'HPC \beta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(2,1).set_text_options('interpreter','tex')
g(2,1).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
%D = [];
%for cshift = 0
%     D = [D; [repmat(0, size(s1.beta)), zt+zt, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zt, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
%end
g(2,2) = gramm('x',zt, 'y', D(:,3), 'color', D(:,1), 'subset', zb<0)
g(2,2).facet_grid([], D(:,4), 'scale', 'free_y')
g(2,2).geom_hline('yintercept',0)
g(2,2).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(2,2).set_text_options('interpreter','tex')
g(2,2).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,3) = gramm('x',zd, 'y', D(:,3), 'color', D(:,1), 'subset', zb>-0.25 & zb<0.25)
g(2,3).facet_grid([],D(:,4), 'scale', 'free_y')
g(2,3).geom_hline('yintercept',0)
g(2,3).set_names('row','', 'column', '', 'x', {'HPC \theta','Z-score'}, 'y', {'\gamma_{hpc} - \gamma_{pfc}','Z-score'})
g(2,3).set_text_options('interpreter','tex')
g(2,3).stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g(2,1).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
g(2,2).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
g(2,3).set_color_options('lightness_range',[0 95],'chroma_range',[0 0], 'chroma', 0);
g.draw()


% --------
zt = zscore(s1.theta,1);
zb = zscore(s1.beta,1);
zg1 =zscore(s1.gamma,1); 
zg2 =zscore(s2.gamma,1); 
bins = 8;
figure;
D = [];
for cshift = -10:10
     D = [D; [repmat(0, size(s1.beta)), zb+zb, zg1, repmat(cshift, size(s1.beta))+repmat(cshift, size(s1.beta))]-[repmat(1, size(s1.beta)), zb, circshift(zg2,cshift), repmat(cshift, size(s1.beta)),]];
end
g = gramm('x',D(:,2), 'y', D(:,3), 'color', D(:,1))
g.facet_grid([],D(:,4), 'scale', 'free_y')
g.set_names('row','', 'column', '')
g.stat_summary('bin_in',bins,'geom',{'area'}, 'type','ci')
g.draw()

% NO DIFFERENCES IN THESE MEAURES
%ratio = zt./zb;
[~,ZT] = sort(zt);
[~,ZB] = sort(zb);
[~,G2] = sort(zg2);
G2 = G2/numel(G2);
[~,G1] = sort(zg2);
G1 = G1/numel(G1);
ZT = ZT/numel(ZT);
ZB = ZB/numel(ZB);
RATIO = ZT./ZB;
DIFFERENCE = ZT-ZB;
PRODDIFFERENCE = G2.*(ZT./ZB);

