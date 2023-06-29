
%% Get averages of bands
%% ---------------------
% Obtain band averages
filter = sieve.bands(data); % Obtain filtration ranges for all meaningful bands
bandavg = @(x) mean(x,2);   % Function for averaging
for band = sieve.bandset(); 
    wpli.(band) = bandavg(data.('wpli')(:,filter.(band))); 
    s1.(band) = bandavg(data.('S1')(:,filter.(band))); 
    s2.(band) = bandavg(data.('S2')(:,filter.(band))); 
end


% --------
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

% --------
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


% PROBLEM: No rewarded outbound or unrewarded inbound ... check this! If true, then wtf.

figure;
g = gramm('x',data.orig_trajdist, 'y', wpli.ripple, 'color', data.reward(:,4))
g.geom_hline('yintercept',0)
g.stat_summary('bin_in',20,'geom',{'area'}, 'type','fitnormalci')
g.draw()

% Ripple comparison, but first have to actually HAVE samples!
figure;
g = gramm('x',data.orig_trajdist>0, 'y', wpli.ripple, 'color', data.reward(:,4))
g.geom_hline('yintercept',0)
g.stat_summary('bin_in',2,'geom',{'bar', 'black_errorbar'}, 'type','fitnormalci')
g.draw()
[p, obsDiff, effectSize] = permutationTest(wpli.ripple(outbound & correct), wpli.ripple(outbound & ~correct), 1000);
sigstar({[0.15, 0.35]}, p)
[p, obsDiff, effectSize] = permutationTest(wpli.ripple(~outbound & correct), wpli.ripple(~outbound & ~correct), 1000, 'exact', 1, 'plotresult', 1);
sigstar({[0.65, 0.85]}, p)

outbound= data.orig_trajdist<=0;
correct = data.reward(:,4) > 0;

