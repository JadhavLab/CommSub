% should use a two-sample t-test, since the co-firing of hpc and pfc are
% independently measured

fig("cofiring per pattern");
h_corrdiff = zeros(1,nPatterns);
p_corrdiff = zeros(1,nPatterns);
for i = 1:3
    subplot(3,1,i);
    hold off;
    cofiring_hh = histogram(withhpc_pairs{i});
    set(cofiring_hh, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    hold on
    cofiring_hp = histogram(withpfc_pairs{i});
    set(cofiring_hp, 'EdgeColor', 'none', 'FaceAlpha', 0.33);
    
    title(patternnames(i))
    ylabel("pairs")
    
    hold on
    avg_hh=line([mean_withhpccorr_pattern(i),mean_withhpccorr_pattern(i)],[0 max(cofiring_hh.Values)]);
    avg_hh.LineStyle = ':'; % Make line dotted
    avg_hh.LineWidth = 2;  % Thicken the line
    avg_hh.Color = 'blue';
    avg.hh.DisplayName = "hpc-hpc corfiring mean";
%     
    hold on
    avg_hp=line([mean_withpfccorr_pattern(i),mean_withpfccorr_pattern(i)],[0 max(cofiring_hh.Values)]);
    avg_hp.LineStyle = ':'; % Make line dotted
    avg_hp.LineWidth = 2;  % Thicken the line
    avg_hp.Color = 'red'; 
    avg.hh.DisplayName = "hpc-pfc corfiring mean";
    
    if source == "hpc"
        legend("hpc-hpc","hpc-pfc")
    else
        legend("pfc-hpc","pfc-pfc")
    end
    
    xlabel("pairwise correlation")
    [h_corrdiff(i),p_corrdiff(i)] = ttest2(double(withhpc_pairs{i}(~isnan(withhpc_pairs{i}))), double(withpfc_pairs{i}(~isnan(withpfc_pairs{i}))));
    xlim([-0.1,0.4])
    
end
        