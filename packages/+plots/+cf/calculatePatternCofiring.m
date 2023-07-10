function cofiring = calculatePatternCofiring(pattern1, pattern2)
% return a linear representation of the co-firing stats between pattern 1 and
% pattern 2
% literally grabs the interaction correlations between two fr vectors/matrices
% nPattern actually is the number of cells in the pattern
    
    [~,nPattern1] = size(pattern1);
    [~,nPattern2] = size(pattern2);
    
    all_cells = [pattern1,pattern2];
    corr_matrix = corrcoef(all_cells);
    t_s = corr_matrix(1:nPattern1,nPattern1+1:nPattern2+nPattern1);
    s_t = corr_matrix(nPattern1+1:nPattern2+nPattern1, 1:nPattern1);
    if munge.detectZscore(pattern1)
        disp("z score detected, applying fischer's z transform");
        % apply fischer's z transform
        t_s = atanh(t_s);
        s_t = atanh(s_t);
    end
    cofiring = [t_s(:)',s_t(:)'];
    
end
