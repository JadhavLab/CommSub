function J_Patterns = JPECC(Patterns, Option)
% Computes lagged JPECC between patterns in two brain araes

%Patterns(i).X_source = [Neurons, Times, Trials]

% Number of cross validation folds.
cvNumFolds = Option.jpecc.cvnum;

for i = 1:40

    disp("----------")
    disp(i)
    disp("----------")

    p = Patterns(i,2,4);

    a = size(p.X_source,1);
    b = size(p.X_target,1);

    tens_curr_source = reshape(p.X_source,[a,60,2000]);
    tens_curr_target = reshape(p.X_target,[b,60,2000]);

    nBins = size(tens_curr_source,2); 

    for t1 = 1:nBins
        disp(t1)
        for t2 = 1:nBins
    
            cs = squeeze(tens_curr_source(:,t1,:))';
            ct = squeeze(tens_curr_target(:,t2,:))';

            nan_rows = any(isnan(cs), 2) | ...
                    any(isnan(ct), 2);

            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan3Samples = sum(~nan_rows) <= 3;
            if lessThan3Samples
                p.jpecc(t1, t2).val1    = nan;
                p.jpecc(t1, t2).p1    = nan;
                p.jpecc(t1, t2).val2    = nan;
                p.jpecc(t1, t2).p2    = nan;
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
        
            % use cross-validation
            U1 = zeros(size(cs,1),1);
            V1 = zeros(size(cs,1),1);
            U2 = zeros(size(cs,1),1);
            V2 = zeros(size(cs,1),1);

            for k = 1:cvNumFolds
            
                % generate canonical dimension of each matrix on the training
                % set  

                [A,B] = canoncorr(...
                    cs(cvp.training(k),:),...
                    ct(cvp.training(k),:));


                % project the test set onto the canonical dimension
                U1(cvp.test(k)) = cs(cvp.test(k),:)*A(:,1);
                V1(cvp.test(k)) = ct(cvp.test(k),:)*B(:,1);

                U2(cvp.test(k)) = cs(cvp.test(k),:)*A(:,2);
                V2(cvp.test(k)) = ct(cvp.test(k),:)*B(:,2);
  
            end

            % correlate the projections. since each test set's projections will
            % be zero mean already, we can just combine them all here
            [rval1, pval1] = corr(U1,V1);
            p.jpecc(t1,t2).val1 = rval1;
            p.jpecc(t1,t2).p1 = pval1;
            [rval2, pval2] = corr(U2,V2);
            p.jpecc(t1,t2).val2 = rval2;
            p.jpecc(t1,t2).p2 = pval2;

        end
    end

    J_Patterns(i) = p;

end


