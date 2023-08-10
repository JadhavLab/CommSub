function Jpecc_Results = JPECC_task(Source, Target, Project, Patterns_overall, Option)
% Computes lagged JPECC between patterns in two brain araes

%Patterns(i).X_source = [Neurons, Times, Trials]

% Number of cross validation folds.
cvNumFolds = Option.jpecc.cvnum;
upper = Option.threshold(2);
lower = Option.threshold(1);

if Project == 1
    a_theta = Patterns_overall(2,1).cca.a;
    b_theta = Patterns_overall(2,1).cca.b;
    a_delta = Patterns_overall(2,2).cca.a;
    b_delta = Patterns_overall(2,2).cca.b;
    a_ripple = Patterns_overall(2,3).cca.a;
    b_ripple = Patterns_overall(2,3).cca.b;
       
    P_a_theta = a_theta*inv(a_theta'*a_theta)*a_theta';
    P_b_theta = b_theta*inv(b_theta'*b_theta)*b_theta';
    P_a_delta = a_delta*inv(a_delta'*a_delta)*a_delta';
    P_b_delta = b_delta*inv(b_delta'*b_delta)*b_delta';
    P_a_ripple = a_ripple*inv(a_ripple'*a_ripple)*a_ripple';
    P_b_ripple = b_ripple*inv(b_ripple'*b_ripple)*b_ripple';

end

for i = 1:size(Source,1)
for j = 1:size(Source,2)

    disp("----------")
    disp([i,j])
    disp("----------")

    S = Source{i,j};

for p = 1:size(Target,1)

    disp("----------")
    disp(p)
    disp("----------")

for q = 1:size(Target,2)
    T = Target{p,q};
    
    %{
    tens_curr_source = S.val;
    tens_curr_target = T.val;
    source_time = S.time;
    target_time = T.time;

    % sampling
    num_samples_source = size(tens_curr_source, 1);
    num_samples_target = size(tens_curr_target, 1);

    if num_samples_source > num_samples_target

        % Sort the time vectors
        [source_time, source_idx] = sort(source_time);
        [target_time, ~] = sort(target_time);

        % Initialize the selected source indices
        selected_source_idx = zeros(1, num_samples_target);

        % For each of the first n target times
        for d = 1:num_samples_target
            % Calculate the time differences between the current target time and all source times
            time_diffs = abs(source_time - target_time(d));
    
            % Find the source time that is closest to the current target time
            [~, closest_source_idx] = min(time_diffs);
    
            % Add the index of the closest source time to the selected source indices
            selected_source_idx(d) = source_idx(closest_source_idx);
    
            % Remove the closest source time from consideration in future iterations
            source_time(closest_source_idx) = [];
            source_idx(closest_source_idx) = [];
        end

        tens_curr_source = tens_curr_source(selected_source_idx, :);

    elseif num_samples_source < num_samples_target

        % Sort the time vectors
        [source_time, ~] = sort(source_time);
        [target_time, target_idx] = sort(target_time);

        % Initialize the selected source indices
        selected_target_idx = zeros(1, num_samples_source);

        % For each of the first n target times
        for d = 1:num_samples_source
            % Calculate the time differences between the current target time and all source times
            time_diffs = abs(target_time - source_time(d));
    
            % Find the source time that is closest to the current target time
            [~, closest_target_idx] = min(time_diffs);
    
            % Add the index of the closest source time to the selected source indices
            selected_target_idx(d) = target_idx(closest_target_idx);
    
            % Remove the closest source time from consideration in future iterations
            target_time(closest_target_idx) = [];
            target_idx(closest_target_idx) = [];
        end

        tens_curr_target = tens_curr_target(selected_target_idx, :);

    end
    %}

    [curr_source, curr_target, is_source] = analysis.Sampling_JPECC_task2(S, T, upper, lower);

    rval1 = cell(1,10);
    pval1 = cell(1,10);
    rval2 = cell(1,10);
    pval2 = cell(1,10);
    rval = cell(1,10);
    if Project == 1
        r_theta_1 = cell(1,10);
        r_delta_1 = cell(1,10);
        r_ripple_1 = cell(1,10);
        r_theta_2 = cell(1,10);
        r_delta_2 = cell(1,10);
        r_ripple_2 = cell(1,10);
    end

    if is_source == 1
        for w = 1:10
            tens_curr_source = curr_source{1,w};
            tens_curr_target = curr_target;
            cs = tens_curr_source;
            ct = tens_curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan4Samples = sum(~nan_rows) <= 4;
            if lessThan4Samples
                rval1{1,w} = nan;
                pval1{1,w} = nan;
                rval2{1,w} = nan;
                pval2{1,w} = nan;
                rval{1,w} = nan;
                if Project == 1
                    r_theta_1{1,w} = nan;
                    r_delta_1{1,w} = nan;
                    r_ripple_1{1,w} = nan;
                    r_theta_2{1,w} = nan;
                    r_delta_2{1,w} = nan;
                    r_ripple_2{1,w} = nan;
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
        
            % use cross-validation
            U1 = zeros(size(cs,1),1);
            V1 = zeros(size(cs,1),1);
            U2 = zeros(size(cs,1),1);
            V2 = zeros(size(cs,1),1);

            if Project == 1
                U1_theta = zeros(size(cs,1),1);
                V1_theta = zeros(size(cs,1),1);
                U1_delta = zeros(size(cs,1),1);
                V1_delta = zeros(size(cs,1),1);
                U1_ripple = zeros(size(cs,1),1);
                V1_ripple = zeros(size(cs,1),1);

                U2_theta = zeros(size(cs,1),1);
                V2_theta = zeros(size(cs,1),1);
                U2_delta = zeros(size(cs,1),1);
                V2_delta = zeros(size(cs,1),1);
                U2_ripple = zeros(size(cs,1),1);
                V2_ripple = zeros(size(cs,1),1);
            end

            for k = 1:cvNumFolds
            
                % generate canonical dimension of each matrix on the training set  

                [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));


                % project the test set onto the canonical dimension
                U1(cvp.test(k)) = cs(cvp.test(k),:)*A(:,1);
                V1(cvp.test(k)) = ct(cvp.test(k),:)*B(:,1);

                U2(cvp.test(k)) = cs(cvp.test(k),:)*A(:,2);
                V2(cvp.test(k)) = ct(cvp.test(k),:)*B(:,2);

                if Project == 1
                    U1_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,1);
                    V1_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,1);
                    U1_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,1);
                    V1_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,1);
                    U1_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,1);
                    V1_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,1);

                    U2_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,2);
                    V2_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,2);
                    U2_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,2);
                    V2_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,2);
                    U2_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,2);
                    V2_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,2);

                end
  
            end

            % correlate the projections. since each test set's projections will
            % be zero mean already, we can just combine them all here
            [rval1{1,w}, pval1{1,w}] = corr(U1,V1);
            [rval2{1,w}, pval2{1,w}] = corr(U2,V2);
            if Project == 1
                [r_theta_1{1,w},~] = corr(U1_theta,V1_theta);
                [r_delta_1{1,w},~] = corr(U1_delta,V1_delta);
                [r_ripple_1{1,w},~] = corr(U1_ripple,V1_ripple);

                [r_theta_2{1,w},~] = corr(U2_theta,V2_theta);
                [r_delta_2{1,w},~] = corr(U2_delta,V2_delta);
                [r_ripple_2{1,w},~] = corr(U2_ripple,V2_ripple);
            end
        end

    elseif is_source ==0
        for w = 1:10
            tens_curr_source = curr_source;
            tens_curr_target = curr_target{1,w};
            cs = tens_curr_source;
            ct = tens_curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan4Samples = sum(~nan_rows) <= 4;
            if lessThan4Samples
                rval1{1,w} = nan;
                pval1{1,w} = nan;
                rval2{1,w} = nan;
                pval2{1,w} = nan;
                rval{1,w} = nan;
                if Project == 1
                    r_theta_1{1,w} = nan;
                    r_delta_1{1,w} = nan;
                    r_ripple_1{1,w} = nan;
                    r_theta_2{1,w} = nan;
                    r_delta_2{1,w} = nan;
                    r_ripple_2{1,w} = nan;
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
        
            % use cross-validation
            U1 = zeros(size(cs,1),1);
            V1 = zeros(size(cs,1),1);
            U2 = zeros(size(cs,1),1);
            V2 = zeros(size(cs,1),1);

            if Project == 1
                U1_theta = zeros(size(cs,1),1);
                V1_theta = zeros(size(cs,1),1);
                U1_delta = zeros(size(cs,1),1);
                V1_delta = zeros(size(cs,1),1);
                U1_ripple = zeros(size(cs,1),1);
                V1_ripple = zeros(size(cs,1),1);

                U2_theta = zeros(size(cs,1),1);
                V2_theta = zeros(size(cs,1),1);
                U2_delta = zeros(size(cs,1),1);
                V2_delta = zeros(size(cs,1),1);
                U2_ripple = zeros(size(cs,1),1);
                V2_ripple = zeros(size(cs,1),1);
            end

            for k = 1:cvNumFolds
            
                % generate canonical dimension of each matrix on the training set  

                [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));


                % project the test set onto the canonical dimension
                U1(cvp.test(k)) = cs(cvp.test(k),:)*A(:,1);
                V1(cvp.test(k)) = ct(cvp.test(k),:)*B(:,1);

                U2(cvp.test(k)) = cs(cvp.test(k),:)*A(:,2);
                V2(cvp.test(k)) = ct(cvp.test(k),:)*B(:,2);
                if Project == 1
                    U1_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,1);
                    V1_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,1);
                    U1_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,1);
                    V1_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,1);
                    U1_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,1);
                    V1_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,1);

                    U2_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,2);
                    V2_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,2);
                    U2_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,2);
                    V2_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,2);
                    U2_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,2);
                    V2_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,2);
                end
  
            end

            % correlate the projections. since each test set's projections will
            % be zero mean already, we can just combine them all here
            [rval1{1,w}, pval1{1,w}] = corr(U1,V1);
            [rval2{1,w}, pval2{1,w}] = corr(U2,V2);
            if Project == 1
                [r_theta_1{1,w},~] = corr(U1_theta,V1_theta);
                [r_delta_1{1,w},~] = corr(U1_delta,V1_delta);
                [r_ripple_1{1,w},~] = corr(U1_ripple,V1_ripple);

                [r_theta_2{1,w},~] = corr(U2_theta,V2_theta);
                [r_delta_2{1,w},~] = corr(U2_delta,V2_delta);
                [r_ripple_2{1,w},~] = corr(U2_ripple,V2_ripple);
            end
        end
    else
        for w = 1:10
            cs = curr_source;
            ct = curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan4Samples = sum(~nan_rows) <= 4;
            if lessThan4Samples
                rval1{1,w} = nan;
                pval1{1,w} = nan;
                rval2{1,w} = nan;
                pval2{1,w} = nan;
                rval{1,w} = nan;
                if Project == 1
                    r_theta{1,w} = nan;
                    r_delta{1,w} = nan;
                    r_ripple{1,w} = nan;
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
        
            % use cross-validation
            U1 = zeros(size(cs,1),1);
            V1 = zeros(size(cs,1),1);
            U2 = zeros(size(cs,1),1);
            V2 = zeros(size(cs,1),1);

            if Project == 1
                U1_theta = zeros(size(cs,1),1);
                V1_theta = zeros(size(cs,1),1);
                U1_delta = zeros(size(cs,1),1);
                V1_delta = zeros(size(cs,1),1);
                U1_ripple = zeros(size(cs,1),1);
                V1_ripple = zeros(size(cs,1),1);

                U2_theta = zeros(size(cs,1),1);
                V2_theta = zeros(size(cs,1),1);
                U2_delta = zeros(size(cs,1),1);
                V2_delta = zeros(size(cs,1),1);
                U2_ripple = zeros(size(cs,1),1);
                V2_ripple = zeros(size(cs,1),1);
            end

            for k = 1:cvNumFolds
            
                % generate canonical dimension of each matrix on the training set  

                [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));


                % project the test set onto the canonical dimension
                U1(cvp.test(k)) = cs(cvp.test(k),:)*A(:,1);
                V1(cvp.test(k)) = ct(cvp.test(k),:)*B(:,1);

                U2(cvp.test(k)) = cs(cvp.test(k),:)*A(:,2);
                V2(cvp.test(k)) = ct(cvp.test(k),:)*B(:,2);

                if Project == 1
                    U1_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,1);
                    V1_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,1);
                    U1_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,1);
                    V1_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,1);
                    U1_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,1);
                    V1_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,1);

                    U2_theta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_theta*A(:,2);
                    V2_theta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_theta*B(:,2);
                    U2_delta(cvp.test(k)) = cs(cvp.test(k),:)*P_a_delta*A(:,2);
                    V2_delta(cvp.test(k)) = ct(cvp.test(k),:)*P_b_delta*B(:,2);
                    U2_ripple(cvp.test(k)) = cs(cvp.test(k),:)*P_a_ripple*A(:,2);
                    V2_ripple(cvp.test(k)) = ct(cvp.test(k),:)*P_b_ripple*B(:,2);
                end
            end

            % correlate the projections. since each test set's projections will
            % be zero mean already, we can just combine them all here
            [rval1{1,w}, pval1{1,w}] = corr(U1,V1);
            [rval2{1,w}, pval2{1,w}] = corr(U2,V2);
            if Project == 1
                [r_theta_1{1,w},~] = corr(U1_theta,V1_theta);
                [r_delta_1{1,w},~] = corr(U1_delta,V1_delta);
                [r_ripple_1{1,w},~] = corr(U1_ripple,V1_ripple);

                [r_theta_2{1,w},~] = corr(U2_theta,V2_theta);
                [r_delta_2{1,w},~] = corr(U2_delta,V2_delta);
                [r_ripple_2{1,w},~] = corr(U2_ripple,V2_ripple);
            end
        end        
    end

    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).val1 = rval1;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).p1 = pval1;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).val2 = rval2;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).p2 = pval2;
    if Project == 1
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).theta1 = r_theta_1;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).delta1 = r_delta_1;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).ripple1 = r_ripple_1;

        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).theta2 = r_theta_2;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).delta2 = r_delta_2;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).ripple2 = r_ripple_2;
    end
end
end

end
end

    %{
    cs = tens_curr_source;
    ct = tens_curr_target;

    nan_rows = any(isnan(cs), 2) | ...
               any(isnan(ct), 2);

    cs = cs(~nan_rows,:);
    ct = ct(~nan_rows,:);

    lessThan3Samples = sum(~nan_rows) <= 3;
    if lessThan3Samples
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val1    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p1    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val2    = nan;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p2    = nan;
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
            
    % generate canonical dimension of each matrix on the training set  

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
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val1 = rval1;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p1 = pval1;
        [rval2, pval2] = corr(U2,V2);
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).val2 = rval2;
        Jpecc_Results.jpecc((i-1)*20+j, (p-1)*20+q).p2 = pval2;

end
end

end
end

    %}
