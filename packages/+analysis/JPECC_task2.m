function Jpecc_Results = JPECC_task2(Source, Target, Project, Patterns_overall, Option)
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

    [curr_source, curr_target, is_source] = analysis.Sampling_JPECC_task2(S, T, upper, lower);

    num_dims = 3; 

    rval_single = cell(num_dims, 30);
    pval_single = cell(num_dims, 30);

    rval_combined = cell(num_dims, 30);
    pval_combined = cell(num_dims, 30);

    if Project == 1
        r_theta_single = cell(num_dims, 30);
        r_delta_single = cell(num_dims, 30);
        r_ripple_single = cell(num_dims, 30);

        r_theta_combined = cell(num_dims, 30);
        r_delta_combined = cell(num_dims, 30);
        r_ripple_combined = cell(num_dims, 30);
    end

    if is_source == 1
        for w = 1:30
            tens_curr_source = curr_source{1,w};
            tens_curr_target = curr_target;
            cs = tens_curr_source;
            ct = tens_curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan5Samples = sum(~nan_rows) <= 5;
            if lessThan5Samples
                for r = 1:3
                rval_single{r,w} = nan;
                pval_single{r,w} = nan;
                rval_combined{r,w} = nan;
                pval_combined{r,w} = nan;
                if Project == 1
                    r_theta_single{r,w} = nan;
                    r_delta_single{r,w} = nan;
                    r_ripple_single{r,w} = nan;
                    r_theta_combined{r,w} = nan;
                    r_delta_combined{r,w} = nan;
                    r_ripple_combined{r,w} = nan;
                end
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
            
            % use cross-validation
            U = zeros(size(cs,1),num_dims);
            V = zeros(size(cs,1),num_dims);

            if Project == 1
                U_theta = zeros(size(cs,1),num_dims);
                V_theta = zeros(size(cs,1),num_dims);
                U_delta = zeros(size(cs,1),num_dims);
                V_delta = zeros(size(cs,1),num_dims);
                U_ripple = zeros(size(cs,1),num_dims);
                V_ripple = zeros(size(cs,1),num_dims);
            end

            for dim = 1:num_dims
                for k = 1:cvNumFolds

                    [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));
                    U(cvp.test(k), dim) = cs(cvp.test(k),:)*A(:,dim);
                    V(cvp.test(k), dim) = ct(cvp.test(k),:)*B(:,dim);

                    if Project == 1
                        U_theta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_theta*A(:,dim);
                        V_theta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_theta*B(:,dim);
                        U_delta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_delta*A(:,dim);
                        V_delta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_delta*B(:,dim);
                        U_ripple(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_ripple*A(:,dim);
                        V_ripple(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_ripple*B(:,dim);
                    end
                end
            end

            for dim = 1:num_dims
                [rval_single{dim, w}, pval_single{dim, w}] = corr(U(:,dim), V(:,dim));
                if Project == 1
                    [r_theta_single{dim, w}, ~] = corr(U_theta(:,dim), V_theta(:,dim));
                    [r_delta_single{dim, w}, ~] = corr(U_delta(:,dim), V_delta(:,dim));
                    [r_ripple_single{dim, w}, ~] = corr(U_ripple(:,dim), V_ripple(:,dim));
                end

                U_combined = reshape(U(:,1:dim), [], 1);
                V_combined = reshape(V(:,1:dim), [], 1);
                [rval_combined{dim, w}, pval_combined{dim, w}] = corr(U_combined, V_combined);
                if Project == 1
                    [r_theta_combined{dim, w}, ~] = corr(reshape(U_theta(:,1:dim), [], 1), reshape(V_theta(:,1:dim), [], 1));
                    [r_delta_combined{dim, w}, ~] = corr(reshape(U_delta(:,1:dim), [], 1), reshape(V_delta(:,1:dim), [], 1));
                    [r_ripple_combined{dim, w}, ~] = corr(reshape(U_ripple(:,1:dim), [], 1), reshape(V_ripple(:,1:dim), [], 1));
                end
            end
        end

    elseif is_source ==0
        for w = 1:30
            tens_curr_source = curr_source;
            tens_curr_target = curr_target{1,w};
            cs = tens_curr_source;
            ct = tens_curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan5Samples = sum(~nan_rows) <= 5;
            if lessThan5Samples
                for r = 1:3
                rval_single{r,w} = nan;
                pval_single{r,w} = nan;
                rval_combined{r,w} = nan;
                pval_combined{r,w} = nan;
                if Project == 1
                    r_theta_single{r,w} = nan;
                    r_delta_single{r,w} = nan;
                    r_ripple_single{r,w} = nan;
                    r_theta_combined{r,w} = nan;
                    r_delta_combined{r,w} = nan;
                    r_ripple_combined{r,w} = nan;
                end
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
            
            % use cross-validation
            U = zeros(size(cs,1),num_dims);
            V = zeros(size(cs,1),num_dims);

            if Project == 1
                U_theta = zeros(size(cs,1),num_dims);
                V_theta = zeros(size(cs,1),num_dims);
                U_delta = zeros(size(cs,1),num_dims);
                V_delta = zeros(size(cs,1),num_dims);
                U_ripple = zeros(size(cs,1),num_dims);
                V_ripple = zeros(size(cs,1),num_dims);
            end

            for dim = 1:num_dims
                for k = 1:cvNumFolds

                    [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));
                    U(cvp.test(k), dim) = cs(cvp.test(k),:)*A(:,dim);
                    V(cvp.test(k), dim) = ct(cvp.test(k),:)*B(:,dim);

                    if Project == 1
                        U_theta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_theta*A(:,dim);
                        V_theta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_theta*B(:,dim);
                        U_delta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_delta*A(:,dim);
                        V_delta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_delta*B(:,dim);
                        U_ripple(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_ripple*A(:,dim);
                        V_ripple(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_ripple*B(:,dim);
                    end
                end
            end

            for dim = 1:num_dims
                [rval_single{dim, w}, pval_single{dim, w}] = corr(U(:,dim), V(:,dim));
                if Project == 1
                    [r_theta_single{dim, w}, ~] = corr(U_theta(:,dim), V_theta(:,dim));
                    [r_delta_single{dim, w}, ~] = corr(U_delta(:,dim), V_delta(:,dim));
                    [r_ripple_single{dim, w}, ~] = corr(U_ripple(:,dim), V_ripple(:,dim));
                end

                U_combined = reshape(U(:,1:dim), [], 1);
                V_combined = reshape(V(:,1:dim), [], 1);
                [rval_combined{dim, w}, pval_combined{dim, w}] = corr(U_combined, V_combined);
                if Project == 1
                    [r_theta_combined{dim, w}, ~] = corr(reshape(U_theta(:,1:dim), [], 1), reshape(V_theta(:,1:dim), [], 1));
                    [r_delta_combined{dim, w}, ~] = corr(reshape(U_delta(:,1:dim), [], 1), reshape(V_delta(:,1:dim), [], 1));
                    [r_ripple_combined{dim, w}, ~] = corr(reshape(U_ripple(:,1:dim), [], 1), reshape(V_ripple(:,1:dim), [], 1));
                end
            end
        end

    else
        for w = 1:30
            cs = curr_source;
            ct = curr_target;

            nan_rows = any(isnan(cs), 2) | ...
                      any(isnan(ct), 2);
            cs = cs(~nan_rows,:);
            ct = ct(~nan_rows,:);

            lessThan5Samples = sum(~nan_rows) <= 5;
            if lessThan5Samples
                for r = 1:3
                rval_single{r,w} = nan;
                pval_single{r,w} = nan;
                rval_combined{r,w} = nan;
                pval_combined{r,w} = nan;
                if Project == 1
                    r_theta_single{r,w} = nan;
                    r_delta_single{r,w} = nan;
                    r_ripple_single{r,w} = nan;
                    r_theta_combined{r,w} = nan;
                    r_delta_combined{r,w} = nan;
                    r_ripple_combined{r,w} = nan;
                end
                end
                continue
            end

            N = size(cs,1); %number of trials
            cvp = cvpartition(N, 'KFold', cvNumFolds);
            
            % use cross-validation
            U = zeros(size(cs,1),num_dims);
            V = zeros(size(cs,1),num_dims);

            if Project == 1
                U_theta = zeros(size(cs,1),num_dims);
                V_theta = zeros(size(cs,1),num_dims);
                U_delta = zeros(size(cs,1),num_dims);
                V_delta = zeros(size(cs,1),num_dims);
                U_ripple = zeros(size(cs,1),num_dims);
                V_ripple = zeros(size(cs,1),num_dims);
            end

            for dim = 1:num_dims
                for k = 1:cvNumFolds

                    [A,B] = canoncorr(...
                        cs(cvp.training(k),:),...
                        ct(cvp.training(k),:));
                    U(cvp.test(k), dim) = cs(cvp.test(k),:)*A(:,dim);
                    V(cvp.test(k), dim) = ct(cvp.test(k),:)*B(:,dim);

                    if Project == 1
                        U_theta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_theta*A(:,dim);
                        V_theta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_theta*B(:,dim);
                        U_delta(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_delta*A(:,dim);
                        V_delta(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_delta*B(:,dim);
                        U_ripple(cvp.test(k), dim) = cs(cvp.test(k),:)*P_a_ripple*A(:,dim);
                        V_ripple(cvp.test(k), dim) = ct(cvp.test(k),:)*P_b_ripple*B(:,dim);
                    end
                end
            end

            for dim = 1:num_dims
                [rval_single{dim, w}, pval_single{dim, w}] = corr(U(:,dim), V(:,dim));
                if Project == 1
                    [r_theta_single{dim, w}, ~] = corr(U_theta(:,dim), V_theta(:,dim));
                    [r_delta_single{dim, w}, ~] = corr(U_delta(:,dim), V_delta(:,dim));
                    [r_ripple_single{dim, w}, ~] = corr(U_ripple(:,dim), V_ripple(:,dim));
                end
                
                U_combined = reshape(U(:,1:dim), [], 1);
                V_combined = reshape(V(:,1:dim), [], 1);
                [rval_combined{dim, w}, pval_combined{dim, w}] = corr(U_combined, V_combined);
                if Project == 1
                    [r_theta_combined{dim, w}, ~] = corr(reshape(U_theta(:,1:dim), [], 1), reshape(V_theta(:,1:dim), [], 1));
                    [r_delta_combined{dim, w}, ~] = corr(reshape(U_delta(:,1:dim), [], 1), reshape(V_delta(:,1:dim), [], 1));
                    [r_ripple_combined{dim, w}, ~] = corr(reshape(U_ripple(:,1:dim), [], 1), reshape(V_ripple(:,1:dim), [], 1));
                end
            end
        end
    end

    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).rsingle = rval_single;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).psingle = pval_single;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).rcombined = rval_combined;
    Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).pcombined = pval_combined;
    if Project == 1
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).theta_single = r_theta_single;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).delta_single = r_delta_single;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).ripple_single = r_ripple_single;

        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).theta_combined = r_theta_combined;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).delta_combined = r_delta_combined;
        Jpecc_Results.jpecc((i-1)*30+j, (p-1)*30+q).ripple_combined = r_ripple_combined;
    end
end
end

end
end

