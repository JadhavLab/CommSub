function [psid_struct, z, B] = calculate(y, Z, nx, ni, hankel_embedding, varargin)
% y - neural data
% z - behavior

ip = inputParser;
ip.addParameter('beh_fields',[]);
ip.addParameter('set_name',"");
ip.addParameter('train',"leave1out"); % leave1out or neighbor
ip.addParameter('bigdata', false);
ip.parse(varargin{:})
Opt = ip.Results;

beh_fields = Opt.beh_fields;

% Process behavior
epochs = Z.epoch;
z = table2array(Z(:, beh_fields));
norm = @(x) (x - min(x))./(max(x) - min(x));
z = norm(z) - 0.5;
[epochGroups, epochs] = findgroups(epochs);
uepochGroups = unique(epochGroups);

psid_struct = struct('nx',nx,...
    'ni',ni,...
    'hankel_embedding',hankel_embedding);

% CV partition
Z_train = {}; 
Y_train = {};
Z_test  = {}; 
Y_test  = {};
if Opt.train == "leave1out"
    for g = uepochGroups'
        gselect = g == epochGroups;
        B{g} = Z(g == epochGroups, :);
        psid_struct.Z_test{g}  = z(gselect, :);
        psid_struct.Z_train{g} = z(~gselect, :);
        psid_struct.Y_test{g}  = y(gselect, :);
        psid_struct.Y_train{g} = y(~gselect, :);
    end
elseif Opt.train == "neighbor"
    previous = uepochGroups - 1;
    previous(previous == 0) = 2; % edge only gets to see the next, not previous
    next     = uepochGroups + 1
    previous(previous == max(uepochGroups)) = max(uepochGroups)-1; % edge only gets to see the previous, not next
    for g = uepochGroups'
        gselect = g == epochGroups;
        trainselect = ismember(epochGroups, [previous(g), next(g)]);
        B{g} = Z(g == epochGroups, :);
        %if Opt.bigdata
        %    psid_struct.Z_test{g}  = single(z(gselect, :));
        %    psid_struct.Z_train{g} = single(z(trainselect, :));
        %    psid_struct.Y_test{g}  = single(y(gselect, :));
        %    psid_struct.Y_train{g} = single(y(trainselect, :));
        %else
            psid_struct.Z_test{g}  = z(gselect, :);
            psid_struct.Z_train{g} = z(trainselect, :);
            psid_struct.Y_test{g}  = y(gselect, :);
            psid_struct.Y_train{g} = y(trainselect, :);
        %end
    end
end

tic;
psid_struct.Z_pred = {};
psid_struct.Y_pred = {};
psid_struct.X_pred = {};
for g = progress(uepochGroups', 'Title', 'Training psid')
    [psid_struct.idSys{g}, psid_struct.WS{g}] = PSID(psid_struct.Y_train{g}, psid_struct.Z_train{g}, nx, ni, hankel_embedding);
    if Opt.bigdata
        psid_struct.Y_train{g} = [];
        psid_struct.Z_train{g} = [];
    end
end

for g = progress(uepochGroups', 'Title', 'Testing psid')
     [Z_pred, Y_pred, X_pred] = PSIDPredict(psid_struct.idSys{g}, psid_struct.Y_test{g});
     if Opt.bigdata
         Z_pred = single(Z_pred);
         Y_pred = single(Y_pred);
         X_pred = single(X_pred);
     end
     if any(isnan(Z_pred),'all') || any(isnan(Y_pred),'all') || any(isnan(X_pred),'all')
         error("Nan result! You may need to change parameters nx, ni, hankel_embedding");
     end
     [psid_struct.Z_pred{g}, psid_struct.Y_pred{g}, psid_struct.X_pred{g}] = deal(Z_pred, Y_pred, X_pred);
    if any(isnan(Z_pred),'all')
        warning('Nan predictions!')
    end
end

psid_struct.Z_diff = {};
for g = progress(uepochGroups', 'Title', 'Diffing psid')
    accuracy = (psid_struct.Z_test{g} - psid_struct.Z_pred{g})./(psid_struct.Z_test{g});
    psid_struct.Z_diff{g} = nanmean(accuracy,'all');
end

for v = 1:2:numel(varargin)
    psid_struct.(varargin{v}) = varargin{v+1};
end

toc

psid_struct.set_name = Opt.set_name;
