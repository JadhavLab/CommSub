function [beh, fnames] = annotate(psid_analysis, beh, varargin)
% Annotates behavior table with psid results

ip = inputParser;
ip.addParameter('as', 'columns');   % how to add the values from psid
ip.addParameter('foi', ["X_pred"]); % fields of interest
ip.addParameter('tag', "");         % what to call this set
ip.parse(varargin{:});
Opt = ip.Results;

Opt.foi = string(Opt.foi);
Opt.tag = string(Opt.tag);
Opt.as  = string(Opt.as);

epochs = findgroups(beh.epoch);
uEpochs = unique(epochs);

if Opt.as == "columns"

    fname = @(tag,field,col) tag + "-" + field + "-" + col;
    fnames = string([]);
    for field = Opt.foi
        columns = 1:size(psid_analysis.(field){1},2);
        for col = columns
            F = fname(Opt.tag, field, col);
            beh.(F) = nan(height(beh), 1);
            fnames = [fnames, F];
        end
        for g = uEpochs'
            for col = columns
                F = fname(Opt.tag, field, col);
                beh.(F)(epochs == g) =  psid_analysis.X_pred{g}(:,col);
            end
        end
    end

elseif Opt.as == "rows"

    fnames = Opt.foi;
    B = [];
    for field = Opt.foi

        tmp = beh;
        tmp.field = repmat(field, [height(beh),1]);

        columns = 1:size(psid_analysis.(field){1},2);
        for col = columns
            for g = uEpochs'
                tmp.col(g==epochs) = col*ones(sum(g==epochs),1);
                tmp.val(g==epochs) = psid_analysis.(field){g}(:,col);
            end
            B = [B; tmp];
        end

    end

    beh = B;

else
    error("Not a recognized 'as' algorithm");
end
