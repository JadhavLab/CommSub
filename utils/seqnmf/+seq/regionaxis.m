function regionaxis(data, field, n)

    if nargin == 2
        n = 5;
    end

    assert(ismember(field,{'W_aglindist', 'W_lindist','W_trajdist','aglindist', 'lindist','trajdist'}));
try

    space = data.(field);
    switch field
    case {'aglindist','lindist'}
        yts = linspace(1, size(space,2), n);
        yticks(yts)
        ytlabels = string(1*((1:n)/n))';
        ytlabels = ytlabels(:);
    case 'trajdist'
        yts = linspace(1, size(space,2), n);
        yticks(yts)
        n = n/2;
        ytlabels = ["Out"] + string(round(10*((n:-1:1)/n)));
        ytlabels = [ytlabels; ["In"] + string(round(10*((1:n)/n)))];
        ytlabels = ytlabels';
        ytlabels = ytlabels(:);
    case {'W_aglindist','W_lindist'}
        yts = linspace(1, size(space,1), n);
        yticks(yts)
        ytlabels = string(1*((1:n)/n))';
        ytlabels = ytlabels(:);
    case 'W_trajdist'
        yts = linspace(1, size(space,1), n);
        yticks(yts)
        n = n/2;
        ytlabels = ["Out"] + string(round(10*((n:-1:1)/n)));
        ytlabels = [ytlabels; ["In"] + string(round(10*((1:n)/n)))];
        ytlabels = ytlabels';
        ytlabels = ytlabels(:);
    end
    yticklabels(ytlabels)

    %ytlabels = cellfun(@(x) sprintf('%2.1f', x), ytlabels, 'UniformOutput', false);
    %ytlabels = string([ytlabels{:}]);

catch ME
    warning('Failed to regionaxis')
end
