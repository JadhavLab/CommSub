function D = seqnmf_packfields(data, fields, varargin)
% SEQNMF_PACKFIELDS packs or concatonates several fields into a matrix and post-processes them for seqNMF

% Parse optional inputs
ip = inputParser;
ip.addParameter('removenan', true);
ip.addParameter('zscore', false);
ip.addParameter('groups', []);
ip.addParameter('minmax', false);
ip.addParameter('floorceil', [-1, 2], @isnumeric);
ip.addParameter('emphasis',false)
ip.parse(varargin{:});
opt = ip.Results;

% If no group assignment, all indices belong to the same gruop
if isempty(opt.groups)
    opt.groups = ones(size(data.t));
end

% Figure out how to preallocate
Dsize = cellfun(@(x) size(data.(x)), fields, 'UniformOutput', false);
Dsize = cat(1,Dsize{:});
assert(all(Dsize(:,1) == Dsize(1,1)));
time = Dsize(1,1);
freq_chunks = Dsize(:,2);
freq_chunks = [0; cumsum(freq_chunks)];
Dsize = [time, freq_chunks(end)];

% Concatonate the data
D = size(Dsize);
for group = unique(opt.groups)
    filter = opt.groups == group;
    counter = 1;
for field = fields

    if isequal(field{1},'wpli') && ~isreal(data.(field{1}))
        data.wpli = abs(data.wpli);
    end

    %% --------Construct D and field specific rescalings--------
    if opt.removenan
        frequencyWiseMeans = nanmean(data.(field{1}),1);
        frequencyWiseMeans(isnan(frequencyWiseMeans)) = 0;
        for frequency = 1:size(data.(field{1}),2)
            onefreq = data.(field{1})(:,frequency);
            data.(field{1})(isnan(onefreq),frequency) = frequencyWiseMeans(frequency);
        end
        try
            assert( ~any(isnan(data.(field{1})),'all') )
        catch ME
            keyboard
        end
    end
    
    % Specific zscore
    if iscell(opt.zscore) && ismember(field{1}, opt.zscore)
        D(filter, freq_chunks(counter)+1:freq_chunks(counter+1)) = ...
            zscore(data.(field{1})(filter,:), 1, 1);
    % Specific minmax rescale
    elseif iscell(opt.minmax) && ismember(field{1}, opt.minmax)
        [minmax_] = minmax(double(data.(field{1})));
        D(filter, freq_chunks(counter)+1:freq_chunks(counter+1)) = ...
            (data.(field{1})(filter,:) - minmax_(1)) / diff(minmax_);
        if opt.floorceil
            potential_range = diff(opt.floorceil);
            D(filter, freq_chunks(counter)+1:freq_chunks(counter+1)) = ...
                D(filter, freq_chunks(counter)+1:freq_chunks(counter+1)) * potential_range ...
                + opt.floorceil(1);
        end
    % Else raw data
    else
        D(filter, freq_chunks(counter)+1:freq_chunks(counter+1)) = data.(field{1})(filter,:);
    end
    counter = counter + 1;
end
end

%% Global rescalings
% Global zscore
if ~iscell(opt.zscore) && opt.zscore
    for group = unique(opt.groups)
        filter = opt.groups == group;
        D(filter,:) = zscore(D(filter,:), 1, 1);
    end
end

%Global minmax
if ~iscell(opt.minmax) && opt.minmax
    for group = unique(opt.groups)
        filter = opt.groups == group;
        [min_,max_] = minmax(D(filter,:));
        D(filter,:) = (D(filter,:)-min_) / (max_ - min_);
    end
end

%Global floorceil
if ~isempty(opt.floorceil)
    D(D < opt.floorceil(1)) = opt.floorceil(1);
    D(D > opt.floorceil(2)) = opt.floorceil(2);
    D = D - opt.floorceil(1);
end

% Remove nan
D(isnan(D)) = 0;

