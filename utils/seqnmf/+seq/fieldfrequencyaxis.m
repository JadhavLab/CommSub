function frequencyaxis(data, n)

if nargin == 1
    n = 5;
end

%Find frequencies to label
freq = data.f(1,:);
spacing = median(diff(freq));
freq = linspace(freq(1)-spacing/2, freq(end)+spacing/2, n);
freq = string(freq);

% Find fields to label
if isfield(data,'fields')
    fields = data.fields;
else
    warning('File did not contain a fields set. Creating one.')
    fields = {};
    maxnum = max(structfun(@numel, data));
    for field = fieldnames(data)'
        if numel(data.(field{1})) == maxnum
            fields{end+1} = field;
        end
    end
end
fields = string(fields);

% Populate enough ticks for labels
yrange = ylim();
yrange = linspace(min(yrange),max(yrange), numel(freq)*numel(fields));
yticks(yrange); 

% Create a matrix of labels for the raveled two dimensions
ytlabels = (fields' + " " + freq)';
ytlabels = ytlabels(end:-1:1); % Unravel into 1d
yticklabels(ytlabels);  % Set labels
