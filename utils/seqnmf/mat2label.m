function label = mat2label(mat, original_label, iso_hotencoded, repelem_hotencoded)
% LABEL2MAT one hot encoder back to labeled data
%
% matrix : hotencoded or pdf of hotencode to be translated back into a matrix  of labels (for a raw hot encoding matrix) or converted into a table where the columns are labeled (for the pdf)
% original_label : the original full time series of data
% iso_hotencoded : the original hot encoding of that full data series

ip = inputParser;
ip.addParameter('repelem_hotencoded', 0);
ip.addParameter('roundPrecisions', []);
ip.parse(varargin{:});
opt = ip.Results;


if isnumeric(original_label)
    summary = @mean;
elseif isstring(original_label) || ischar(original_label)
    summary = @mode;
    original_label = categorical(original_label)
elseif iscategorical(original_label)
    summary = @mode;
end

if nargin >= 3
    keyboard
    if opt.repelem_hotencoded
        iso_hotencoded = iso_hotencoded(:, 1:repelem_hotencoded:end);
    end
    [label_hotencoded, ~] = find(iso_hotencoded');
    uLabel_hotencoded = unique(label_hotencoded);
    rosettaStone = zeros(size(uLabel_hotencoded));
    for lab = uLabel_hotencoded'
        rosettaStone(lab) = summary(original_label(lab == label_hotencoded));
    end
    if opt.roundPrecision
        rosettaStone = round(rosettaStone, opt.roundPrecision);
    end
end

% For standard hot enconding of 1s and 0s
if all(ismember(mat(:),[1 0]))
    [~,label] = find(mat);
    if nargin >= 3
        label = rosettaStone(label);
    end
% For probability distrubtion of hot encoding
% ... in this case, we just convert the mat to
% a table and label the columns
else
    label = num2cell(mat, [1 size(mat,2)]);
    if nargin >= 3
        label = table(label{:}, 'VariableNames', rosettaStone);
    else
        label = table(label{:});
    end
end



%size = length(mat);
%
%for i = 1:size
%    if mat(i) == 1
%        label = i;
%        break;
%    end
%end
