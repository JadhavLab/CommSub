function varargout = behaviorColor(behavior_name, varargin)

switch behavior_name
    case 'directional_lindist'
        colors = crameri('bukavu', varargin{:});
        colors_midpoint = round(size(colors,1)/2);
        colors(colors_midpoint+1:end) = colors(end:-1:colors_midpoint);
    case 'trajbound'
        colors = crameri('vik', varargin{:});
        colors = binaryColor(colors, 0.25);
    case 'rewarded'
        colors = crameri('bam', varargin{:});
        colors = binaryColor(colors, 0.25);
    case 'lindist'
        colors = crameri('tofino', varargin{:});
    case 'leftright'
        colors = crameri('tokyo', varargin{:});
        colors = binaryColor(colors, 0.2);
    case 'X'
        colors = crameri('tofino', varargin{:});
    case 'Y'
        colors = crameri('turku', varargin{:});
    otherwise
        error('Unknown');
end

if nargout == 0
    colormap(colors);
else
    varargout{1} = colors;
end

function colors = binaryColor(colors, fraction)

    colors_midpoint = round(size(colors,1)/2);
    fraction_up = round(size(colors,1) * fraction);
    fraction_down = round(size(colors,1) * (1-fraction));

    C = colors;
    N = size(colors(colors_midpoint:end, :),1);
    colors(colors_midpoint:end, :) = repmat(C(fraction_up,:),N,1);
    N = size(colors(1:colors_midpoint-1, :),1);
    colors(1:colors_midpoint-1, :) = repmat(C(fraction_down,:),N,1);
