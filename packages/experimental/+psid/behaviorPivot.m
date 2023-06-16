function pivot = behaviorPivot(behavior_name, varargin)

switch behavior_name
    case 'directional_lindist'
        pivot = 0;
    case 'trajbound'
        pivot = 0.5;
    case 'rewarded'
        pivot = 0.5;
    case 'lindist'
        pivot = 0.5;
    case 'leftright'
        pivot = 0.5;
    case 'X'
        pivot = @median;
    case 'Y'
        pivot = @median;
    otherwise
        error('Unknown');
end

if numel(varargin) > 0 && ~isnumeric(pivot)
    pivot = pivot(varargin{1});
end
