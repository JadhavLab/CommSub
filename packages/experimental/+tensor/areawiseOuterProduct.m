function Z = areawiseOuterProduct(X, Y,varargin)
%
% creates a tensor outer product of two area-specific tensors, where time slices
% of a given output tensor are the outer product pieces that would be averaged
% to obtain a B regression matrix.

ip = inputParser;
ip.addParameter('big',false)
ip.parse(varargin{:})
Opt = ip.Results;

%keyboard
if isstring(Y) && ~iscell(X)

    if isstruct(X)
        X = ful(X);
    end

    uVals = unique(Y);
    Z = cell(numel(uVals),1);
    for i = 1:numel(uVals)
        Z{i} = X(Y == uVals(i),:,:);
    end
    X = Z{1};
    Y = Z{2};
end

if ~iscell(X) && ~iscell(Y)
    notACell = true;
    X = {X};
    Y = {Y};
else
    notACell = false;
end

center = true;

Z = cell(size(X));
for pattern = 1:numel(X)
    x = X{pattern};
    y = Y{pattern};
    if isstruct(x)
        x = ful(x);
    end
    if isstruct(y)
        y = ful(y);
    end
    if center % remove mean and divide by sigma
        x = (x - nanmean(x(:,:),2))./nanstd(x(:,:),1,2);
        y = (y - nanmean(y(:,:),2))./nanstd(y(:,:),1,2);
        x(isinf(x)) = 0;
        y(isinf(y)) = 0;
    end
    nxNeuron = size(x,1);
    nyNeuron = size(y,1);
    nTime = size(x,2);
    nTrial = size(x,3);
    if ~Opt.big
        z = zeros(nxNeuron, nyNeuron, nTime, nTrial, 'single');
        for time = 1:nTime
            for trial = 1:nTrial
                z(:,:,time,trial) = x(:,time,trial) * y(:,time,trial)';
            end
        end
    else
        x=single(x);
        y=single(y);
        x=shiftdim(x,-1);
        y=shiftdim(y,-1);
        x=fmt(permute(x,[2,1,3,4]));
        y=fmt(permute(y,[1,2,3,4]));
        %one  of
        z = single(bsxfun(@times, x, y));
    end
    z = fmt(z); % Place tensor in tensorlab's format
    Z{pattern} = z;

end

if notACell
    Z = Z{1};
end
