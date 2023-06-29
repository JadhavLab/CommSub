function [newdata, newfields, varargout] = splitdatatofields(data, fields, varargin)

ip = inputParser;
ip.addParameter('tensor',false);
ip.addParameter('structfield','data');
ip.addParameter('noPhi',false);
ip.addParameter('repelemMapping', containers.Map)
ip.parse(varargin{:});
opt = ip.Results;

fieldcount = numel(fields);
phiSinInFields = find(ismember(fields, 'phi_sin'));
phiCosInFields = find(ismember(fields, 'phi_cos'));
if isempty(phiSinInFields); phiSinInFields = 0; phiCosInFields = 0; end


if isstruct(data)
    return_struct = true;
    D = data.(opt.structfield);
else
    return_struct = false;
    D = data;
end

if ndims(D) == 2
    [N,T] = size(D);
    % Data chunks
    if isstruct(data)
        n_sizes = cellfun(@(x) size(data.(x)), fields, 'UniformOutput', false);
        n_sizes = cat(1,n_sizes{:});
        % Idetify the N dimension
        n_sizes = n_sizes(:,2);
        n_sizes = n_sizes(:);
    else
        n_sizes = repmat(N/numel(fields), numel(fields), 1);
    end
    % Split chunks into varargout
    varargout= mat2cell(D, n_sizes, T);
elseif ndims(D) == 3
    [N,K,T] = size(D);
    % Data chunks
    if isstruct(data)
        n_sizes = cellfun(@(x) size(data.(x)), fields, 'UniformOutput', false);
        n_sizes = cat(1,n_sizes{:});
        % Idetify the N dimension
        n_sizes = n_sizes(:,2);
        n_sizes = n_sizes(:);
    else
        n_sizes = repmat(N/numel(fields), numel(fields), 1);
    end
    % Split chunks into varargout
    varargout= mat2cell(D, n_sizes, K, T);
else
    error('Invalid dimensions of data')
end

% Map it from cos/sin back into circular domain
if phiSinInFields ~= -1 && ~opt.noPhi
    phi_sin = varargout{phiSinInFields};
    phi_cos = varargout{phiCosInFields};
    % Normalize by radius just in case we're dealing with a W which might have pushed the sin and cos components arbirtarily off the unit circle
    phi_sin(isnan(phi_sin)) = 0;
    phi_cos(isnan(phi_cos)) = 0;
    R = sqrt(phi_sin.^2  + phi_cos.^2);
    R = nanmax(R,[],[1,2]);
    phi_sin = bsxfun(@rdivide, phi_sin, R);
    phi_cos = bsxfun(@rdivide,phi_cos, R);
    phi = atan(phi_sin./phi_cos);
    phi(phi_cos < 0 & phi_sin > 0) = phi(phi_cos < 0 & phi_sin > 0) + pi; % map anything in the negativ cosine quadrants to pi from the result...atan in (-pi/2, pi/2)
    phi(phi_cos < 0 & phi_sin < 0) = phi(phi_cos < 0 & phi_sin < 0) - pi; % map anything in the negativ cosine quadrants to pi from the result...atan in (-pi/2, pi/2)
    phi(phi_cos == 0 & phi_sin > 0) = pi/2;
    phi(phi_cos == 0 & phi_sin < 0) = -pi/2;

    varargout([phiSinInFields, phiCosInFields]) = [];
    varargout{end+1} = phi;
    if isstruct(data)
        newdata = data;
        newdata.data = cat(1, varargout{:});
    else
        newdata = cat(1, varargout{:});
    end
    if isstruct(newdata)
        newdata.fields = fields;
    end
    newfields = fields;
    newfields([phiSinInFields, phiCosInFields]) = [];
    newfields{end+1} = 'phi';
else
    newdata = data;
    newfields = fields;
end

% If return a tensor type!
if opt.tensor
    if return_struct
        newdata.data = cat(3, varargout{:}); % Return tensor instead of matrix
    else
        newdata = cat(3, varargout{:}); % Return tensor instead of matrix
    end
end
