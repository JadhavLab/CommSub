function [data, newfields] = seqnmf_unpackfields(data, fields, varargin)
% Returns data mapped back into it's natural domain

% ---------------------
% Parse optional inputs
% ---------------------
ip = inputParser();
ip.KeepUnmatched = true;
ip.addParameter('undoPackTransform', false);
ip.addParameter('noPhi', true);
ip.parse(varargin{:});
pack_fields = ip.Unmatched;
opt = ip.Results;

fieldnum = numel(fields);
if ismember('phi_cos', fields) && ~opt.noPhi
    fieldnum=fieldnum-1;
end

% ---------
% Unpack W
% --------
D = cell(1,fieldnum);
[Wprime, newfields, D{:}] = splitdatatofields(data, fields, 'structfield', 'W', 'noPhi', opt.noPhi);
if opt.undoPackTransform
    [Wprime, newfields, D{:}] = undopacktransform(D, newfields, pack_fields{:});
end
dcnt=0;
for field = newfields
    dcnt = dcnt + 1;
    field = field{1};
    data.(['W_' field]) = D{dcnt};
end
%data.Wprime = Wprime;

% ---------
% Unpack WH
% ---------
if ~isfield('WH', data)
    data.WH = helper.reconstruct(data.W, data.H)
end
[WHprime, newfields, D{:}] = splitdatatofields(data, fields, 'structfield', 'WH', 'noPhi', opt.noPhi);
if opt.undoPackTransform
    [WHprime, newfields, D{:}] = undopacktransform(D, newfields, pack_fields{:});
end
dcnt = 0;
for field = newfields
    dcnt = dcnt + 1;
    field = field{1};
    data.(['WH_' field]) = D{dcnt};
end
%data.WHprime = WHprime;
