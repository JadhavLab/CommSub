function fieldcolor(field, color)

switch field
case {'S1','S2','C'}
    if nargin == 2
        cmocean(color)
    else
        cmocean('thermal')
    end
case {'phi'}
    cmocean('phase')
end
