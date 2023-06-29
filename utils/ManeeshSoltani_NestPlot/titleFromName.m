function titleFromName(F, g)
%TITLEFROMNAME Assigns .Name of the figure as an SGtitle

if nargin == 0
    F = gcf;
end
if nargin == 2 && ~isempty(g)
    useGramm = true;
else
    useGramm = false;
end

name = string(get(gcf, 'Name'));
if useGramm
    g.set_title(name);
else
    sgtitle(F, name);
end
