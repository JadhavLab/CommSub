
% ===================================================
% OPTION STRUCT encoding properties of the script run
% ===================================================
% see +option.default() to set default options
if ~exist('Option','var')
    Option = option.defaults(); 
else
    Option = option.setdefaults(Option);
end

Option.analysis.rankRegress    = false;
Option.analysis.factorAnalysis = false;
Option.analysis.timeVarying    = false;
Option.analysis.checks         = true;
Option.save                    = false;

const = option.constants()
for animal = const.all_animals
    Option.animal = animal;
    TheScript
end
