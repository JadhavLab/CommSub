paramsets = { {'run',   10, 2, 'traj',       ""},             ... attempt linearized traj via seqnmf
    {'run',   10, 3, 'traj',       "Events-based"}, ... 3 patterns at 10 sec
    {'run',   4,  2, 'traj',       "Events-based"},             ... attempt linearized traj via seqnmf {'run',   4,  5, 'traj',       "Events-based"}, ... 5 patterns {'run',   4,  5, 'speccoh',    "Events-based"}, ... 5 patterns
    {'run',   4,  9, 'speccoh',    "Events-based"}, ... 5 patterns
    {'run',   4,  9, 'trajNoWPLI', "Events-based"}, ... 8 patterns : can i keep finding past a "reasonable" point
    {'run',   1,  9, 'speccoh',    "Events-based"}, ... 8 patterns : can i keep finding past a "reasonable" point
    {'run',   1,  9, 'trajNoWPLI', "Events-based"}, ... 8 patterns : can i keep finding past a "reasonable" point
    {'sleep', 3,  4, 'speccoh',   "Events-based"} ...
    {'run', 4,  9, 'traj',   "Events-based"} ...
    {'run', 4,  9, 'speccoh',   "Events-based"} ...
    {'run',   0.5,  9, 'speccoh',    "Events-based"}, ... 8 patterns : can i keep finding past a "reasonable" point
    };
paramsetstable = cat(1,paramsets{:});
paramsetstable = num2cell(paramsetstable,1);
paramsetstable = table(paramsetstable{:}, 'VariableNames', {'epoch_type','timescale','K','type','seqStyle'});
for variable = paramsetstable.Properties.VariableNames
    try
        paramsetstable.(variable{1}) = cell2mat(paramsetstable.(variable{1}));
    catch ME
        paramsetstable.(variable{1}) = string(paramsetstable.(variable{1}));
    end
end

