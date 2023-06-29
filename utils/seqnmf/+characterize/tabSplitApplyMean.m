function T = tabSplitApplyMean(T, G)
% Applies mean/mode to groups of T
%
% Ryan Y.

types = arrayfun(@(x) class(table2array(T(1,x))) , 1:size(T,2),'UniformOutput',false); 
names = T.Properties.VariableNames;
t = cell(1,size(T,2));
t = ry_splitapply(size(T,2), @characterize.tab_meanAndMode, T, G);
% for c = 1:size(T,2)
%     try
%         t{c} = cast(t{c}, types{c});
%     catch ME
%         if isequal(types{c}, 'string')
%             t{c} = string(double(t{c}));
%         else
%             error('Fuck')
%         end
%     end
%     assert(~iscategorical(t{c}))
%     %if iscategorical(t{c})
%     %    t{c} = removecats(t{c}, categories(t{c}));
%     %end
% end
t = arrayfun( @(x) cat(1, t{:,x}) ,1:size(T,2), 'UniformOutput', false);
T = table(t{:}, 'VariableNames', names);
