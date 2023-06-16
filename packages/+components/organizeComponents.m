function [reorganized] = organizeComponents(components, dimSelected, numRhythm, organize_method)

% this function takes the default component matrix (pattern*desired # of
% dimensions) and reorganize it into the sequence we might want to
% visualize

% it could be high/low lfp pattern, such that theta components are followed
% by low-theta components...etc ("stength")
% it could also be ordered by dimension number, such that 1st dim of all
% patterns will be grouped together, then 2nd.... ("dim")

% input: the default component matrix, the dimensions selected for each
% pattern, the desired organize method

% output: the reorganzied component matrix

reorganized = [];
if organize_method == "strength"
    % split into high and low
    numComponents = size(components,1);
    strong_rhythm = components(1:numComponents/2,:);
    weak_rhythm   = components(numComponents/2+1:numComponents, :);
    
    for i = 1:numRhythm
        reorganized = [reorganized; strong_rhythm((i-1)*dimSelected+1:(i-1)*dimSelected+dimSelected,:)];
        reorganized = [reorganized; weak_rhythm((i-1)*dimSelected+1:(i-1)*dimSelected+dimSelected,:)];
    end
    
elseif organize_method == "dim"
    for i = 1:dimSelected
        for j = 1:numRhythm*2
            
            reorganized = [reorganized; components(i+(j-1)*dimSelected,:)];
        end
    end
end


end

