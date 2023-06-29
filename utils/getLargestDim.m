function largest = getLargestDim(Patterns)

% Function to get the largest OptDim within the Patterns struct

[nDirection, nPattern] = size(Patterns);
largest = 0;
for i = 1:nDirection
    for j = 1:nPattern
        curr = Patterns(i,j).rankRegress;
        currDim = curr.optDimReducedRankRegress;
        if currDim>largest
            largest = currDim;
        end
    end
end
end

