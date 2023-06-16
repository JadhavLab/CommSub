% function cell_belongs = findBelongingPartition(Patterns, areaPerNeuron )

% this function takes in a pattern struct array and outputs a table of
% cells and the partitions they are involved in

% [cellRegion, cellIDintheRegion, iMethod (method name), partitionnumber_in_this_method]
% table of the above tuple would be produced to facilitate mapping of cells
% to the partitions they are in

% Input: Patterns, the Pattern struct (after assignment of cells)
% Output: the table

[nMethods, nPartition, nDirection, ~] = size(Patterns);


cell_belongs = [];

for m = 1:nMethods
    pfcIndexTracker = 1;
    hpcIndexTracker = 1;
    currMethod = Option.generateH{m};
    for n = 1:numel(areaPerNeuron)
        
        % tracks which neuron is it in the same manner as when assigning
        % partitions
        if areaPerNeuron(n) == "CA1"
            currRegion = "hpc";
            currIndex = hpcIndexTracker;
            hpcIndexTracker = hpcIndexTracker+1;
        else
            currRegion = "pfc";
            currIndex = pfcIndexTracker;
            pfcIndexTracker = pfcIndexTracker+1;
        end
        
        % data structure keeping track of which partition this neuron
        % serves as part of the source or target population
        asSource = cell(nDirection, nPatterns*2);
        asTarget = cell(nDirection, nPatterns*2);
        
        % traverse through the partitions and
        for p = 1:nPartition
            for d = 1:nDirection
                for i = 1*nPatterns*2
                    currPattern = Patterns(m,p,d,i);
                    temp = currPattern.directionality.split("-");
                    currSource = temp(1);
                    currTarget = temp(2);
                    index_source = currPattern.index_source;
                    index_target = currPattern.index_target;
                    
                    if currRegion == currSource
                        if ismember(currIndex, index_source)
                            asSource{d,i} = [asSource{d,i}, currIndex];
                            continue; % because it cannot be the target cell in this particular interaction
                        end
                    end
                    
                    % for a given interaction, a cell cannot be served as
                    % both source and target, even if source and target are
                    % different regions and have overlapping indices.
                    
                    % hpc- in the target, pfc - in the target
                    if currRegion == currTarget
                        if ismember(currIndex, index_target)
                            asTarget{d,i} = [asTarget{d,i}, currIndex];
                        end
                    end
                end
                
                
            end
            
            if currRegion ~= currSource
                assert(isempty([asSource{:}]))
            end
        end
        
        %         should also have: partitions participated in as source,
        %         partitions participated in as target
        tableRow = table(areaPerNeuron(n), m, currMethod, currRegion, currIndex);
        cell_belongs = [cell_belongs; tableRow];
    end
end
% end

