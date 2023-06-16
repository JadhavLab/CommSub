function [pfcFR, hpcFR] = separateFiringRate(avgFR, areaPerNeuron)
% Separates out the average firing rate of each neuron

pfcFR = [];
hpcFR = [];

for i = 1:numel(areaPerNeuron)
    if areaPerNeuron(i) == "PFC"
        pfcFR = [pfcFR, avgFR(i)];
    else
        hpcFR = [hpcFR, avgFR(i)];
    end
end

end

