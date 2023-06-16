function T = abbreviateGenH(T)
% Abbreviates generateH method names

T.generateH = T.generateH.replace('fromFilteredEEG','hilbert').replace('fromWpli','wpli').replace('  ','').replace('fromRipTime','+riptimes').replace('fromSpectra','fft').replace('fromCoherence','coherence');
if all(contains(T.generateH,'riptimes'))
    T.generateH = T.generateH.replace('+riptimes','');
end
