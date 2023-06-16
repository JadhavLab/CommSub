function X = generateH(X)
% Creates shortcut phraes for the long-windede genH strings

if iscell(X)
    X = string(X);
end

X = X.replace('fromFilteredEEG','hilbert').replace('fromWpli','wpli').replace('  ','').replace('fromRipTimes','+riptimes').replace('fromSpectra','fft').replace('fromCoherence','coherence');
if all(contains(X,'riptimes'))
    X = X.replace('+riptimes','');
end
