function F = band(data,band)
% SIEVE.LFPBAND picks related lfp bands
%
% Inputs
% ------
% data : the seqnmf datat struct
% band : string/char of the band you would like to index out

switch string(band)
case "delta"
    freqs = [0.1 4];
case "theta"
    freqs = [6 12];
case "beta"
    freqs = [14 20];
case "ripplegamma"
    freqs = [20 40];
case "lowgamma"
    freqs = [30 60];
case "highgamma"
    freqs = [60 100];
case "gamma"
    freqs = [30 100];
case "epsilon"
    freqs = [100 140];
case "ripple"
    freqs = [150 200];
otherwise
    error("Band=" + band + " is not recognized");
end

F = data.f <= freqs(2) & data.f >= freqs(1);
