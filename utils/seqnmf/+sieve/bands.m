function filter = bands(data)
% SIEVE.BANDS returns a filter structure for every posisble band (complete)

filter = [];
for band = sieve.bandset()
    filter.(char(band)) = sieve.band(data,band);
end
