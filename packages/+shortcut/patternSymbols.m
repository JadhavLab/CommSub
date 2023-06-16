function X = patternSymbols(X, level)

if nargin == 1
    level = 1;
end

X = X.replace('theta','\theta').replace('delta','\delta').replace('ripple','SPW-R');

if level > 1
    X = X.replace('-control',' C');
end
