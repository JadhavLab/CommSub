function X = directionality(X)
% shortcut directionatlity names

if iscellstr(X)
    X = string(X);
end

X = X.replace('hpc-hpc','HH').replace('hpc-pfc','HP');
