function T = string2categorical(T)

for field = string(T.Properties.VariableNames)
    if isstring(T.(field))
        T.(field) = categorical(T.(field));
    end
end
