
T = query.getPatternTable(Patterns);
patternType = T.patternType;
patternType( ~contains(patternType, '-' ), : ) = patternType( ~contains(patternType, '-' ), : ) + "-";
columns = patternType.split('-');
columns(columns(:,2)=="",2) = "pattern activity";
T.patternAbstract = columns(:,1);
T.control = columns(:,2);

clf
fig('Dimensionality comparison')
subset = T.directionality=="hpc-pfc";
