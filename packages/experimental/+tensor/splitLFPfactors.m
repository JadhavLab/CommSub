function Usplit = splitLFPfactors(U, nType)

[nFrequencyByType, nFactors] =  size(U);
Usplit = mat2cell(U, nType, 1);
