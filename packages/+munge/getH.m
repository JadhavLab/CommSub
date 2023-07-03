function [Data, Option] = getH(Data, Option, genh)
    genH = nd.fieldGet(Option, 'generateH');
    genH = replace(genH, "fromCoherence  fromRipTimes", "coh");
    genH = replace(genH, "fromSpectra  fromRipTimes", "spec");
    genH = replace(genH, "fromWpli  fromRipTimes", "wpli");
    Data = Data(genH == genh, :, :, :, :, :, :);
    Option = Option(genH == genh);
end
