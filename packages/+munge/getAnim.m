function [Data, Option] = getAnim(Data, Option, anim)
    animal = nd.fieldGet(Option, 'animal');
    Data = Data(animal == anim, :, :, :, :, :, :);
    Option = Option(animal == anim);
end
