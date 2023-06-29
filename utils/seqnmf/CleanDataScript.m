for animal = {'HPa','HPb','HPc'}
    pushd(animal{1} + "_data");
    %Clean cellinfo
    cleanData(pwd, 'datatype', 'cellinfo', 'acceptableTags', {'spikewidth','meanrate','numspikes','area'})
    popd;
end
