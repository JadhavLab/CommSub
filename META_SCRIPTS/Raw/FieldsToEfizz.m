% Grab fields from cgramcnew and place into existing efizz structure
% (rather than from scratch)

grabfields = ["Cavg", "Ctoppair"];

const = option.constants();
animal_list = const.all_animals;
animal = animal_list{1};
for i = 1:length(animal_list)
    animal = animal_list{i};
    efizz = matfile(animal + "spectralBehavior.mat", 'Writable', true);
    assert(ndbFile.exist(animal, 'cgramcnew'), "cgramcnew does not exist for " + animal);
    % cgramcnew = ndb.load(animal, 'cgramcnew', 'inds', [1,2]);
    cgramcnew = ndb.load(animal, 'cgramcnew');
    cgramcnew = nd.cat(cgramcnew, 1, [], 'removeEmpty', true);
    
    tc=cgramcnew.t;
    te=efizz.t;
    assert(all(tc==te), "time vectors do not match for " + animal);

    for j = 1:length(grabfields)
        field = grabfields(j);
        efizz.(field) = cgramcnew.(field);
    end
    
end


