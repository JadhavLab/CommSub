function behaviorMatrix = generateBehavior(animal)

% this function generates the behavior matrix of an animal {time, x, y}

load (animal+"pos01.mat");
pos = ndb.toNd(pos);
postime = [];
pos_x = [];
pos_y = [];

for i = 2:2:numel(pos)
    postime = [postime, pos(i).data(:,1)'];
    pos_x = [pos_x, pos(i).data(:,2)'];
    pos_y = [pos_y, pos(i).data(:,3)'];
end

behaviorMatrix = [postime', pos_x', pos_y'];

end

