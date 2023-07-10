addpath(genpath('~/Code/analysis/'))
addpath(genpath(codedefine()))

const  = option.constants();
animal = const.all_animals(1);
for animal = "ZT2"
    animal = string(animal);
    disp(animal)
    ndb.load(animal, 'cellinfo')
    ndb.load(animal, 'tetinfo')
    cellinfo = ndb.toNd(cellinfo);
    tetinfo  = ndb.toNd(tetinfo);
    c=coding.table.info(animal, 'cellinfo');
    t=coding.table.info(animal, 'tetinfo');
    % PICK CA1 tet
    % find all rip tets
    rip_tets = t(t.descrip == "riptet", :);
    % groupby tetrode and sum numcells
    rip_tets = groupsummary(rip_tets, {'tetrode'}, 'mean', 'numcells');
    % fine top 3 tetrodes with most cells
    top3 = sortrows(rip_tets, 'mean_numcells', 'descend');
    top3 = top3(1:3, :);
    % PICK PFC tet
    % find all pfc tets
    pfc_tets = t(t.area == "PFC", :);
    % groupby tetrode and sum numcells
    pfc_tets = groupsummary(pfc_tets, {'tetrode'}, 'mean', 'numcells');
    % fine top 3 tetrodes with most cells
    top3_pfc = sortrows(pfc_tets, 'mean_numcells', 'descend');
    % Loadup the eeg data for each of these hits and average them
    ca1_tet = top3.tetrode(1:3);
    pf_tet = top3_pfc.tetrode(1:3);
    e      = doavg(animal, ca1_tet, pf_tet, 'eeg');
    assert(isfield(e, 'D'))
    theta  = doavg(animal, ca1_tet, pf_tet, 'theta');
    delta  = doavg(animal, ca1_tet, pf_tet, 'delta');
    ripple = doavg(animal, ca1_tet, pf_tet, 'ripple');
    assert(isfield(ripple, 'D'))
    for i = 1:numel(e)
        e(i).theta = theta(i);
        e(i).delta = delta(i);
        e(i).ripple = ripple(i);
    end
    % load the eeg data
    e = nd.toNdBranch(e);
    % save the eeg data
    ndb.save(e, animal, 'avgeeg', 'indices', 0);
end

clear e theta delta ripple avgeeg eeg

% animal = const.all_animals(1);
% %% load the eeg data
% for animal = const.all_animals
%     m=load(animal + "avgeeg.mat");
%     m.avgeeg = ndb.toNd(m.avgeeg);
%     for i = 1:numel(m.avgeeg)
%         m.avgeeg(i).D = squeeze(m.avgeeg(i).D);
%     end
%     utils.type.struct(m)
%     save(1animal + "avgeeg.mat", '-struct', 'm', '-v7.3', '-nocompression')
% end

% %% load the eeg data
% and add area field
for animal = progress(const.all_animals, 'Title', 'adding area field')
    m=load(animal + "avgeeg01.mat");
    m.avgeeg = ndb.toNd(m.avgeeg);
    inds = ndb.indicesMatrixForm(m.avgeeg);
    for ind = inds'
        item = ndb.get(m.avgeeg, ind);
        if ind(end) == 1
            item.area = "CA1";
        else
            item.area = "PFC";
        end
        m.avgeeg=ndb.set(m.avgeeg, ind, item);
    end
    ndb.save(m.avgeeg, animal, 'avgeeg', 'indices', 0);
end

!pushover-cli "done with avgeeg"

function e =  doavg(animal, ca1_tet, pf_tet, eegprop)
    assert(ndbFile.exist(animal, eegprop))
    eeg = ndb.load(animal, eegprop, 'indices', 1);
    eeg = ndb.toNd(eeg);
    ca1 = eeg(:,:,ca1_tet);
    pfc = eeg(:,:,pf_tet);
    clear eeg
    % average the eeg data
    e = ca1(1,:,1);
    for i = 1:size(ca1,2)
        data = {ca1(1,i,:).data};
        data = cat(3, data{:});
        mudata = mean(data, 3);
        e(:,i,1).data = mudata;
        e(:,i,1).D = squeeze(data);
        e(:,i,1).tet = ca1_tet(i);
    end
    % get e for pfc
    % average the pfc data
    e2 = pfc(1,:,1);
    for i = 1:size(pfc,2)
        data = {pfc(1,i,:).data};
        data = cat(3, data{:});
        mudata = mean(data, 3);
        e2(:,i,1).data = mudata;
        e2(:,i,1).D = squeeze(data);
        e2(:,i,1).tet = pf_tet(i);
    end
    % % plot the eeg data
    % figure
    % tiledlayout(2,1)
    % nexttile
    % plot(squeeze(e(:,1,:).D), 'LineWidth', 0.5, 'LineStyle', ':')
    % hold on
    % plot(e(1).data, 'LineWidth', 2, 'Color', 'k')
    % title(sprintf('tet %d', ca1_tet))
    % % plot the pfc data
    % nexttile
    % plot(squeeze(e2(:,1,:).D), 'LineWidth', 0.5, 'LineStyle', ':')
    % hold on
    % plot(e2(1).data, 'LineWidth', 2, 'Color', 'k')
    % title(sprintf('tet %d', pf_tet))
    % linkaxes(findall(gcf, 'type', 'axes'), 'x')
    % if ~exist(figuredefine("avgeeg"), 'dir')
    %     mkdir(figuredefine("avgeeg"))
    % end
    % savefig(fullfile(figuredefine("avgeeg"), sprintf('%s', animal)))
    e = cat(3, e, e2);
    e = util.type.castefficient(e, 'compressReals', true);
end


