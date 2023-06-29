function cleanData(folder, varargin)

ip = inputParser;
ip.addParameter('datatype','cellinfo')
ip.addParameter('acceptableTags',{'spikewidth','meanrate','numspikes','area'});
ip.parse(varargin{:});
opt = ip.Results;
acceptableTags = opt.acceptableTags;

pushd(folder)


for D = dir(['*' opt.datatype '*'])'

    load(D.name);


    days   = 1:8;
    epochs = 1:8;
    tets   = 1:32;
    cells  = 1:100;


    p = ProgressBar(numel(days)*numel(epochs)*numel(tets)*numel(cells), 'Title', ['Cleaning data in ' D.folder]);
    BF = {};
    for d = days
    for e = epochs
    for t = tets
    for c = cells
        try
            cdata = cellinfo{d}{e}{t}{c};
            if isempty(cdata)
                continue
            end
            fields = fieldnames(cdata);
            badfields = setdiff(fields, acceptableTags);
            if ~isempty(badfields) % Bad fields found
                % Remove bad fields
                cellinfo{d}{e}{t}{c} = rmfield(cdata, badfields);
                BF = union(badfields,BF);
            end
        catch ME
        end
        p.step([],[],[])
    end
    end
    end
    end
    disp("Bad fields detected = > ");
    disp(strvcat(BF))
    p.release()

    % Show changed cell info
    open cellinfo

    %keyboard

    % Save cell info
    save(D.name,'cellinfo');
    disp('Saved')
end

% Go back to folder
popd;
