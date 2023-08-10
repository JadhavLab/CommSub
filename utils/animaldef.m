function animalinfo = animaldef(animalname)
	
if ispc
    rootfolder = 'F:/ComSub/';
else
    % Get the user name! and determine user specific root folder
	[~,username] = unix('whoami');
	switch deblank(username)
		case 'ss'				% on laptop
			if exist('/Volumes/Cerebellum/','dir')
				rootfolder = '/Volumes/Cerebellum/';
            elseif exist('/Volumes/Colliculus/','dir')
                rootfolder = '/Volumes/Colliculus/';
            else
				rootfolder = '/Volumes/mordorDATA/';
			end
		case {'jadhavlab','ryoung','root'}		% lab computer
			% particular lab computer must be rsynced such that its ~/Code
			% and ~/Data folders matches Ryan's
            if exist('~/Data/Raw/','dir')
                rootfolder = '~/Data/Raw/';
            elseif exist('/Volumes/Cerebellum/','dir')
                rootfolder = '/Volumes/Cerebellum/';
            elseif exist('/Volumes/Colliculus/','dir')
                rootfolder = '/Volumes/Colliculus/';
            elseif exist('/Volumes/mordorDATA/','dir')
                rootfolder = '/Volumes/mordorDATA/';
            else
                error('animaldef: Cannot find folder on system!');
            end
        otherwise % guess at the folder, and it doesn't exist, then tell the user about the probable issue: animaldef hasn't been setup for this user
            if exist('~/Data/Raw/','dir')
                rootfolder = '~/Data/Raw/';
            elseif ispc && isequal(deblank(username),'users\ryoung')
                rootfolder = 'Z:\';
            else
                error(sprintf(['Haven''nt properly constructured animaldef to find path '...
				'to animal data for username %s'],username));
            end
    end
end


switch animalname

    % mingxin animals
    case 'SM8', animalinfo={'SM8', '/media/ryoung/Colliculus/SM8_direct', 'SM8'};
    case 'MD7', animalinfo={'MD7', '/media/ryoung/Colliculus/MD7_direct', 'MD7'};

    % Goal Maze Animals
    case 'RY7', animalinfo={'RY7',['/media/ryoung/Thalamus' filesep 'ry_GoalCoding_Project/RY7_experiment/RY7_direct'],'RY7'};
    case 'RY9', animalinfo={'RY9',['/Volumes/fornix/RY9_experiment/RY9_direct/'],'RY9'};
    case 'RY16', 
        %animalinfo={'RY16',['/Volumes/FastData/ry_GoalCoding_Project/RY16_experiment' filesep 'RY16_direct'],'RY16'};
        animalinfo={'RY16',['/Volumes/GenuDrive/RY16_direct/'],'RY16'};
    case 'RY22'
        animalinfo={animalname, ['/media/ryoung/Ark/RY22_direct/'], animalname};
    case 'RY20'
        animalinfo={animalname, '/media/ryoung/fornix/RY20_direct', animalname};
    % Phase disrupt animals
    case 'CL1', animalinfo={'CL1',[rootfolder 'Phase/CL1_direct/'],'CL1'};
    case 'SG7', animalinfo={'SG7',[rootfolder 'OdorPlace/SG7Expt/SG7_direct/'],'SG7'};

    % CA3-CA1 animals
    % Right now, just animals we navigated tetrodes down to CA3 to figure out how to get there
    case 'EG1'
	animalinfo = {'EG1', [rootfolder 'CA3-CA1/Efizz/EG1_direct/'] , 'EG1'};

        
    % Single day animals
    case 'ZT2'
        animalinfo = {'ZT2', [rootfolder 'SingleDayExpt/ZT2_direct/'],'ZT2'};
    case 'YD6'
        animalinfo = {'YD6', [rootfolder 'W-Track_WellTrained_EPHYS/YD6_direct/'],'YD6'};
    case 'ER1'
        animalinfo = {'ER1', [rootfolder 'SingleDayExpt/ER1_NEW_direct/'], 'ER1'};
    case 'KL8'
        animalinfo = {'KL8', [rootfolder 'SingleDayExpt/KL8_direct/'], 'KL8'};
    case 'JS7'
        animalinfo = {'JS7', [rootfolder 'W-track_SingleDay/Efizz/JS7_direct/'],'JS7'};
    case 'SJ6'
        animalinfo = {'SJ6', [rootfolder 'W-track_SingleDay/Efizz/SJ6_direct/'],'SJ6'};
    case 'SJ1'
        animalinfo = {'SJ1', [rootfolder 'W-track_SingleDay/Beh/SJ1_direct/'],'SJ1'};
    case 'SJ2'
        animalinfo = {'SJ2', [rootfolder 'W-track_SingleDay/Beh/SJ2_direct/'],'SJ2'};
    case 'SJ3'
        animalinfo = {'SJ3', [rootfolder 'W-track_SingleDay/Beh/SJ3_direct/'],'SJ3'};
    case 'SJ4'
        animalinfo = {'SJ4', [rootfolder  'W-track_SingleDay/Beh/SJ4_direct/'],'SJ4'};
    case 'JS12'
        animalinfo = {'JS12', [rootfolder 'SingleDayExpt/JS12_direct/'],'JS12'};
    case 'JS13'
        animalinfo = {'JS13', [rootfolder 'SingleDayExpt/JS13_direct/'],'JS13'};
    case 'JS14'
        animalinfo = {'JS14', [rootfolder 'SingleDayExpt/JS14_direct/'],'JS14'};
    case 'JS15'
        animalinfo = {'JS15', [rootfolder 'SingleDayExpt/JS15_direct/'],'JS15'};
    case 'JS17'
        animalinfo = {'JS17', [rootfolder 'SingleDayExpt/JS17_direct/'],'JS17'};
    case 'JS21'
        animalinfo = {'JS21', [rootfolder 'SingleDayExpt/JS21_direct/'],'JS21'};
        
    % New 8day HC-PFC animals
    case 'SJ5'
        animalinfo = {'SJ5', [rootfolder 'HP_8dayExpt/SJ5_direct/'],'SJ5'};    
        
    % Ripple interruption, old animals
    case 'sjc'
        animalinfo = {'sjc', [rootfolder 'RippleDisruption_all/RippleDisruption/sjc_direct/'], 'sjc'};
    case 'RE1'
        animalinfo = {'RE1', [rootfolder 'RippleDisruption_all/RippleDisruption/RE1_direct/'], 'RE1'};
    case 'RNa'
        animalinfo = {'RNa', [rootfolder 'RippleDisruption_all/RippleDisruption/RNa_direct/'], 'RNa'};
    case 'RNb'
        animalinfo = {'RNb', [rootfolder 'RippleDisruption_all/RippleDisruption/RNb_direct/'], 'RNb'};
    case 'RNc'
        animalinfo = {'RNc', [rootfolder 'RippleDisruption_all/RippleDisruption/RNc_direct/'], 'RNc'};
    case 'RNd'
        animalinfo = {'RNd', [rootfolder 'RippleDisruption_all/RippleDisruption/RNd_direct/'], 'RNd'};
    case 'RCa'
        animalinfo = {'RCa', [rootfolder 'RippleDisruption_all/RippleDisruption/RCa_direct/'], 'RCa'};
    case 'RCb'
        animalinfo = {'RCb', [rootfolder 'RippleDisruption_all/RippleDisruption/RCb_direct/'], 'RCb'};
    case 'RCc'
        animalinfo = {'RCc', [rootfolder 'RippleDisruption_all/RippleDisruption/RCc_direct/'], 'RCc'};
    case 'RCd'
        animalinfo = {'RCd', [rootfolder 'RippleDisruption_all/RippleDisruption/RCd_direct/'], 'RCd'};
    case 'REc'
        animalinfo = {'REc', [rootfolder 'RippleDisruption_all/RippleDisruption/REc_direct/'], 'REc'};
    case 'REd'
        animalinfo = {'REd', [rootfolder 'RippleDisruption_all/RippleDisruption/REd_direct/'], 'REd'};
    case 'REe'
        animalinfo = {'REe', [rootfolder 'RippleDisruption_all/RippleDisruption/REe_direct/'], 'REe'};
    case 'REf'
        animalinfo = {'REf', [rootfolder 'RippleDisruption_all/RippleDisruption/REf_direct/'], 'REf'};
    case 'REg'
        animalinfo = {'REg', [rootfolder 'RippleDisruption_all/RippleDisruption/REg_direct/'], 'REg'};
    case 'REh'
        animalinfo = {'REh', [rootfolder 'RippleDisruption_all/RippleDisruption/REh_direct/'], 'REh'};

    % Hippocampal-prefrontal animals
    case 'HPa'
        animalinfo = {'HPa', [rootfolder 'HP_8dayExpt/HPExpt/HPa_direct/'], 'HPa'};
    case 'HPb'
        animalinfo = {'HPb', [rootfolder 'HP_8dayExpt/HPExpt/HPb_direct/'], 'HPb'};
    case 'HPc'
        animalinfo = {'HPc', [rootfolder 'HP_8dayExpt/HPExpt/HPc_direct/'], 'HPc'};
    case 'Nadal'
        animalinfo = {'Ndl', [rootfolder 'HP_8dayExpt/HPExpt/Ndl_direct/'], 'Ndl'};
    case 'Rosenthal'
        animalinfo = {'Rtl', [rootfolder 'HP_8dayExpt/HPExpt/Rtl_direct/'], 'Rtl'};
    case 'Borg'
        animalinfo = {'Brg', [rootfolder 'HP_8dayExpt/HPExpt/Brg_direct/'], 'Brg'};
        
    % UH  OH -- animal name not recognized ...
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end

if ~exist(animalinfo{2},'dir')
    warning(['Directory ' animalinfo{2} ' not found! Creating...']);
    mkdir(animalinfo{2});
end

%fprintf('Using animal info at %s \n', animalinfo{2});
