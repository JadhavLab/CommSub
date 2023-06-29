function figure_ = seqnmf_plotdata(data, varargin)
% SEQNMF_PLOTDATA Yields my collection of fancy plots for seqnmf using the data
% in my seqnmf struct
%
% Inputs
% ------
%
% Optional inputs
% ---------------
%
% Output
% ------


ip = inputParser;
ip.addParameter('fields',{}, @iscell)
ip.addParameter('plotType', "rawdata")
ip.addParameter('plotOpt', struct())
ip.parse(varargin{:})
opt = ip.Results;

setting = @(x) isfield(opt.plotOpt, x) && plotOpt.x == true;
if isempty(opt.fields)
    opt.fields = data.fields;
end
nFields = numel(opt.fields);

% Corrections
if contains(opt.fields,'phi')
    data = seqnmf_unpackfields(data, opt.fields, 'noPhi', true);
end

for plotType = opt.plotType
    switch opt.plotType

    %%%%%%%%%%%%%%%%% RAW DATA PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "rawdata"
        figure_ = fig('transformeddata');
        if fieldtest('separate')
            ax = subplot(1,1,1);
            nestable on;
            nfield = 0;
            for field = data.fields
                nfield = nfield + 1;
                axx = nestplot(nFields,1,nfield);

                imagesc(data.t, data.f, data.(field{1}));
                colorbar
                fieldcolor(field{1});

                % Xaxis
                if nfield < nFields
                    xticks([]);
                else
                    xlabel('Time (s)')
                end

                % Yaxis
                switch field{1}
                case {'S1','S2','C','wpli'}
                    ylabel('Frequency (hz)')
                    seq.frequencyaxis(data)
                case {'lindist','trajdist'}
                    ylabel('Region')
                    seq.regionaxis(data, field{1})
                end

            end
        else
            imagesc(data.data);
        end
        set(gca,'ydir','normal')

    %%%%%%%%%%%%%%%%% Ws Visualization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'wplot'

        figure_ = fig('wplot');
        clf

		% PLOT Ws
        ax = subplot(1,2,1);
        nestable('on',ax);
        nfield = 0;
        faxes = [];
        raxes = [];
        for field = fieldnames(data)
            if iscell(field)
                field = field{1};
            end
            if contains(field, 'W_')
                wfield = field(strfind(field,'W_')+2:end);
                disp(field)
            else
                continue
            end
            nfield=nfield+1;
            axx = nestplot(nFields, 1, nfield);
            if ~isequal(wfield,'phi')
                imagesc(1:size(data.(field), 3), data.f(1,:), imgaussfilt(flipud(squeeze(data.(field)(:,1,:))),2));
            else
                imagesc(1:size(data.(field), 3), data.f(1,:), flipud(squeeze(data.(field)(:,1,:))));
            end
            set(gca,'ydir','normal')
            xlim([1,size(data.(field),3)])

            % Xaxis
            if nfield < nFields
                xticks(axx, []);
            else
                xlabel(axx, 'Time (s)')
            end

            % Yaxis
            switch wfield
            case {'S1','S2','C','wpli','pfcSPca1LFP_C','ca1SPpfcLFP_C','pfcSPca1LFP_phi','ca1SPpfcLFP_phi','phi'}
                ylabel(axx,sprintf('%s\n Frequency (hz)', wfield));
                seq.frequencyaxis(data)
                fieldcolor(wfield)
                colorbar
                faxes = [faxes, axx];
            case {'lindist','trajdist'}
                ylabel(axx,sprintf('%s\n Regiont', wfield));
                seq.regionaxis(data, wfield)
                colorbar
                raxes = [raxes, axx];
            end
            ylim([-inf, inf])
        keyboard
        end
        if ~isempty(faxes); linkaxes(faxes, 'y'); end
        if ~isempty(raxes); linkaxes(raxes, 'y'); end
        if ~isempty([raxes, faxes]); linkaxes([faxes, raxes], 'x'); end
        nestable('off',ax)

		% -- PLOT FULL DATA --
        subplot(1,2,2)
		%Versus full data
        imagesc(1:size(data.W, 3), data.f(1,:), flipud(squeeze(data.W(:,1,:))));
		seq.fieldfrequencyaxis(data, 5)
        cmocean('haline')

    %%%%%%%%%%%%%%%%% Ws Visualization verus DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	case 'W-data'

		figure_ = seqnmf_plotdata(data, 'fields', data.fields, 'plotType', 'wplot');
		subplot(1,2,2)
        imagesc(data.t, data.f(1,:), flipud(squeeze(data.data(:,:))));
		seq.fieldfrequencyaxis(data, 5);
		X = xlim;
		xlim([range(X)*0.12+X(1), range(X)*0.1205+X(1)]);
        cmocean('-gray')

    end

end
