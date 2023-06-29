function hotencodingInspect(data, hotencodingfields)
% HOTENCODINGINSPECT inspects how the hot encoding performed

% Find relevent fields

if ischar(hotencodingfields)
    hotencodingfields = {hotencodingfields};
end

% Validation plots
nFields = numel(hotencodingfields);
clf;
fcnt = 0;
colors = cmocean('haline',nFields);
for field = hotencodingfields
    fcnt = fcnt + 1;
    plainField = field{1};
    num = 1 + (fcnt-1)*2;
    ax(num) = subplot(nFields, 2, num);
    plot(data.(['orig_' plainField]), 'Color', colors(fcnt,:))
    title([plainField ' raw'])
    num = 2 + (fcnt-1)*2;
    ax(num) = subplot(nFields, 2, num);
    imagesc(data.(plainField)');
    set(gca,'ydir','normal')
    title([plainField ' hot encoded'])
end
linkaxes(findobj(gcf,'type','axes'), 'x')
