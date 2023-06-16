function animateComponent(beh, component_matrix, k, varargin)

ip = inputParser;
ip.addParameter('kws',{'linestyle','none', 'markersize', 12});
Opt = ip.Results;

fig('animated component');

component_matrix = (component_matrix - min(component_matrix))./(max(component_matrix)-min(component_matrix));
component_matrix = round(1000*component_matrix)

plot(beh.x, beh.y, ':k');
hold on;

colors = crameri('bam', 1000);
a = animatedline(Opt.kws{:});
for i = 1:height(beh)
    clearpoints(a);
    value = component_matrix(k,i);
    set(a,'MarkerFaceColor', color(value,:));
    addpoints(a, beh.x(i), beh.y(i));
    drawnow limitrate
end

