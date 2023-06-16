function [stable_component, learning_component, stable_corr, learning_corr] = ....
    calculateLearningByComponent(component, throwout_times, learning_rate, ...
    stable_rate, learning_times, stable_times)


% calculates the component strength 
% how strong the component strength correlate with learning rate and stable
% rate
%% CURRENTLY UNAVAIAVLE

component = component(~throwout_times,:);

stable_component = component(stable_times,:);
learning_component = component(learning_times,:);

% how to get this?

% norm = @(x) (x - median(x));
% std_norm = @(x) bsxfun(@rdivide, norm(x), std(norm(x)));
% stable_corr = (1/size(stable_times)) * std_norm(stable_rate) * std_norm(stable_component')';
% learning_corr = (1/size(learning_times)) * std_norm(learning_rate) * std_norm(learning_component')';

stable_corr = stable_rate*stable_component;
learning_corr = learning_rate*learning_component;
end
% 
