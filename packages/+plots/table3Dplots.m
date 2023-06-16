function  table3Dplots(x,y,z,varargin)

% this function plots the relation among three 


ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('demarcatingLine', true, @islogical);
ip.parse(varargin{:});
opt = ip.Results;

scatter3(x,y,z)
end

