function [outputArg1,outputArg2] = participatingCells(B_, varargin)
%PARTICIPATINGCELLS Summary of this function goes here
%   Detailed explanation goes here

ip = inputParser;
ip.addParameter('color', "black", @isstring);
ip.addParameter('averaged', true, @islogical);
ip.addParameter('do_plot', true, @islogical);
ip.addParameter('normalized', true, @islogical);

ip.parse(varargin{:});
opt = ip.Results;


end

