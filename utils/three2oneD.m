function oneD = three2oneD(x,y,z,x_max,y_max,z_max)

% converts a 3D indexing into 1D based on MATLAB indexing convention
% (1-based)

oneD = x-1 + (y-1)*x_max + (z-1)*y_max*x_max +1;

end

