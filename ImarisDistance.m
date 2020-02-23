function [output_distances] = ImarisDistance(varargin)
[file,path,idx] = uigetfile('.csv','Choose the csv data file from imars');
DataTable = readtable([path,file]);
x = DataTable.PositionX;
y = DataTable.PositionY;
z = DataTable.PositionZ;

x_dif = x-x';
y_dif = y-y';
z_dif = z-z';

distances=(x_dif.^2+y_dif.^2+z_dif.^2).^0.5;
distances_nan = (distances./distances).*distances;

output_distances = min(distances_nan,[],2);

if length(varargin) == 1
    histogram(output_distances, varargin{1})
else
    histogram(output_distances)
end

xlabel('distances µm')
ylabel('counts')



