function [lon_l,lat_l] = distcoordinate(lon,lat)
%DISTCOORDINATE Create distance coordinate system from lat/lon grid points
%   [lon_l,lat_l] = distcoordinate(lon,lat)
%
%   Creates distance coordinate system in meters from longitude/latitude grid
%   for use in gradient, divergence, and curl calculations
%
%   Inputs:
%       lon - longitude range (n*1 vector)
%       lat - latitude range (n*1 vector)
%
%   Outputs:
%       lon_l - longitudinal distance coordinates (meters)
%       lat_l - latitudinal distance coordinates (meters)

nlon = size(lon,1);
nlat = size(lat,1);
distlat = zeros(nlon,nlat-1); 
distlon = zeros(nlon-1,nlat);

% Calculate latitudinal distances between grid points
for i = 1:nlat-1
    [distlat(:,i), phaseangle] = sw_dist([lat(i,1) lat(i+1,1)], [lon(1,1) lon(1,1)], 'km');
end

% Calculate longitudinal distances between grid points  
for i = 1:nlat
    [distlon(:,i), phaseangle] = sw_dist([lat(i,1) lat(i,1)], [lon(1,1) lon(2,1)], 'km');
end

% Initialize coordinate matrices
lat_l = zeros(nlon,nlat);
lon_l = zeros(nlon,nlat);

% Build cumulative latitudinal distance coordinates
for i = 1:nlon
    lat_l(i,:) = cumsum(cat(2, 0, distlat(i,:))) * 1000;  % Convert km to meters
end

% Build cumulative longitudinal distance coordinates
for i = 1:nlat
    lon_l(:,i) = cumsum(cat(1, 0, distlon(:,i))) * 1000;  % Convert km to meters
end

end