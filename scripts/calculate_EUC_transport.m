%% Calculate EUC transport (170E-80W, 341:560)
% We define Equatorial Undercurrent (EUC) transport as the eastward velocity 
% between 【±2.5°】 and 【0–350 m】, similar to previous studies.
%
% Citation:
% Stellema, A., Sen Gupta, A., Taschetto, A. S., & Feng, M. (2022). 
% Pacific Equatorial Undercurrent: Mean state, sources, and future changes across models. 
% Frontiers in Climate, 4, 933091. https://doi.org/10.3389/fclim.2022.933091

% Load basin data
load('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat', 'BASIN');

% Load ECCO eastward velocity data (unit: m/s)
load('Sox-Oxygen-Subduction\data\input\velocity_eastward_350m.mat');

% Define calculation region
lonrg = 341:560; % 170E-80W (341:560)
latrg = 176:185; % ±2.5° (176:185)

%% Calculate EUC transport

% Extract velocity data for region
evel = EVEL(lonrg, latrg, :, :);

% Longitude coordinates
long = oceanlon(lonrg, 1);

% Calculate distance coordinates (unit: meters)
[lon_l, lat_l] = distcoordinate(oceanlon(:,1), oceanlat(1,:)');
ans = diff(lon_l, 1, 1);
lon_l = [lon_l(1,:) - ans(1,:); lon_l];
delt_lon = diff(lon_l, 1, 1);
delt_lat = diff(lat_l, 1, 2);
delt_lat = [delt_lat delt_lat(:,1)];
delt_lat = delt_lat(lonrg, latrg);
clear lon_l lat_l delt_lon

% Calculate vertical distances (unit: meters)
delt_z = diff([0; dep], 1, 1);

% Calculate cross-sectional area (m²)
delt_lat_3D = repmat(delt_lat, [1 1 length(dep)]);
Area_3D = delt_lat_3D .* reshape(delt_z, [1 1 length(dep)]);

% Calculate transport (m³/s)
Trspt = squeeze(sum(Area_3D .* evel, [2 3]));

% Convert to Sverdrups (Sv)
Trspt = Trspt / 10^6;

save('Sox-Oxygen-Subduction\data\processed\transport_euc.mat', 'Trspt', 'long');