%% Calculate Wind Stress Curl (WSC)

% Load longitude/latitude data
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');

% Load eastward and northward wind stress components
load('Sox-Oxygen-Subduction\data\input\physical_variables.mat', 'taue', 'taun');

%% Wind stress curl calculation [Δtaun/Δx - Δtaue/Δy] [unit: N^2/m^3]
[lon_l, lat_l] = distcoordinate(oceanlon(:,1), oceanlat(1,:)');

for m = 1:168
    % Calculate wind stress curl using gradient method
    [~, gXy, gYx] = curldffs(taue(:,:,m), lon_l, taun(:,:,m), lat_l, 'gradient');
    
    WC(:,:,m) = gYx1 - gXy1;
end

save('Sox-Oxygen-Subduction\data\processed\wind_stress_curl.mat', 'WC');

%% Calculate zonal integral of 12-month climatology

% Calculate monthly climatology
WCm = mean(reshape(WC, [720 360 12 14]), 4, 'omitnan');
% Apply Pacific Ocean basin mask
WCm(repmat(BASIN, [1 1 12]) ~= 2) = nan;

% Define latitude bins
latstr = '[-89:2:89]';
eval(['lats = ' latstr ';']);

% Reshape and calculate zonal mean
WCm = reshape(WCm, [720*4 90 12]);
WCm = squeeze(mean(WCm, 1, 'omitnan'));
WCm(WCm == 0) = nan;

save('Sox-Oxygen-Subduction\data\processed\wind_stress_curl_zonal_climatology.mat', 'WCm');