%% FIGURE 2: Four components of Sox in equatorial region (170E-80W, ±2.5°) across seasons
% Load area data
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');

% Load oxygen subduction components
m = matfile('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat');
% Calculate seasonal means for each component
OS(:,:,:,1) = squeeze(mean(reshape(cat(3, m.OSp(:,:,168), m.OSp(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));
OS(:,:,:,2) = squeeze(mean(reshape(cat(3, m.OSlavoX(:,:,168), m.OSlavoX(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));
OS(:,:,:,3) = squeeze(mean(reshape(cat(3, m.OSlavoY(:,:,168), m.OSlavoY(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));
OS(:,:,:,4) = squeeze(mean(reshape(cat(3, m.OSwb(:,:,168), m.OSwb(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));
OS(:,:,:,5) = squeeze(mean(reshape(cat(3, m.OSeddy(:,:,168), m.OSeddy(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));

% Define equatorial region (170E-80W, ±2.5°)
latrg = 176:185;  % ±2.5° latitude
lonrg = 341:560;  % 170E-80W longitude

% Calculate area-integrated values for equatorial region
os = squeeze(sum(OS(lonrg, latrg, :, :) .* Areaall(lonrg, latrg), [1 2], 'omitnan') ./ 10^6);

% Load EUC transport data
load('Sox-Oxygen-Subduction\data\processes\transport_euc.mat', 'Trspt');
Trspt = mean(reshape(Trspt, [220 12 14]), [1 3], 'omitnan');

% Calculate monthly climatology of total Sox
load('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat','OSp');
OSp = mean(reshape(OSp, [720 360 12 14]), 4, 'omitnan');
os1 = squeeze(sum(OSp(lonrg, latrg, :) .* Areaall(lonrg, latrg), [1 2], 'omitnan') ./ 10^6);

clear OS OSp

%% Subplot a - Seasonal bar plot of Sox components
cls(1,:) = [.4 .4 .4];        % Total
cls(2,:) = [0.7800 0.6600 0.6300]; % Zonal lateral
cls(3,:) = [0.7600 0.7900 0.4600]; % Meridional lateral  
cls(4,:) = [1.0000 0.8300 0.4500]; % Vertical velocity
cls(5,:) = [0.4700 0.8000 0.8800]; % Eddy

h = figure;

subplot(3,1,1);
hold on; box on; grid on;
yyaxis left
b = bar([1:4], os, 'BarWidth', 1);
p = plot([0.5+1/6:1/3:4.5-1/6]', cat(1, os1(12), os1(1:11)), ...
    'color', [.3 .3 .3], 'LineStyle', '-', 'Marker', '*', 'LineWidth', 1);
set(gca, 'YColor', 'k', 'xlim', [0.5 4.5], 'xtick', [1:1:4], ...
    'XTickLabel', {'DJF', 'MAM', 'JJA', 'SON'}, 'ylim', [-40 20], 'fontsize', 11);
xlabel('Season');
ylabel({['Integrated S^o^x'];[' (Tmol mon^-^1)']});

for c = 1:5
    b(c).FaceColor = cls(c,:);
    b(c).EdgeColor = cls(c,:);
end

% Add seasonal separation lines
for x = 1.5:1:3.5
    plot([x x], [-40 20], 'k--');
end

yyaxis right
b(6) = plot([0.5+1/6:1/3:4.5-1/6], cat(2, Trspt(12), Trspt(1:11)), ...
    'color', [0.8000 0.1200 0.5300], 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'YColor', [0.8000 0.1200 0.5300], 'xlim', [0.5 4.5], 'ylim', [0 50], 'fontsize', 11);
ylabel({['EUC Transport'];[' (Sv)']});

legend([b], {'Total', 'Zonal La.', 'Merid. La.', 'Vert. Vel.', 'Edd.'}, ...
    'Orientation', 'horizontal', 'Location', 'north');

%% Clear temporary variables
clear Areaall b c cls latrg lonrg long os OS os1 OSmonclim OSseaclim p Trspt x

%% Subplots bc - Hovmöller diagrams of Sox and vertical velocity
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');

% Load data
load('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat','OSp');
load('Sox-Oxygen-Subduction\data\input\physical_variables.mat', 'wb');

wb = mean(reshape(wb, [720 360 12 14]), 4, 'omitnan');

% Select equatorial region
latrg = 176:185;  % ±2.5° latitude
lonrg = 341:560;  % 170E-80W longitude

OSp = OSp(lonrg, latrg, :);
Areaall = Areaall(lonrg, latrg);
Areaall(isnan(OSp(:,:,1))) = nan;

% Calculate zonal-integrated Sox
OS = squeeze(sum(OSp .* Areaall, 2, 'omitnan')) / 10^6;
OS = mean(reshape(OS, [220 12 14]), 3, 'omitnan');
OSa = OS - mean(OS, 2, 'omitnan');  % Anomalies
OS = repmat(OS, [1 2]);             % Duplicate for continuous cycle
OSa = repmat(OSa, [1 2]);

% Calculate area-weighted vertical velocity
wb = squeeze(sum(wb(lonrg, latrg, :) .* Areaall, 2, 'omitnan') ./ sum(Areaall, 2, 'omitnan'));
wba = wb - mean(wb, 2, 'omitnan');  % Anomalies
wb = repmat(wb, [1 2]);             % Duplicate for continuous cycle
wba = repmat(wba, [1 2]);

% Create coordinate grids
[X, Y] = meshgrid(oceanlon(341:560, 1), 1:24);

% Apply 1° longitude smoothing for visualization
OS = reshape(repmat(mean(reshape(OS, [2 110 24]), 1), [2 1 1]), [220 24]);
OSa = reshape(repmat(mean(reshape(OSa, [2 110 24]), 1), [2 1 1]), [220 24]);
wb = reshape(repmat(mean(reshape(wb, [2 110 24]), 1), [2 1 1]), [220 24]);
wba = reshape(repmat(mean(reshape(wba, [2 110 24]), 1), [2 1 1]), [220 24]);

%% Subplots bc - Plot Hovmöller diagrams
addpath F:\Sox-Oxygen-Subduction\toolbox\Mapcolors

h;
% Subplot 3: Climatological Hovmöller
c1 = subplot(3,2,3);
hold on; box on;
pcolor(X, Y, OS');
shading interp
% Add vertical velocity contours
contour(X, Y, wb', [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wb', [0.5:0.5:3]*12.*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb', -[0.5:0.5:3]*12.*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

caxis([-35 35]*10^(-2));
colormap(c1, flipud(othercolor('RdYlBu10', 30)));
cc1 = colorbar('Ticks', [-0.3 -0.15 0 0.15 0.3]);
cc1.Label.String = 'Tmol mon^-^1';
cc1.FontSize = 11;
xlabel('Longitude');
ylabel('Month');
set(gca, 'fontsize', 11, 'TickDir', 'both', 'YTick', [2:3:24], ...
    'YTickLabel', {'Feb', 'May', 'Aug', 'Nov', 'Feb', 'May', 'Aug', 'Nov'}, ...
    'ylim', [1 24], 'xlim', [170 280], 'XTick', [180:30:280], ...
    'XTickLabel', {'180°', '150°W', '120°W', '90°W'});

% Subplot 4: Anomaly Hovmöller
c2 = subplot(3,2,4);
hold on; box on;
pcolor(X, Y, OSa');
shading interp
% Add vertical velocity anomaly contours
contour(X, Y, wba', [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wba', [0.25:0.25:3]*12.*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wba', -[0.25:0.25:4]*12.*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

% Add propagation lines
plot([269.25 174.75], [3 9], 'k', 'LineStyle', '-', 'Marker', 'none');
plot([267.25 174.75], [8 15], 'k', 'LineStyle', '-', 'Marker', 'none');

caxis([-25 25]*10^(-2));
colormap(c2, polarmap(30));
cc1 = colorbar('Ticks', [-0.2 -0.1 0 0.1 0.2]);
cc1.Label.String = 'Tmol mon^-^1';
cc1.FontSize = 11;
xlabel('Longitude');
ylabel('Month');
set(gca, 'fontsize', 11, 'TickDir', 'both', 'YTick', [2:3:24], ...
    'YTickLabel', {'Feb', 'May', 'Aug', 'Nov', 'Feb', 'May', 'Aug', 'Nov'}, ...
    'ylim', [1 24], 'xlim', [170 280], 'XTick', [180:30:280], ...
    'XTickLabel', {'180°', '150°W', '120°W', '90°W'});

%% Clear temporary variables
clear Areaall c1 c2 cc1 latrg lonrg oceanlat oceanlon OS OSa OSlavoX_new OSp wb wba X Y yrs

%% Subplot d - ENSO analysis based on Nino3.4 index (2004-2017)
% Identify El Niño and La Niña events using 5-month running mean threshold of ±0.5°C

% Nino3.4 index data from: https://psl.noaa.gov/data/timeseries/monthly/NINO34/
% Period: 2004-2017 with 2 additional months at both ends
nino34 = [0.33 0.43 0.26 0.17 -0.10 0.06 0.10 0.14 0.41 0.66 0.67 0.73 0.62 0.71...
     0.56 0.26 0.28 0.28 0.30 0.22 -0.01 -0.04 -0.08 -0.15 -0.44 -0.75...
-0.98 -0.71 -0.73 -0.30 -0.11 0.09 0.03 0.37 0.63 0.76 0.98 1.10...
   0.59 0.12 -0.15 -0.16 -0.39 -0.16 -0.37 -0.57 -1.04 -1.40 -1.58 -1.61...
   -1.79 -1.70 -1.17 -0.89 -0.64 -0.44 -0.04 -0.04 -0.28 -0.30 -0.37 -0.90...
   -1.00 -0.71 -0.72 -0.25 0.17 0.49 0.69 0.62 0.68 0.96 1.49 1.81...
    1.43 1.18 1.07 0.56 -0.15 -0.62 -0.89 -1.33 -1.56 -1.65 -1.57 -1.63...
   -1.70 -1.26 -0.98 -0.74 -0.53 -0.25 -0.23 -0.66 -0.76 -0.93 -1.09 -1.05...
   -0.93 -0.61 -0.48 -0.29 -0.18 0.14 0.44 0.66 0.44 0.23 0.33 -0.13...
   -0.42 -0.40 -0.14 -0.08 -0.28 -0.33 -0.28 -0.29 -0.09 -0.24 -0.02 -0.09...
  -0.42 -0.45 -0.07 0.28 0.45 0.48 0.13 0.14 0.37 0.48 0.89 0.77...
    0.59 0.57 0.48 0.90 1.04 1.28 1.56 1.87 2.01 2.21 2.57 2.56...
   2.56 2.11 1.60 1.05 0.45 0.06 -0.25 -0.48 -0.46 -0.75 -0.63 -0.51...
   -0.34 -0.01 -0.09 0.22 0.30 0.22 0.22 -0.18 -0.56 -0.52 -0.84 -0.85 -0.98 -0.78];

%% Identify ENSO events using 5-month running threshold
ifenso = zeros(172,1); % 1 - El Niño; -1 - La Niña; 0 - Normal
for i = 3:length(nino34)-2
    data = [nino34(i-2) nino34(i-1) nino34(i) nino34(i+1) nino34(i+2)];
    
    if sum(data > 0.5) == 5
        ifenso(i-2:i+2) = 1;
    elseif sum(data < -0.5) == 5
        ifenso(i-2:i+2) = -1;
    else
        ifenso(i) = 0;
    end
end

ifenso = ifenso(3:170,:);
nino34 = nino34(:,3:170);

%% Calculate seasonally detrended Sox anomalies

% Load basin data for equatorial Pacific
load('F:\Sox-Oxygen-Subduction\data\PO_Basins.mat', 'BASIN');
load('F:\Sox-Oxygen-Subduction\data\Areaall.mat');
load('F:\Sox-Oxygen-Subduction\data\LonLat.mat');
load('F:\Sox-Oxygen-Subduction\data\OSresults1_monthly.mat', 'OSp');

% Calculate integrated Sox for equatorial Pacific (2.5S-2.5N, 170E-80W)
for i = 1:168
    os = OSp(:,:,i);
    area = Areaall;
    
    os(BASIN ~= 2) = nan;
    os(oceanlat > 2.5 | oceanlat < -2.5 | oceanlon < 170 | oceanlon > 280) = nan;
    area(isnan(os)) = nan;
    os(isnan(area)) = nan;
    
    fullseries(i,1) = sum(os .* area, 'all', 'omitnan') / 10^6;
end

% Remove seasonal cycle
repeatseries = repmat(mean(reshape(fullseries, [12 14]), 2), [14 1]);
deSox = fullseries - repeatseries;

%% Design low-pass filter for ENSO signal

cutoff_period = 24; % Cutoff period (months): corresponds to ENSO timescale (2-7 years)
fs_monthly = 12;    % Sampling frequency (samples/year)
cutoff_frequency = 1 / (cutoff_period / 12); % Cutoff frequency (year⁻¹)
% Calculate normalized cutoff frequency
Wn = cutoff_frequency / (fs_monthly / 2);
fprintf('Actual filter parameters:\n');
fprintf('Wn = %.6f\n', Wn);

% Design Butterworth low-pass filter
order = 4;
[b, a] = butter(order, Wn, 'low');

% Apply filter to detrended Sox
x = deSox;
y = filtfilt(b, a, x);

% Calculate correlation with Nino3.4 index
[a, b] = corrcoef(y, nino34);
R2 = a(1,2);

%% Plot ENSO analysis
h;
subplot(3,1,3);
year = [2004:1/12:2018-1/12];
hold on; box on;
yyaxis left
plot([2004 2018], [0 0], 'k--');

% Shade ENSO events
for j = 1:length(ifenso)
    v2 = [year(j)-1/24 -30; year(j)-1/24 30; year(j)+1/24 30; year(j)+1/24 -30];
    f2 = [1 2 3 4];
    
    if ifenso(j) == 1
        patch('Faces', f2, 'Vertices', v2, 'FaceColor', 'r', 'FaceAlpha', .3, 'EdgeColor', 'none');
    elseif ifenso(j) == -1
        patch('Faces', f2, 'Vertices', v2, 'FaceColor', 'b', 'FaceAlpha', .3, 'EdgeColor', 'none');
    end
end

p1(1) = plot([2004:1/12:2017+11/12], y, 'k-', 'LineWidth', 1);
ylabel({['Seasonally Detrended S^o^x'];['(Tmol mon^-^1)']});
set(gca, 'YColor', 'k', 'ylim', [-10 10], 'xlim', [2004 2018], 'fontsize', 9);

yyaxis right
p1(2) = plot(2004:1/12:2018-1/12, nino34, 'r-', 'LineWidth', 0.8);
ylabel(['Nino3.4 Index']);
set(gca, 'YColor', 'r', 'ylim', [-3 3], 'xlim', [2004 2018], 'fontsize', 9);

title({['r = ' num2str(R2, '%.2f') ' , p < 0.01']});