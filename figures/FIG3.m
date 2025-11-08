%% FIGURE 3: Zonal mean analysis of Sox, vertical velocity, and wind stress curl
addpath Sox-Oxygen-Subduction\toolbox\Mapcolors

load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');
load('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat','BASIN');

% Load wind stress curl data
load('Sox-Oxygen-Subduction\data\processed\wind_stress_curl_zonal_climatology.mat','WCm');
WCm = repmat(WCm, [1 2]);          % Duplicate for continuous annual cycle
WC_a = WCm - mean(WCm, 2, 'omitnan'); % Anomalies

%% Calculate zonal mean Sox and vertical velocity
% Load area data
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');
Areaall(BASIN ~= 2) = nan;
Areaall = repmat(Areaall, [1 1 12]);

% Load and process Sox data
load('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat','OSp');
Sox = mean(reshape(OSp, [720 360 12 14]), 4, 'omitnan');
Basin = repmat(BASIN, [1 1 12]);
Sox(Basin ~= 2) = nan;
clear OSp
Areaall(isnan(Sox)) = nan;
Sox(isnan(Areaall)) = nan;

% Calculate zonal mean Sox
latstr = '[-89:2:89]';
eval(['lats = ' latstr ';']);
sox = reshape(Sox, [720*4 90 12]);
Areaall = reshape(Areaall, [720*4 90 12]);
sox = squeeze(sum(sox .* Areaall, 1, 'omitnan') ./ 10^6);
sox(sox == 0) = nan;
sox = repmat(sox, [1 2]);          % Duplicate for continuous annual cycle
sox_anom = sox - mean(sox, 2, 'omitnan'); % Anomalies

% Load and process vertical velocity data
load('Sox-Oxygen-Subduction\data\input\physical_variables.mat', 'wb');
wb = mean(reshape(wb, [720 360 12 14]), 4, 'omitnan');
Basin = repmat(BASIN, [1 1 12]);
wb(Basin ~= 2) = nan;
wb(isnan(Sox)) = nan;

% Calculate zonal mean vertical velocity
wb = reshape(wb, [720*4 90 12]);
wb = squeeze(mean(wb, 1, 'omitnan'));
wb(wb == 0) = nan;
wb = repmat(wb, [1 2]);            % Duplicate for continuous annual cycle
wb_anom = wb - mean(wb, 2, 'omitnan'); % Anomalies

% Create coordinate grid
[X, Y] = meshgrid([1:24], lats);

%% Plot Sox + vertical velocity contours
h = figure;

% Subplot 1: Sox climatology with vertical velocity contours
c1 = subplot(1,4,1);
hold on; box on;
pcolor(X, Y, sox);
shading interp
% Add vertical velocity contours
contour(X, Y, wb, [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wb, [1:1:3].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb, [3:2:10].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb, -[1:1:3].*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

ylabel('Latitude');
set(gca, 'fontsize', 12, 'ylim', [-65 60], 'xlim', [1 24], 'TickDir', 'both', ...
    'YTick', [-60:5:60], 'YTickLabel', {'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'}, ...
    'XTick', [2:3:24], 'XTickLabel', {'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

caxis([-15 15]);
colormap(c1, flipud(othercolor('RdYlBu10', 18)));
cc1 = colorbar('location', 'southoutside', 'Orientation', 'horizontal', ...
    'Ticks', [-15:5:15], 'Limits', [-15 15]);
cc1.Label.String = 'S^o^x (Tmol mon^-^1)';
grid on; box on;
cc1.FontSize = 12;

%% Subplot 2: Sox anomalies with vertical velocity anomaly contours
c2 = subplot(1,4,2);
hold on; box on;
pcolor(X, Y, sox_anom);
shading interp
% Add vertical velocity anomaly contours
contour(X, Y, wb_anom, [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wb_anom, [0.5:0.5:3].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb_anom, -[0.5:0.5:4].*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

set(gca, 'fontsize', 12, 'ylim', [-65 60], 'xlim', [1 24], 'TickDir', 'both', ...
    'YTick', [-60:5:60], 'YTickLabel', {'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'}, ...
    'XTick', [2:3:24], 'XTickLabel', {'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

caxis([-7 7]);
colormap(c2, polarmap(28));
cc2 = colorbar('location', 'southoutside', 'Orientation', 'horizontal', ...
    'Ticks', [-6:2:6], 'Limits', [-7 7]);
cc2.Label.String = 'S^o^x Anomalies (Tmol mon^-^1)';
cc2.FontSize = 12;

%% Subplot 3: Wind stress curl with vertical velocity contours
c3 = subplot(1,4,3);
hold on; box on;
pcolor(X, Y, WCm);
shading interp
% Add vertical velocity contours
contour(X, Y, wb, [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wb, [1:1:3].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb, [3:2:10].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb, -[1:1:3].*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

set(gca, 'fontsize', 12, 'ylim', [-65 60], 'xlim', [1 24], 'TickDir', 'both', ...
    'YTick', [-60:5:60], 'YTickLabel', {'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'}, ...
    'XTick', [2:3:24], 'XTickLabel', {'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

caxis([-3 3]*10^(-7));
colormap(c3, flipud(othercolor('BrBG5', 14)));
cc3 = colorbar('location', 'southoutside', 'Orientation', 'horizontal', ...
    'Ticks', [-3:1.5:3]*10^(-7), 'Limits', [-3 3]*10^(-7));
cc3.Label.String = 'WSC (N^2 m^-^3)';
cc3.FontSize = 12;

%% Subplot 4: Wind stress curl anomalies with vertical velocity anomaly contours
c4 = subplot(1,4,4);
hold on;
pcolor(X, Y, WC_a);
shading interp
% Add vertical velocity anomaly contours
contour(X, Y, wb_anom, [0 0], 'color', [.7 .7 .7], 'Linewidth', 1, 'ShowText', 'off');
contour(X, Y, wb_anom, [0.5:0.5:4].*10^(-6), 'LineStyle', '-', 'color', [.7 .7 .7]);
contour(X, Y, wb_anom, -[0.5:0.5:4].*10^(-6), 'LineStyle', '--', 'color', [.7 .7 .7]);

set(gca, 'fontsize', 12, 'ylim', [-65 60], 'xlim', [1 24], 'TickDir', 'both', ...
    'YTick', [-60:5:60], 'YTickLabel', {'60°S','','','45°S','','','30°S','','','15°S','','','0','','','15°N','','','30°N','','','45°N','','','60°N'}, ...
    'XTick', [2:3:24], 'XTickLabel', {'Feb','May','Aug','Nov','Feb','May','Aug','Nov'});

caxis([-1.5 1.5]*10^(-7));
colormap(c4, polarmap(28));
cc4 = colorbar('location', 'southoutside', 'Orientation', 'horizontal', ...
    'Ticks', [-1.5:0.5:1.5]*10^(-7), 'Limits', [-1.5 1.5]*10^(-7));
cc4.Label.String = {['WSC Anomalies'];[' (N^2 m^-^3)']};
cc4.Label.FontSize = 12;
grid on; box on;

%% Add seasonal region boxes to all subplots
h;
% Red boxes - subtropical winter regions
for c = 1:4
    subplot(1,4,c);
    hold on;
    % Southern Hemisphere: May-November, lat: -11° to -33°
    line([5 5], [-11 -33], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([5 11], [-11 -11], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([5 11], [-33 -33], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([11 11], [-11 -33], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    
    % Northern Hemisphere: November-April, lat: 11° to 35°
    line([11 11], [11 35], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([11 16], [11 11], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([11 16], [35 35], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
    line([16 16], [11 35], 'color', 'r', 'linestyle', '-', 'linewidth', 1);
end

% Blue boxes - subpolar summer regions
for c = 1:4
    subplot(1,4,c);
    hold on;
    % Northern Hemisphere: May-November, lat: 35° to 60°
    line([5 5], [35 60], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([5 11], [35 35], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([5 11], [60 60], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([11 11], [35 60], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    
    % Southern Hemisphere: November-April, lat: -33° to -49°
    line([11 11], [-33 -49], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([11 16], [-33 -33], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([11 16], [-49 -49], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
    line([16 16], [-33 -49], 'color', 'b', 'linestyle', '-', 'linewidth', 1);
end