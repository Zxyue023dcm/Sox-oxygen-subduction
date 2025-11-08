%% FIGURE 4: Eddy-dominated regions analysis and seasonal cycles
addpath Sox-Oxygen-Subduction\toolbox\m_map1.4.tar\m_map1.4\m_map;
addpath Sox-Oxygen-Subduction\toolbox\Mapcolors

% Load eddy region definitions
mf = matfile('Sox-Oxygen-Subduction\data\results\analysis_eddy_regions.mat');
EddyRegion = mf.EddyRegion;
BASIN = mf.BASIN_ed;

% Extract region data
for r = 1:length(EddyRegion)
    EddyLine(:,:,r) = mf.EddyLine(:,:,r);
    Num{r} = EddyRegion{r};
    OSbar(:,:,r) = mf.OSbar_ed(:,:,r);
    CpnPrp(:,:,r) = mf.CpnPrp_ed(:,:,r);
end

% Load longitude/latitude data
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');

% Load oxygen subduction data for relative percentage calculation
clear OS
m = matfile('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat');

OS(:,:,1) = abs(mean(m.OSeddy, 3, 'omitnan'));
OS(:,:,2) = abs(mean(m.OSeddy, 3, 'omitnan')) + abs(mean(m.OSwb, 3, 'omitnan')) + ...
             abs(mean(m.OSlavoX, 3, 'omitnan')) + abs(mean(m.OSlavoY, 3, 'omitnan'));

clear OSallclim

% Calculate relative percentage of eddy contribution
for i = 1:720
    for j = 1:360
        CpnsProp(i,j) = OS(i,j,1) ./ OS(i,j,2);
    end
end

% Apply Pacific Ocean basin mask
m = matfile('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat');
CpnsProp(m.BASIN ~= 2) = nan;

%% Map showing eddy-dominated regions
h = figure;

% Subplot 1: Relative percentage map
subplot(2,4,1);
m_proj('miller', 'lat', [-60 62], 'long', [122 292]);
hold on
m_pcolor(oceanlon, oceanlat, CpnsProp*100);
shading interp
caxis([0 100]);
colormap(othercolor('Blues9', 10));

c = colorbar('location', 'southoutside', 'Orientation', 'horizontal', ...
    'Ticks', [0 50 100], 'TickLabels', {'0', '50%', '100%'});
c.Label.String = 'Relative Percentage (RP) of Eddy Contribution';
c.FontSize = 12;

m_coast('patch', [.7 .7 .7], 'edgecolor', [.3 .3 .3]);
m_grid('linestyle', 'none', 'tickdir', 'both', 'fontsize', 11, 'gridcolor', [.8 .8 .8], ...
    'linewidth', 1, 'xtick', [135:45:290], 'ytick', [-60 -30 0 30 60]);

title({[' ']; ['  ']; ['  ']}, 'fontsize', 12);

% Add region boundaries to map
% Expand boundaries by 0.5° for better visibility
a = 0.5;
EddyLine([2 3], 2, :) = EddyLine([2 3], 2, :) - a;
EddyLine([1 4 5], 2, :) = EddyLine([1 4 5], 2, :) + a;
EddyLine([1 2 5], 1, :) = EddyLine([1 2 5], 1, :) - a;
EddyLine([3 4], 1, :) = EddyLine([3 4], 1, :) + a;

h;
for n = 1:3
    hold on;
    for i = 1:4
        m_line(EddyLine(i:i+1, 1, n), EddyLine(i:i+1, 2, n), 'color', 'k', 'linewidth', 1);
    end
    m_text(EddyLine(1, 1, n)-8, EddyLine(1, 2, n)-3, Num{n}, 'color', 'k', 'fontsize', 10);
end

%% Plot seasonal cycles of Sox components for each eddy region
cls(1,:) = [.4 .4 .4];        % Total
cls(2,:) = [0.7800 0.6600 0.6300]; % Zonal lateral
cls(3,:) = [0.7600 0.7900 0.4600]; % Meridional lateral
cls(4,:) = [1.0000 0.8300 0.4500]; % Vertical velocity
cls(5,:) = [0.4700 0.8000 0.8800]; % Eddy

mon = [1:12];

h;
for g = 1:size(EddyLine,3)
    % Upper row: Absolute values
    subplot(4,4,g+1);
    hold on; box on;
    
    hb = bar(mon, OSbar(2:5,:,g));
    for c = 1:4
        hb(c).FaceColor = cls(c+1,:);
        hb(c).EdgeColor = cls(c+1,:);
    end
    h1 = plot(mon, OSbar(1,:,g), "Color", cls(1,:), 'LineStyle', '-', 'Marker', '*', 'LineWidth', 1);
    
    title({[Num{g}]});
    set(gca, 'xlim', [0.5 12.5], 'XTick', [1:1:12], 'TickDir', 'both', ...
        'XTickLabel', {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}, ...
        'fontsize', 11);
    
    if g == 1 || g == 4
        ylabel({['S^o^x (Tmol mon^-^1)']});
    end
    
    % Lower row: Relative percentages
    subplot(4,4,g+1+4);
    hold on; box on;
    hbp = bar(mon, CpnPrp(:,:,g), 'stacked');
    for c = 1:4
        hbp(c).FaceColor = cls(c+1,:);
        hbp(c).EdgeColor = cls(c+1,:);
    end
    
    if g == 1 || g == 4
        ylabel({['RP (%)']});
    end
    
    title({[' ']});
    set(gca, 'xlim', [0.5 12.5], 'ylim', [0 100], 'TickDir', 'both', 'XTick', [1:1:12], ...
        'XTickLabel', {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}, ...
        'fontsize', 11);
end

% Add legend
hl = legend([h1 hb], {'Total', 'Zonal La.', 'Merid. La.', 'Vert. Vel.', 'Edd.'}, ...
    'Orientation', 'horizontal');