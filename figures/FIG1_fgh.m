
%% FIGURE 1 continued: [Total OS] DJF - JJA seasonal analysis
% Load longitude/latitude data
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');

% Load basin definitions
load('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat','BASIN1','BASIN','groupnum');
groupnum = flipud(groupnum);

% Load area data
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');

% Load oxygen subduction data
load('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat','OSp');
OSp(repmat(BASIN,[1 1 168]) ~= 2) = nan;

% Calculate seasonal climatology (DJF, MAM, JJA, SON)
OS = squeeze(mean(reshape(cat(3, OSp(:,:,168), OSp(:,:,1:167)), [720 360 3 4 14]), [3 5], 'omitnan'));

%% Calculate latitudinal distributions for DJF-JJA comparison
% Metrics: [Subduction area fraction, Total subduction rate, Subduction rate per unit area]

% Define latitude bins
latstr{2,1} = '[-89:2:89]';
latstr{2,2} = [720*4 90 1];
eval(['lats = ' latstr{2,1} ';']);

% Calculate latitudinal distributions for each season
for s = 1:4
    os = OS(:,:,s);
    os(isnan(BASIN1)) = nan;
    Area = Areaall;
    
    oslat = reshape(os, latstr{2,2});
    area = reshape(Area, latstr{2,2});
    
    for j = 1:size(lats,2)
        % Area fraction calculation
        line1 = oslat(:,j); 
        line2 = area(:,j);
        line2(isnan(line1)) = nan;
        
        % Subduction area fraction
        subarea = sum(line2(line1 > 0), 'omitnan') ./ sum(line2, 'omitnan');
        obarea = sum(line2(line1 < 0), 'omitnan') ./ sum(line2, 'omitnan');
        OSArea_lat(j,s) = subarea; 
        
        % Total subduction amount
        subsum = sum(line2 .* line1, 'omitnan') ./ 10^6; 
        OSSum_lat(j,s) = subsum; % Tmol mon^-1
        
        % Subduction rate per unit area
        subs = sum(line2(line1 > 0) .* line1(line1 > 0), 'omitnan') ./ sum(line2(line1 > 0), 'omitnan');
        OSV_S_lat(j,s) = subs;  % mol m^-2 mon^-1
        subo = sum(line2(line1 < 0) .* line1(line1 < 0), 'omitnan') ./ sum(line2(line1 < 0), 'omitnan');
        OSV_O_lat(j,s) = subo;  % mol m^-2 mon^-1
    end
end

%% Calculate regional statistics for 5 Pacific basins

for g = 1:5
    for s = 1:4
        os = OS(:,:,s);
        Area = Areaall;
        
        % Apply basin mask
        os(BASIN1 ~= groupnum{g,1}) = nan;
        Area(BASIN1 ~= groupnum{g,1}) = nan;
        os(isnan(Area)) = nan;
        Area(isnan(os)) = nan;
        
        % Area fraction
        subarea = sum(Area(os > 0), 'all', 'omitnan') ./ sum(Area, 'all', 'omitnan');
        obarea = sum(Area(os < 0), 'all', 'omitnan') ./ sum(Area, 'all', 'omitnan');
        OSArea(g,s) = subarea; 
        
        % Total subduction and obduction amounts
        subsum = sum(Area(os > 0) .* os(os > 0), 'all', 'omitnan') ./ 10^6;
        OSSum_S(g,s) = subsum; % Tmol mon^-1
        obsum = sum(Area(os < 0) .* os(os < 0), 'all', 'omitnan') ./ 10^6;
        OSSum_O(g,s) = obsum; % Tmol mon^-1
        OSSum(g,s) = sum(Area .* os, 'all', 'omitnan') ./ 10^6;
        
        % Subduction rate per unit area
        subs = sum(Area(os > 0) .* os(os > 0), 'all', 'omitnan') ./ sum(Area(os > 0), 'all', 'omitnan');
        OSV_S(g,s) = subs;  % mol m^-2 mon^-1
        subo = sum(Area(os < 0) .* os(os < 0), 'all', 'omitnan') ./ sum(Area(os < 0), 'all', 'omitnan');
        OSV_O(g,s) = subo;
        OSV(g,s) = sum(Area .* os, 'all', 'omitnan') ./ sum(Area, 'all', 'omitnan');
    end
end

%% Create horizontal bar plots for 5 regions [DJF-JJA comparison]

% Define colors for DJF and JJA
cls(1,:) = [0.16 0.62 0.56]; % DJF color
cls(2,:) = [0.91 0.44 0.32]; % JJA color

% Font size settings
sz1 = 11;

h = figure;

% ========================== Fraction of Subduction Area ======================
subplot(3,1,1);
hold on; box on;
yyaxis left
p = plot(lats, OSArea_lat(:,[1 3])*100, 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'xlim', [-62 62], 'ylim', [0 120], 'fontsize', sz1, ...
    'xtick', [-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel', {'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'}, ...
    'YColor', 'k', 'ytick', [0 25 50 75 100], 'fontsize', sz1);
ylabel({['Fraction of'];['Subduction Area'];['(%)']}, 'fontsize', sz1);

yyaxis right
bh = bar([-47.5 -22.5 0 25 50], OSArea(:,[1 3])*100);
for c = 1:2
   bh(c).FaceColor = cls(c,:);
   bh(c).EdgeColor = cls(c,:); 
   bh(c).FaceAlpha = 0.85;
   p(c).Color = [cls(c,:) 1];
end
set(gca, 'xlim', [-62 62], 'ylim', [0 120], 'YColor', 'k', 'fontsize', sz1, ...
    'ytick', [0 25 50 75 100]);

% ========================== S^o^x and O^o^x rate ======================
subplot(3,1,2);
hold on; box on;
yyaxis left
plot([-65 65], [0 0], 'k-');
p = plot(lats, OSV_S_lat(:,[1 3]), 'LineStyle', '-', 'Marker', 'none', 'LineWidth', 1);
p1 = plot(lats, OSV_O_lat(:,[1 3]), 'LineStyle', '--', 'Marker', 'none', 'LineWidth', 1);

for c = 1:2
    p1(c).Color = [cls(c,:) 1];
    p(c).Color = [cls(c,:) 1];
end
ax1 = gca;
set(ax1, 'xlim', [-62 62], 'ylim', [-6.5 6.5], 'fontsize', sz1, ...
   'xtick', [-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel', {'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'}, ...
    'YColor', 'k');
ylabel({['S^o^x and O^o^x rate'];['(mol m^-^2 mon^-^1)']}, 'fontsize', sz1);

yyaxis right
bh = bar([-47.5 -22.5 0 25 50], OSV_S(:,[1 3]));
b1 = bar([-47.5 -22.5 0 25 50], OSV_O(:,[1 3]));
for c = 1:2
    bh(c).FaceColor = cls(c,:);
    bh(c).EdgeColor = cls(c,:); 
    bh(c).FaceAlpha = 0.85;
    b1(c).FaceColor = 'none';
    b1(c).EdgeColor = cls(c,:); 
end
set(gca, 'xlim', [-62 62], 'ylim', [-6.5 6.5], 'fontsize', sz1, 'YColor', 'k');

% ========================== Zonally integrated S^o^x ======================
subplot(3,1,3);
hold on; box on;
yyaxis left
plot([-65 65], [0 0], 'k-');
p = plot(lats, OSSum_lat(:,[1 3]), 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'xlim', [-62 62], 'ylim', [-16 16], 'fontsize', sz1, ...
   'xtick', [-60 -47.5 -35 -22.5 -10 0 10 25 40 50 60], ...
    'XTickLabel', {'60°S',' ','35°S',' ','10°S','0','10°N',' ','40°N',' ','60°N'}, ...
    'YColor', 'k');
xlabel({['Latitude']}, 'fontsize', sz1);
ylabel({['Zonally integrated S^o^x'];['(Tmol mon^-^1)']}, 'fontsize', sz1);

yyaxis right
bh = bar([-47.5 -22.5 0 25 50], OSSum(:,[1 3]));
for c = 1:2
    bh(c).FaceColor = cls(c,:);
    bh(c).EdgeColor = cls(c,:); 
    bh(c).FaceAlpha = 0.85;
    p(c).Color = [cls(c,:) 1];
end
ylabel({['Sub-regional S^o^x'];['(Tmol mon^-^1)']}, 'fontsize', sz1);
set(gca, 'xlim', [-62 62], 'ylim', [-25 25], 'fontsize', sz1, 'YColor', 'k');

%% Add regional boundaries and labels
h;
for i = 1:3
    subplot(3,1,i);
    hold on;
    % Add vertical lines marking regional boundaries
    plot([-35 -35], [-50 120], 'k--');
    plot([-10 -10], [-50 120], 'k--');
    plot([10 10], [-50 120], 'k--');
    plot([40 40], [-50 120], 'k--');
end

% Add region names to first subplot
subplot(3,1,1);
hold on;
for g = 1:5
    text(-50, 100, groupnum{g,5}, "FontSize", 12, "FontName", 'Times New Roman');
end

% Add legend
legend(bh, {'DJF', 'JJA'}, 'fontsize', 12, 'Orientation', 'horizontal'); 