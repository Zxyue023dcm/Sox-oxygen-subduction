%% FIGURE 1 (6 subplots): [climatological Sox, STD, Rseason] * [spatial distribution + latitudinal distribution]
addpath Sox-Oxygen-Subduction\toolbox\m_map1.4.tar\m_map1.4\m_map;
addpath Sox-Oxygen-Subduction\toolbox\Mapcolors

load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');
load('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat','BASIN');
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');

% Load net oxygen subduction data
load('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat','OSp');

% Apply Pacific Ocean basin mask
OSp(repmat(BASIN,[1 1 168]) ~= 2) = nan;

% Calculate climatological mean Sox
OSc = mean(OSp, 3, 'omitnan');

%% Spatial distribution of climatological Sox

h = figure;

% Subplot 1: Climatological mean
c2 = subplot(3,1,1);
m_proj('miller','lat',[-63 62],'long',[122 294]);
hold on
m_pcolor(oceanlon, oceanlat, OSc); 
shading interp
caxis([-7 7]); 
colormap(c2, flipud(othercolor('RdYlBu10',14)));
c = colorbar('location','westoutside');
c.Label.String = {['O^o^x                                   S^o^x'];['(mol m^-^2 mon^-^1)']};
c.FontSize = 11;

m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',10);

% Add latitude reference lines
m_line([122 292],[-35 -35],'linewidth',1,'linestyle','--','color','k');
m_line([122 292],[-10 -10],'linewidth',1,'linestyle','--','color','k');
m_line([122 292],[10 10],'linewidth',1,'linestyle','--','color','k');
m_line([122 292],[40 40],'linewidth',1,'linestyle','--','color','k');

%% Spatial distribution of seasonal standard deviation

% Calculate monthly climatology
OSstd = mean(reshape(OSp,[720 360 12 14]), 4, 'omitnan');
stdp = std(OSstd, 0, 3, 'omitnan');

% Subplot 2: Seasonal variability
h;
c4 = subplot(3,1,2);
m_proj('miller','lat',[-63 62],'long',[122 294]);
hold on
m_pcolor(oceanlon, oceanlat, stdp);
colormap(c4, parula(12));
caxis([0 3]);
c = colorbar('location','westoutside');
c.Label.String = {['σ_s_e_a_s_o_n'];['(mol m^-^2 mon^-^1)']};
c.FontSize = 11;

m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',10);

%% Spatial distribution of seasonal correlation

% Calculate interannual variability time series
OSia = OSp; % monthly data

% Create repeated seasonal cycle time series
OSsc = repmat(mean(reshape(OSia,[720 360 12 14]), 4, 'omitnan'), [1 1 14]);

% Calculate correlation coefficients
for i = 1:720
    for j = 1:360
        x = squeeze(OSia(i,j,:));
        y = squeeze(OSsc(i,j,:));
        [a, b] = corrcoef(x', y');  
        R(i,j) = a(1,2);
        pValue(i,j) = b(1,2);
    end
end

% Apply statistical significance filtering
R(BASIN ~= 2) = nan;
R(pValue > 0.05) = nan; % Remove statistically insignificant points (p > 0.05)

% Identify non-significant points for scatter plot
graydot = [oceanlon(pValue > 0.05 & BASIN == 2) oceanlat(pValue > 0.05 & BASIN == 2)];

%% Plot seasonal correlation
h;
c6 = subplot(3,1,3);
m_proj('miller','lat',[-63 62],'long',[122 294]);
hold on
% Plot non-significant points in gray
m_scatter(graydot(:,1), graydot(:,2), 5, [.2 .2 .2], 'filled');
% Plot R² values for significant correlations
m_pcolor(oceanlon, oceanlat, R.^2);
c = colorbar('location','westoutside');
c.Label.String = 'r^2_s_e_a_s_o_n';
c.FontSize = 11;
caxis([0 1]);
colormap(c6, othercolor('PuBuGn7'));

m_coast('patch',[.7 .7 .7],'edgecolor',[.7 .7 .7]);
m_grid('linestyle','none','tickdir',[],'gridcolor',[.8 .8 .8],...
     'linewidth',1,'xtick',[135 180:45:300],'ytick',[-60 -35 -10 10 40 60],'fontsize',10);

%% Add panel labels and region annotations
h;
Num = {'(a)','(b)','(c)'};
for n = 1:3
    subplot(3,1,n);
    hold on; box on
    text(-1.5, 1.1, Num{n}, 'FontSize', 12, 'FontWeight', 'bold');
end

% Add region labels to first subplot
h;
subplot(3,1,1);
hold on; 
text(1.5, 1.1, 'SSP', 'FontSize', 12, 'fontName', 'Times New Roman');
text(1.5, 1.1, 'NSTP', 'FontSize', 12, 'fontName', 'Times New Roman');
text(1.5, 1.1, 'EP', 'FontSize', 12, 'fontName', 'Times New Roman');
text(1.5, 1.1, 'SSTP', 'FontSize', 12, 'fontName', 'Times New Roman');
text(1.5, 1.1, 'SAP', 'FontSize', 12, 'fontName', 'Times New Roman');