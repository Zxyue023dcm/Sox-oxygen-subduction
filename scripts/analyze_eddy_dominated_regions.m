%% Analysis of Eddy-Dominated Regions
% Load area data
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');
% Load longitude/latitude data
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');

% Load oxygen subduction data
clear OS
m = matfile('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat');
OS(:,:,1) = abs(mean(m.OSeddy,3,'omitnan'));
OS(:,:,2) = abs(mean(m.OSeddy,3,'omitnan')) + abs(mean(m.OSwb,3,'omitnan')) + ...
             abs(mean(m.OSlavoX,3,'omitnan')) + abs(mean(m.OSlavoY,3,'omitnan'));

%% Select three example regions | Mark BASIN
EddyRegion = {'Kuroshi Extension', 'North Equatorial Pacific', 'Coastal Peru'};

% Define polygon boundaries for each eddy region
% Kuroshi Extension
EddyLine(:,:,1) = [144 45; 144 33; 170 33; 170 45; 144 45];
% North Equatorial Pacific  
EddyLine(:,:,2) = [155 14; 155 5.5; 270 5.5; 270 14; 155 14];
% Coastal Peru
EddyLine(:,:,3) = [260 -6; 260 -32; 290 -32; 290 -6; 260 -6];

%% Mark only regions where eddy contribution RP > 50%

% Calculate Relative Percentage (RP)
for i = 1:720
    for j = 1:360
        CpnsProp(i,j) = OS(i,j,1) ./ OS(i,j,2);
    end
end

% Apply Pacific Ocean basin mask
mb = matfile('Sox-Oxygen-Subduction\data\grid\basin_mask_pacific.mat');
CpnsProp(mb.BASIN ~= 2) = nan;

% Mark grid points with strong eddy contribution
BASIN = nan(720,360);
for n = 1:3
    BASIN(oceanlon >= EddyLine(1,1,n) & oceanlon <= EddyLine(3,1,n) & ...
          oceanlat >= EddyLine(2,2,n) & oceanlat <= EddyLine(1,2,n)) = n;
end

BASIN_ed = BASIN;
for n = 1:3
    BASIN_ed(BASIN == n & CpnsProp < 0.5) = nan;
end

save('Sox-Oxygen-Subduction\data\results\analysis_eddy_regions.mat', ...
     'BASIN_ed', 'EddyRegion', 'EddyLine', '-v7.3');

%% Seasonal time variation of Sox components in eddy regions [Net subduction, RP]

% Load oxygen subduction components
clear OS
OS(:,:,:,1) = m.OSp;
OS(:,:,:,2) = m.OSlavoX; 
OS(:,:,:,3) = m.OSlavoY;
OS(:,:,:,4) = m.OSwb;
OS(:,:,:,5) = m.OSeddy;
OS = reshape(OS, [720 360 12 14 5]);

% Calculate net oxygen subduction for each term and region
for i = 1:5
    for g = 1:size(EddyRegion,2)
        bnum = g;
        for y = 1:14
            for s = 1:12
                os = OS(:,:,s,y,i);
                data = [];
                for b = 1:size(bnum,1)
                    data = [data; os(BASIN_ed == bnum(b)) ...
                            Areaall(BASIN_ed == bnum(b)) ...
                            BASIN_ed(BASIN_ed == bnum(b))];
                end
                line1 = data(:,1); 
                line2 = data(:,2);
                line1(isnan(line2)) = nan; 
                line2(isnan(line1)) = nan;

                % Calculate subduction and obduction
                subsum = sum(line2(line1 > 0) .* line1(line1 > 0), 'omitnan') ./ 10^6;
                obsum = sum(line2(line1 < 0) .* line1(line1 < 0), 'omitnan') ./ 10^6;
                % Net oxygen subduction (T mol /month)
                POterms_NET(i,s,g,y) = (subsum + obsum) ./ 12;
            end
        end
    end
end

% Calculate multi-year monthly mean
OSbar_ed = mean(POterms_NET, 4, 'omitnan');

% Calculate relative percentage (RP) of magnitude for each term
CpnPrp_ed = abs(OSbar_ed(2:5,:,:)) ./ sum(abs(OSbar_ed(2:5,:,:)), 1, 'omitnan') * 100;

save('Sox-Oxygen-Subduction\data\results\analysis_eddy_regions.mat', ...
     'OSbar_ed', 'CpnPrp_ed', '-append');

