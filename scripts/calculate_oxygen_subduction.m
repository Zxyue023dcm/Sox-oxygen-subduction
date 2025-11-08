%% Load data

% Oxygen data
load('Sox-Oxygen-Subduction\data\input\oxygen_concentration.mat');
% Physical data
load('Sox-Oxygen-Subduction\data\input\physical_variables.mat');
% Area data
load('Sox-Oxygen-Subduction\data\grid\grid_area_pacific.mat');
% Longitude/latitude data
load('Sox-Oxygen-Subduction\data\grid\grid_coordinates.mat');

%% Term 1 & 2 - Lateral induction in eastward(X) and northward(Y) directions

% Calculate distance coordinates
[lon_l, lat_l] = distcoordinate(oceanlon(:,1), oceanlat(1,:)');
ans = diff(lon_l, 1, 1);
lon_l = [lon_l(1,:) - ans(1,:); lon_l];
lat_l = [lat_l(1,:); lat_l];
lon_l = [nan(size(lon_l,1), 1) lon_l];
lat_l = [nan(size(lat_l,1), 1) lat_l];

% Calculate water mass subduction rate - lateral induction
for i = 1:168
    % Prepare mixed layer depth with boundary extension
    hml = [Hml(end,:,i); Hml(:,:,i)];
    hml = [nan(size(hml,1), 1) hml];
    [~, gHmlx, gHmly] = divdffs(hml, lon_l, hml, lat_l, 'backward');

    % Prepare stratification depth with boundary extension  
    hml = [Hst(end,:,i); Hst(:,:,i)];
    hml = [nan(size(hml,1), 1) hml];
    [~, gHstx, gHsty] = divdffs(hml, lon_l, hml, lat_l, 'backward');

    % Calculate lateral induction components
    Slavo_x(:,:,i) = gHmlx1 .* evel1(:,:,i) + gHstx1 .* evel2(:,:,i);
    Slavo_y(:,:,i) = gHmly1 .* nvel1(:,:,i) + gHsty1 .* nvel2(:,:,i);
end

% Calculate oxygen subduction rate Sox = S*O2 using upwind scheme
OSlavoX = Slavo_x;
% Positive velocities (use oxygen concentration at 10m below mixed layer base)
OSlavoX(Slavo_x > 0) = -Slavo_x(Slavo_x > 0) .* oxy10(Slavo_x > 0) .* ...
                        rhob10(Slavo_x > 0) * 10^(-6) * 3600 * 24 * 365 / 12;
% Negative velocities (use oxygen concentration at mixed layer base)
OSlavoX(Slavo_x < 0) = -Slavo_x(Slavo_x < 0) .* oxy(Slavo_x < 0) .* ...
                        rho(Slavo_x < 0) * 10^(-6) * 3600 * 24 * 365 / 12;

OSlavoY = Slavo_y;
% Positive velocities (use oxygen concentration at 10m below mixed layer base)
OSlavoY(Slavo_y > 0) = -Slavo_y(Slavo_y > 0) .* oxy10(Slavo_y > 0) .* ...
                        rhob10(Slavo_y > 0) * 10^(-6) * 3600 * 24 * 365 / 12;
% Negative velocities (use oxygen concentration at mixed layer base)
OSlavoY(Slavo_y < 0) = -Slavo_y(Slavo_y < 0) .* oxy(Slavo_y < 0) .* ...
                        rho(Slavo_y < 0) * 10^(-6) * 3600 * 24 * 365 / 12;

OSlavo = OSLavoX + OSLavoY;
save('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat', ...
     'OSlavoX', 'OSlavoY', 'OSlavo', '-v7.3');

%% Term 3 - Vertical velocity contribution

OSwb = wb;
% Positive vertical velocities (use oxygen concentration at 10m below mixed layer base)
OSwb(wb > 0) = -wb(wb > 0) .* oxy10(wb > 0) .* rhob10(wb > 0) * ...
                10^(-6) * 3600 * 24 * 365 / 12;
% Negative vertical velocities (use oxygen concentration at mixed layer base)
OSwb(wb < 0) = -wb(wb < 0) .* oxy(wb < 0) .* rho(wb < 0) * ...
                10^(-6) * 3600 * 24 * 365 / 12;

save('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat', 'OSwb', '-append');

%% Term 4 - Eddy-induced contribution

% Calculate water mass subduction rate - eddy component
for i = 1:length(year)
    % Calculate eddy transport components
    estarh = estar1(:,:,i) .* Hml(:,:,i) + estar2(:,:,i) .* Hst(:,:,i);
    nstarh = nstar1(:,:,i) .* Hml(:,:,i) + nstar2(:,:,i) .* Hst(:,:,i);

    % Apply boundary conditions
    es = [estarh(end,:); estarh; estarh(1,:)];
    ns = [nstarh(end,:); nstarh; nstarh(1,:)];
    
    % Calculate eddy-induced subduction
    [Seddy(:,:,i), ~, ~] = divdffs(es, lon_l, ns, lat_l, 'gradient');
end

% Calculate oxygen subduction rate for eddy component
OSeddy = Seddy;
% Positive eddy subduction (use oxygen concentration at 10m below mixed layer base)
OSeddy(Seddy > 0) = -Seddy(Seddy > 0) .* oxy10(Seddy > 0) .* ...
                     rhob10(Seddy > 0) * 10^(-6) * 3600 * 24 * 365 / 12;
% Negative eddy subduction (use oxygen concentration at mixed layer base)
OSeddy(Seddy < 0) = -Seddy(Seddy < 0) .* oxy(Seddy < 0) .* ...
                     rho(Seddy < 0) * 10^(-6) * 3600 * 24 * 365 / 12;

save('Sox-Oxygen-Subduction\data\results\oxygen_subduction_results.mat', 'OSeddy', '-append');