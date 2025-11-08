
function [CURL,gXy,gYx] = curldffs(datax,distx,datay,disty,method)
%CURLDFFS    Calculate curl using central difference method
%   [CURL,gXy,gYx] = curldffs(datax,distx,datay,disty,method)
%
%   Input grid data, variables are 2D matrices
%   Differentiation method: gradient 
%   Reduces two lengths in the first dimension
%   Note: For global calculations, matrix must be organized with 
%         first dimension as longitude, second dimension as latitude
%
%   Inputs:
%       datax   - x-component data matrix
%       distx   - x-direction distance coordinates  
%       datay   - y-component data matrix
%       disty   - y-direction distance coordinates
%       method  - differentiation method (currently supports 'gradient')
%
%   Outputs:
%       CURL    - calculated curl field
%       gXy     - derivative of x-component in y-direction
%       gYx     - derivative of y-component in x-direction

lon_l = distx;
lat_l = disty;

if method == 'gradient'
    % Extend longitude boundaries using forward/backward differences
    ans = diff(lon_l,1,1);
    lon_l = [lon_l(1,:)-ans(1,:); lon_l; lon_l(end,:)+ans(1,:)];
    
    % Extend latitude boundaries (replicate first and last rows)
    lat_l = [lat_l(1,:); lat_l; lat_l(end,:)];
    
    % Extend data boundaries (periodic boundary conditions)
    datax = [datax(end,:); datax; datax(1,:)];
    datay = [datay(end,:); datay; datay(1,:)];

    % Calculate gradients for each latitude point
    for p = 1:size(lat_l,2)
        hx = lon_l(:,p); 
        hy = lat_l(1,:); 
        
        % Compute gradients using MATLAB's gradient function
        [py,~] = gradient(datax, hy, hx); 
        [~,qx] = gradient(datay, hy, hx); 
        
        % Store gradient components
        gXy1(:,p) = py(:,p);
        gYx1(:,p) = qx(:,p);
    end
    
    % Remove extended boundaries from gradient results
    gXy = gXy1(2:size(lat_l,1)-1,:);
    gYx = gYx1(2:size(lat_l,1)-1,:);
end

% Calculate curl as the difference between derivatives
CURL = gYx - gXy;

end



