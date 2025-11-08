
function [DIV,gX,gY] = divdffs(datax,distx,datay,disty,method)
%DIVDFFS Calculate divergence using backward/central difference methods
%   [DIV,gX,gY] = divdffs(datax,distx,datay,disty,method)
%
%   Input grid data as 2D matrices
%   Backward/forward methods: reduce by one length in both dimensions
%   Gradient method: reduce by two lengths in first dimension
%   Note: For global calculations, matrix must be organized with 
%         first dimension as longitude, second dimension as latitude
%
%   Inputs:
%       datax   - x-component data matrix
%       distx   - x-direction distance coordinates
%       datay   - y-component data matrix  
%       disty   - y-direction distance coordinates
%       method  - differentiation method ('backward' or 'gradient')
%
%   Outputs:
%       DIV     - calculated divergence field
%       gX      - x-component derivative
%       gY      - y-component derivative

lon_l = distx;
lat_l = disty;

if method == 'backward'
    % Backward difference method
    
    % Calculate x-direction differences for data
    for i = 2:size(datax,1)
        delX(i-1,:) = (datax(i,:) - datax(i-1,:));
    end
    
    % Calculate x-direction differences for coordinates
    for i = 2:size(lon_l,1)
        delx(i-1,:) = (lon_l(i,:) - lon_l(i-1,:));
    end
    
    % Compute x-gradient
    gX = delX(:,2:end) ./ delx(:,2:end);
    
    % Calculate y-direction differences for data
    for i = 2:size(datay,2)
        delY(:,i-1) = (datay(:,i) - datay(:,i-1));
    end
    
    % Calculate y-direction differences for coordinates
    for i = 2:size(lat_l,2)
        dely(:,i-1) = (lat_l(:,i) - lat_l(:,i-1));
    end
    
    % Compute y-gradient
    gY = delY(2:end,:) ./ dely(2:end,:);
    
    % Adjust coordinate matrices to match gradient dimensions
    lon_l = lon_l(2:end,2:end);
    lat_l = lat_l(2:end,2:end);
    
    % Smooth x-gradient by averaging adjacent longitude bands
    for a = 1:size(gX,1)/12
        gX(2*a-1:2*a,:) = repmat(mean(gX(2*a-1:2*a,:), 1, 'omitnan'), [2,1]);
    end
    for a = 1:size(gX,1)/12-1
        gX(2*a:2*a+1,:) = repmat(mean(gX(2*a:2*a+1,:), 1, 'omitnan'), [2,1]);
    end
    for a = 1:size(gX,1)/12
        gX(2*a-1:2*a,:) = repmat(mean(gX(2*a-1:2*a,:), 1, 'omitnan'), [2,1]);
    end
    
    % Smooth y-gradient by averaging adjacent latitude bands
    for b = 1:size(gY,2)/12
        gY(:,2*b-1:2*b) = repmat(mean(gY(:,2*b-1:2*b), 2, 'omitnan'), [1,2]);
    end
    for b = 1:size(gY,2)/12-1
        gY(:,2*b:2*b+1) = repmat(mean(gY(:,2*b:2*b+1), 2, 'omitnan'), [1,2]);
    end
    for b = 1:size(gY,2)/12
        gY(:,2*b-1:2*b) = repmat(mean(gY(:,2*b-1:2*b), 2, 'omitnan'), [1,2]);
    end
    
elseif method == 'gradient'
    % Gradient method with boundary extension
    for p = 1:size(lat_l,2)
        hx = lon_l(:,p); 
        hy = lat_l(1,:); 
        
        % Compute gradients using MATLAB's gradient function
        [~,px] = gradient(datax, hy, hx); 
        [qy,~] = gradient(datay, hy, hx); 
        
        % Store gradient components
        guHmlx1(:,p) = px(:,p);
        guHmly1(:,p) = qy(:,p);
    end
    
    % Remove extended boundaries from gradient results
    gX = guHmlx1(2:size(lat_l,1)-1,:);
    gY = guHmly1(2:size(lat_l,1)-1,:);
end

% Calculate divergence as sum of derivatives
DIV = gX + gY;

end

