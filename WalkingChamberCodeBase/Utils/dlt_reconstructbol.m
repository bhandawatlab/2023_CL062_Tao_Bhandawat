function [xyz,rmse] = dlt_reconstructbol(camPts,c)

% function [xyz,rmse] = dlt_reconstruct(camPts,c)
%
% This function reconstructs the 3D position of a coordinate based on a set
% of DLT coefficients and [u,v] pixel coordinates from 2 or more cameras
%
% Inputs:
%  c - 11 DLT coefficients for all n cameras, [11,n] array
%  camPts - [u,v] pixel coordinates from all n cameras over f frames,
%   [f,2*n] array
%
% Outputs:
%  xyz - the xyz location in each frame, an [f,3] array
%  rmse - the root mean square error for each xyz point, and [f,1] array,
%   units are [u,v] i.e. camera coordinates or pixels
%
% Ty Hedrick
%
% This function is an implementaion of DLT reconstruction as described at http://www.kwon3d.com/theory/dlt/dlt.html
% However it does not use the weight matrix in the iteration.

% number of frames; this is the number of input pixels (u, v) we want to map into object space (x, y, z)
nFrames = size(camPts, 1);

% number of cameras
nCams = round(size(camPts, 2) / 2);

% setup output variables
xyz(1:nFrames, 1:3) = NaN;
rmse(1:nFrames, 1) = NaN;

L = c; % Makes the following equations easier to understand IMO because it follows the source doc

% process each frame; for each input pixel we must calculate a different 
for i = 1:nFrames
    
    % get a list of cameras with non-NaN [u,v]
    m = find(isnan(camPts(i, 1:2:nCams * 2)) == false);
    numM = numel(m);
    
    % If we don't have at least 2 good cameras, skip this point and continue with the next frame
    if numM < 2
        continue
    end
    
    % We need to iterate to solve the reconstruction problem because R is a function of x, y, z. See Eq. 23
    counter = 1;
    
    while true
        % initialize least-square solution matrices; there are 2 equations per camera. See Eq. 22
        X = zeros(numM * 2, 3);
        Y = zeros(numM * 2, 1);

        % Calculate R
        % Because R is a function of x, y, z, you have to iterate over each individual pixel multiple times in order to 
        % get the actual (converged) x, y, z position for that given pixel. 
        R = ones(1, nCams);
        %if counter > 1
        %    R(m) = L(9, m) .* xyz(i, 1) + L(10, m) .* xyz(i, 2) + L(11, m) .* xyz(i, 3) + 1;
        %end

        X(1:2:numM * 2, 1) = camPts(i, m * 2 - 1) .* L(9, m) - L(1, m) ./ R(m);
        X(1:2:numM * 2, 2) = camPts(i, m * 2 - 1) .* L(10, m) - L(2, m) ./ R(m);
        X(1:2:numM * 2, 3) = camPts(i, m * 2 - 1) .* L(11, m) - L(3, m) ./ R(m);

        X(2:2:numM * 2, 1) = camPts(i, m * 2) .* L(9, m) - L(5, m) ./ R(m); 
        X(2:2:numM * 2, 2) = camPts(i, m * 2) .* L(10, m) - L(6, m) ./ R(m);
        X(2:2:numM * 2, 3) = camPts(i, m * 2) .* L(11, m) - L(7, m) ./ R(m);


        Y(1:2:numM * 2) = L(4, m) - camPts(i, m * 2 - 1) ./ R(m);
        Y(2:2:numM * 2) = L(8, m) - camPts(i, m * 2) ./ R(m);

        % get the least squares solution to the reconstruction
        xyz(i, 1:3) = linsolve(X, Y);
        
        % Stop iterating when the x, y, z values have converged or 1000 iterations have passed
        if counter > 1000 || counter > 1 && sqrt(sum((savedXYZ - xyz(i, 1:3)) .^ 2)) / 3 < 0.1
            break
        end
        savedXYZ = xyz(i, 1:3);
        counter = counter + 1;
    end
    
    % Display the number of iterations. It will have a number of iterations for each frame i.e. each input pixel
    % fprintf('There were %d iterations.\n', counter - 1);
    
    % For a more informative error, you should calculate it between the output x, y, z points and some measured x, y, z
    % points, not between the matrix Y and X * (the output x, y, z points). That will only tell you how much error is in
    % the matrix solver (linsolve) component of this process. It is not the overall reconstruction error.
    % See the script test_reconstruction.m for details on estimating the reconstruction error.

    % Compute the estimated [u,v] given the solution to the reconstruction
    uv = X * xyz(i, 1:3)';

    % compute the number of degrees of freedom in the reconstruction
    dof = numel(Y) - 3;

    % estimate the root mean square reconstruction error
    rmse(i, 1) = (sum((Y - uv) .^ 2) / dof) ^ 0.5;
end