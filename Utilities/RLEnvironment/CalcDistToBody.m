function [bodyDist] = CalcDistToBody(s, varargin)
%[bodyDist] = CalcHPDirect(s, varargin)
%   Calculates distance to body using the fast marching method

grid_size       = s.wrld.size;
target_voxels   = [s.clc.startSR ; s.clc.startSC ; s.clc.startSZ]';

% Initialize a large value for all grid cells
distance = inf * ones(grid_size);

% Set distance of target voxels to 0
for i = 1:size(target_voxels, 1)
    distance(target_voxels(i, 1), target_voxels(i, 2), target_voxels(i, 3)) = 0;
end

% Create a list of "active" voxels (initially, the target voxels)
active_voxels = target_voxels;

% Loop until there are no more active voxels
while ~isempty(active_voxels)
    new_active_voxels = [];

    % Iterate over active voxels
    for i = 1:size(active_voxels, 1)
        voxel = active_voxels(i, :);

        % Iterate over neighboring voxels
        for dx = -1:1
            for dy = -1:1
                for dz = -1:1
                    neighbor = voxel + [dx, dy, dz];

                    % Check if the neighbor is inside the grid
                    if all(neighbor > 0) && all(neighbor <= grid_size)

                        % Calculate the distance from the current voxel to its neighbor
                        d = sqrt(dx^2 + dy^2 + dz^2);
                        new_distance = distance(voxel(1), voxel(2), voxel(3)) + d;

                        % Update the distance if it's smaller than the current value
                        if new_distance < distance(neighbor(1), neighbor(2), neighbor(3))
                            distance(neighbor(1), neighbor(2), neighbor(3)) = new_distance;

                            % Add the neighbor to the new list of active voxels
                            new_active_voxels = [new_active_voxels; neighbor];
                        end
                    end
                end
            end
        end
    end

    % Update the list of active voxels
    active_voxels = unique(new_active_voxels, 'rows');
end

bodyDist = distance;

end

