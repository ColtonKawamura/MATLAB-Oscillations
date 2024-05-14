function [r, g_r] = compute_rdf_with_neighbors(x_all, y_all, neighbor_list_all, time_vector, r_max, dr)
% Purpose - Computes the radial distribution function (RDF) for a system of particles using a neighbor list.
%
% Format:   [r, g_r] = compute_rdf_with_neighbors(x_all, y_all, neighbor_list_all, r_max, dr)
%
% Input:    x_all              - NxM array of x positions (N particles, M time steps)
%           y_all              - NxM array of y positions (N particles, M time steps)
%           neighbor_list_all  - Cell array where each cell contains the indices of neighboring particles for each particle
%           r_max              - Maximum distance to consider for the RDF
%           dr                 - Bin width for distances, makes the curve more preceise 
%
% Output:   r                  - Array of distance bins
%           g_r                - RDF values corresponding to each distance bin
%

% Number of particles and time steps
num_particles = length(neighbor_list_all);
num_timesteps = length(time_vector);

% Area of the circle with radius r_max (for normalization in 2D)
area_max = pi * r_max^2;
particle_density = num_particles / area_max; % Number density

% Define distance bins
r = 0:dr:r_max; % determines the radial sections that r_max is divided into
g_r = zeros(size(r)); % initalizes g_r and it's sections
    
% Compute pair distances using neighbor list for each time step
for t = 1:num_timesteps
    for nn = 1:num_particles
        particle_pos = [x_all(nn, t), y_all(nn, t)];
        neighbors = neighbor_list_all{nn};
        
        for j = neighbors
            if j > nn % To avoid double-counting pairs, works since we are starting from nn = 1
                neighbor_pos = [x_all(j, t), y_all(j, t)];
                distance_to_neighbor = norm(particle_pos - neighbor_pos);
                if distance_to_neighbor < r_max
                    bin_index = floor(distance_to_neighbor / dr) + 1; % Need to +1 because MATLAB starts at index = 1, so index=0 needs to shift to 1.
                    g_r(bin_index) = g_r(bin_index) + 2; % Each pair contributes to the count 
                end
            end
        end
    end
end

% Normalize RDF
for k = 1:length(r)-1 % goes through each r-bin
    r_inner = r(k);
    r_outer = r(k+1);
    shell_area = pi * (r_outer^2 - r_inner^2);
    g_r(k) = g_r(k) / (num_particles * num_timesteps * particle_density * shell_area); % g(k)/num_particles = average overall particles, 
end

% Plot RDF
figure;
plot(r, g_r, 'LineWidth', 2);
xlabel('Distance (Particle Diameters)');
ylabel('g(r)');
title('Radial Distribution Function');
grid on;