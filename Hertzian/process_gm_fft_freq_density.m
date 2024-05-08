function process_gm_fft_freq_density(time_vector, index_particles, index_oscillating_wall, driving_frequency, driving_amplitude, position_particles, initial_distance_from_oscillation)
    % Purpose - finds attenuation, wavenumber, and wave speed for a ganular mechanics simulation.
    %
    % Format:   [fitted_attenuation, wavenumber, wavespeed] = ...
    %            process_gm_fft(plot_title, time_vector, index_particles, index_oscillating_wall, driving_frequency, driving_amplitude, position_particles)
    %
    % Input:    time_vector            - your time-series for that have the same length of position series
    %           index_particles        - index for each particle to be analyzed
    %           index_oscillating_wall - the index for particles that make up the oscillating wall
    %           driving_frequency      - frequency of the oscillating wall in units of inverse time
    %           driving_amplitude      - amplitude of wall in units of distance
    %           position_particles     - position-series for particles for whatever axis you want to analyze
    %
    % Output:   Plot showing the density of frequency strength as a function of distance
    %
    % Note:     This is very computationally expensive. Change iskip, frequency_cutoff, and position_cutoff as needed

iskip = 1;
frequency_cutoff = 2;  % Define the frequency cutoff
position_cutoff = 6; % Distance cutoff

% Pre Allocate for Speed
average_dt = mean(diff(time_vector));
sampling_freq = 1 / average_dt;
nyquist_freq = sampling_freq / 2;  % Nyquist frequency 
freq_vector = linspace(0, 1, fix(length(time_vector)/2)+1) * nyquist_freq;
index_vector = 1:numel(freq_vector);

% Filter to only include frequencies below the cutoff
freq_vector = freq_vector(freq_vector <= frequency_cutoff);
index_vector = index_vector(freq_vector <= frequency_cutoff);

% Initialize the plot
figure;
hold on;
cmap = jet(64);  % Colormap from blue to red
amp_min = inf;
amp_max = 0;

for nn = index_particles(1:iskip:end)  % Incremental index processing
    if ~index_oscillating_wall(nn)  % Ensure the particle is not on the oscillating wall
        position_nn = position_particles(nn, :);  % Extract position time-series for particle nn
        
        if length(unique(position_nn)) > 10  % Process only if significant movement
            if initial_distance_from_oscillation(nn) < position_cutoff
                % Center and normalize the data
                centered_data = position_nn - mean(position_nn);
                normalized_fft_data = fft(centered_data) / length(time_vector);
                normalized_fft_data_single_sided_nn = abs(normalized_fft_data(index_vector)) * 2; % Doubled to change change from double sided (about f=0) to single sided
                
                % Update amplitude range for color mapping
                amp_min = min(amp_min, min(normalized_fft_data_single_sided_nn));
                amp_max = max(amp_max, max(normalized_fft_data_single_sided_nn));
                
                % Plot each frequency component
                distance_from_oscillation = initial_distance_from_oscillation(nn);  % Initial position of the particle
                for mM = 1:length(freq_vector)
                    % fprintf('Working on freq_vector %d', mM) % Testing purposes to make sure my computer isn't bricked

                    % Normalize amplitude to [1, 64] for colormap indexing. ceil() needed because MATLAB only lines integer values
                    amp_index = max(1, ceil(64 * (normalized_fft_data_single_sided_nn(mM) - amp_min) / (amp_max - amp_min))); % Shift amplitude to base at zero and normalize. Then scale by 64 to match color scale.
                    plot(distance_from_oscillation, freq_vector(mM), 'o', 'Color', cmap(amp_index, :), 'MarkerFaceColor', cmap(amp_index, :));
                end
            end
        end
    end
end

% Highlight the driving frequency
% line(xlim, [driving_frequency driving_frequency], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
% text(mean(xlim), driving_frequency, sprintf(' Driving Frequency: %.2f Hz', driving_frequency), ...
%      'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'black');

xlabel('Initial Position from Oscillating Wall');
ylabel('Frequency (Inverse Time)');
title('Frequency Spectrum of Particles');
colorbar;  % Shows the color scale
colormap(cmap);  % Ensures the colorbar uses the same colormap
caxis([amp_min, amp_max]);  % Set the colorbar's amplitude range
% xlim([0 max(initial_distance_from_oscillation)]);
ylim([0,max(freq_vector)]);
grid on;
hold off;
