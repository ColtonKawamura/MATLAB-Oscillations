function [fitted_attenuation, wavenumber, wavespeed] = process_gm_fft(time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, figure_handle, initial_distance_from_oscillation)
% Purpose - finds attenuation, wavenumber, and wave speed for a ganular mechanics simulation.
%
% Format:   [fitted_attenuation, wavenumber, wavespeed] = ...
%            process_gm_fft(plot_title, time_vector, index_particles, index_oscillating_wall, driving_frequency, driving_amplitude, position_particles)
%
% Input:    time_vector                          - your time-series for that have the same length of position series
%           index_particles                      - index for each particle to be analyzed
%           index_oscillating_wall               - the index for particles that make up the oscillating wall
%           driving_frequency                    - frequency of the oscillating wall in units of inverse time
%           figure_handle                        - which axis you are looking at, for figure saving purposes
%           position_particles                   - position-series for particles for whatever axis you want to analyze
%           initial_distance_from_oscillation    - distance to the particle on the axis normal to oscillating wall (if oscillating in x, this is each particles initial x-position)
%
% Output:   fitted_attenuation       - attenuation that was fit from the amplitude-distance relationship
%           wavenumber               - wavenumber fit from the phase-distance relationship
%           wavespeed                - wavespeed from wavespeed = frequency / wavenumber
%           {}_attenuation_plot.fig  - plot of the attenuation where {} is whatever you put for your figure_handle input
%
% Note:     Dial in your criteria for clean data with below parameters

% Define the threshold frequency and flag
threshold_frequency = .05*driving_frequency; % disregards frequencies less than 10 percent of driving frequency
ignore_below_threshold = true; % Set to false to disable filtering

% Define minimum peak amplitude / prominence to count as "good"
min_peak_amplitude = 5E-5;


%-------------------------------------------------------------------------------------------------------------------------
% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];

iskip = 1;
freq_match_tolerance = 0.05;


% Pre Allocate for Speed
average_dt = mean(diff(time_vector));
sampling_freq = 1 / average_dt;
nyquist_freq = sampling_freq / 2;  % Nyquist frequency 
freq_vector = linspace(0, 1, fix(length(time_vector)/2)+1) * nyquist_freq;
index_vector = 1:numel(freq_vector);
number_elements_time = numel(time_vector);

for nn = index_particles(1:iskip:end)
    if(~index_oscillating_wall(nn))

        position_nn = position_particles(nn,:); 
       
        if length(unique(position_nn))>10

            centered_data = position_nn-mean(position_nn); %Center the data on zero for mean
            normalized_fft_data = fft(centered_data)/number_elements_time; 

            % If you want to ignore frequencies below a certain threshold, only take those at or above the threshold
            if ignore_below_threshold

                valid_indices = freq_vector >= threshold_frequency; % Shave off the indices that are below the threshold
                [peak_amplitudes, peak_locations] = findpeaks(abs(normalized_fft_data(valid_indices) * 2), 'MinPeakProminence', min_peak_amplitude); % Find the peaks within some Minimum Prominence
                % Now need to shift found peaks back to their original index before threshold
                peak_locations = peak_locations + find(valid_indices, 1, 'first') - 1; % find(..) tells you where your filtered data starts relative to the original full dataset. -1 because Matlab starts at 1 index
            else
                [peak_amplitudes, peak_locations] = findpeaks(abs(normalized_fft_data * 2), 'MinPeakProminence', 0); % Adjust 'MinPeakProminence' based on noise
            end

            if ~isempty(peak_amplitudes)
                [max_particle_amplitude, idx_max] = max(peak_amplitudes);
                dominant_frequency = freq_vector(peak_locations(idx_max));

                % Find the index of the closest frequency to the desired frequency
                [~, idx_driving_freq] = min(abs(freq_vector - driving_frequency));

                % makes sure driving frequency is within the range of spectrum
                if idx_driving_freq > 1 && idx_driving_freq < numel(freq_vector) 

                    % Only passes data that has a domininat frequency close to the driving frequency
                    if abs(dominant_frequency - driving_frequency) < freq_match_tolerance
                        amplitude_vector = [amplitude_vector, max_particle_amplitude]; 
                        initial_position_vector = [initial_position_vector, initial_distance_from_oscillation(nn)];
                        phase_vector = [phase_vector, angle(normalized_fft_data(idx_driving_freq))];
                        cleaned_particle_index = [cleaned_particle_index, nn];
                    end
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attenuation Fitting and Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform linear fit
coefficients = polyfit(initial_position_vector, log(abs(amplitude_vector)), 1);

% Extract fitted_attenuation and intercept
fitted_attenuation = coefficients(1);
intercept_attenuation = coefficients(2);

% Create a linear fit line
fit_line = exp(intercept_attenuation) * exp(initial_position_vector.*fitted_attenuation);

% Plot original data and linear fit
handle_figure = figure;
semilogy(initial_position_vector, abs(amplitude_vector), 'bo', 'DisplayName', 'Data');
hold on;
semilogy(initial_position_vector, fit_line, 'r-', 'DisplayName', 'Linear Fit');
xlabel('Distance from Oscillation (Particle Diameters)');
ylabel('Particle Oscillation Amplitude');
legend('show');
grid on; 
figure_name = [figure_handle, '_attenuation_plot.fig'];
savefig(handle_figure, figure_name);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wavenumber and Speed Fitting and Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unwrapped_phase_vector = unwrap(phase_vector);

% Fit a line to the data
p = polyfit(initial_position_vector, unwrapped_phase_vector, 1);
fitted_line = polyval(p, initial_position_vector);

% Store the slope of the line as wavenumber
wavenumber = p(1);
wavespeed = driving_frequency/wavenumber;

% Plot initial position vs. phase as dots
figure;
scatter(initial_position_vector, unwrapped_phase_vector, 'o');
grid on;
hold on;  % Keep the plot for adding the fitted line

% Plot the fitted line
plot(initial_position_vector, fitted_line, '-r');

% Label the axes
xlabel('Distance from Oscillation (Particle Diameters)');
ylabel('\Delta\phi');

% Customizing y-axis to show multiples of pi
y_max = max(unwrapped_phase_vector);  % Get the maximum y value
y_min = min(unwrapped_phase_vector);  % Get the minimum y value
yticks = [ceil(y_min/pi)*pi:pi:floor(y_max/pi)*pi];  % Define y-ticks in steps of pi
yticklabels = arrayfun(@(x) sprintf('%.2f\\pi', x/pi), yticks, 'UniformOutput', false);  % Create custom y-tick labels
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);  % Apply custom ticks and labels

% Hold off to finish the plotting
hold off;
