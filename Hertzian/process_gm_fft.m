function [fitted_attenuation, wavenumber, wavespeed] = process_gm_fft(time_vector, index_particles, index_oscillating_wall, driving_frequency, position_particles, figure_handle)
% Purpose - finds attenuation, wavenumber, and wave speed for a ganular mechanics simulation.
%
% Format:   [fitted_attenuation, wavenumber, wavespeed] = ...
%            process_gm_fft(plot_title, time_vector, index_particles, index_oscillating_wall, driving_frequency, driving_amplitude, position_particles)
%
% Input:    time_vector            - your time-series for that have the same length of position series
%           index_particles        - index for each particle to be analyzed
%           index_oscillating_wall - the index for particles that make up the oscillating wall
%           driving_frequency      - frequency of the oscillating wall in units of inverse time
%           figure_handle          - which axsis you are looking at, for figure saving purposes
%           position_particles     - position-series for particles for whatever axis you want to analyze
%
% Output:   fitted_attenuation      - attenuation that was fit from the amplitude-distance relationship
%           wavenumber              - wavenumber fit from the phase-distance relationship
%           wavespeed              - wavespeed from wavespeed = frequency / wavenumber
%
% Note:     

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];

iskip = 1;
freq_match_tolerance = 1;


% Pre Allocate for Speed
average_dt = mean(diff(time_vector));
sampling_freq = 1 / average_dt;
nyquist_freq = sampling_freq / 2;  % Nyquist frequency 
freq_vector = linspace(0, 1, fix(length(time_vector)/2)+1) * nyquist_freq;
index_vector = 1:numel(freq_vector);
number_elements_time = numel(time_vector);

for nn = index_particles(1:iskip:end) % Sorts them by increments of iskip...for iskip>1, speeds things up
    if(~index_oscillating_wall(nn)) % Executes the following block only if the particle indexed by nn is not on the left wall.
        position_nn = position_particles(nn,:); % extracts and stores the time_vector-series positions of the particle indexed by nn from the array x_all AKA looks at and picks out 
       
        if length(unique(position_nn))>10 % Checks if position_nn has more than 100 unique values AKA only process data that moves
            centered_data = position_nn-mean(position_nn); %Center the data on zero for mean
            normalized_fft_data = fft(centered_data)/number_elements_time; 

            % Find the dominant frequency and its max amplitude
            %   Need to double it because when signal is centered, power is
            %   distributed in both positive and negative. double the abs accounts for this
            [max_particle_amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
            dominant_frequency = freq_vector(idx_max);

            % Find the index of the closest frequency to the desired frequency
            [~, idx_driving_freq] = min(abs(freq_vector - driving_frequency));

            % Check if there is a peak around the desired frequency and amplitude is greater than 
            if idx_driving_freq > 1 && idx_driving_freq < numel(freq_vector) 
                %  fprintf('*** Checking for Slope  ***\n');
                % Calculate the sign of the slope before and after the desired frequency
                sign_slope_before = sign(normalized_fft_data(idx_driving_freq) - normalized_fft_data(idx_driving_freq - 1));
                sign_slope_after = sign(normalized_fft_data(idx_driving_freq + 1) - normalized_fft_data(idx_driving_freq));
                
                % Check if the signs of the slopes are different and if the values on both sides are greater than the value at the desired frequency
                if sign_slope_before ~= sign_slope_after && abs(normalized_fft_data(idx_driving_freq - 1)) < abs(normalized_fft_data(idx_driving_freq)) && abs(normalized_fft_data(idx_driving_freq + 1)) < abs(normalized_fft_data(idx_driving_freq)) && abs(dominant_frequency - driving_frequency) < freq_match_tolerance
                    %  fprintf('Peak found around the driving frequency. Storing data\n');

                    amplitude_vector = [amplitude_vector, max_particle_amplitude]; % Pulls amplitude from fft calculation
                    initial_position_vector = [initial_position_vector, position_nn(1)];
                    phase_vector = [phase_vector, angle(normalized_fft_data(idx_driving_freq))];
                    cleaned_particle_index = [cleaned_particle_index, nn];
                else
                    %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');
                end
            else
                %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');
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
xlabel('Distance');
ylabel('Particle Oscillation Amplitude');
legend('show');
grid on; 
figure_name = [figure_handle, '_attenuation_plot.fig']
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
xlabel('Initial Position');
ylabel('\Delta\phi');

% Customizing y-axis to show multiples of pi
y_max = max(unwrapped_phase_vector);  % Get the maximum y value
y_min = min(unwrapped_phase_vector);  % Get the minimum y value
yticks = [ceil(y_min/pi)*pi:pi:floor(y_max/pi)*pi];  % Define y-ticks in steps of pi
yticklabels = arrayfun(@(x) sprintf('%.2f\\pi', x/pi), yticks, 'UniformOutput', false);  % Create custom y-tick labels
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);  % Apply custom ticks and labels

% Hold off to finish the plotting
hold off;
