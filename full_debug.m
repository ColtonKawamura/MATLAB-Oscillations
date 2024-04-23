[~,idx]=min(abs(x0-44))

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
valid_probe_numbers = [];
initial_phase_offset = 0;
tolerance = 0.5

for nn = isort(1:iskip:end) # Sorts them by increments of iskip...for iskip>1, speeds things up
    if(~left_wall_list(nn)) % Executes the following block only if the particle indexed by nn is not on the left wall.
        x_temp = x_all(nn,:); % Extracts and stores the time-series positions of the particle indexed by nn from the array x_all

        if length(unique(x_temp)) > 10 % Checks if x_temp has more than 100 unique values AKA only process data that moves
            % fprintf('*** Doing fft fit  ***\n');
            probe_data = x_temp;
            average_dt = mean(diff(time));
            sampling_freq = 1 / average_dt;
            Fn = sampling_freq / 2; % Nyquist frequency

            % Determine the index to start from the last third of the data
            start_index = floor(2 * numel(probe_data) / 3) + 1;

            % Select only the last third of the probe_data
            last_third_data = probe_data(start_index:end);

            % Adjust the number of elements to reflect the last third of data
            third_number_elements_time = numel(last_third_data);

            % Center the data on zero for mean
            centered_data = last_third_data - mean(last_third_data);

            % Perform FFT normalized by the number of elements
            normalized_fft_data = fft(centered_data) / third_number_elements_time;

            % Frequency vector calculation adjusted for the new data length
            freq_vector = linspace(0, 1, fix(third_number_elements_time / 2) + 1) * Fn;

            % Index vector for the frequencies
            index_vector = 1:numel(freq_vector);

            % Find the dominant frequency and its amplitude
            [amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
            dominant_frequency = freq_vector(idx_max);

            % Find the index of the frequency closest to driving frequency
            desired_frequency = driving_frequency;

            % Find the index of the closest frequency to the desired frequency
            [~, idx_desired] = min(abs(freq_vector - desired_frequency));

            % Check if there is a peak around the desired frequency and amplitude is greater than 
            if idx_desired > 1 && idx_desired < numel(freq_vector) && amplitude > cutoff_amplitude && abs(dominant_frequency - desired_frequency) < tolerance
                % fprintf('*** Checking for Slope  ***\n');
                % Calculate the sign of the slope before and after the desired frequency
                sign_slope_before = sign(normalized_fft_data(idx_desired) - normalized_fft_data(idx_desired - 1));
                sign_slope_after = sign(normalized_fft_data(idx_desired + 1) - normalized_fft_data(idx_desired));
                
                % Check if the signs of the slopes are different and if the values on both sides are greater than the value at the desired frequency
                if sign_slope_before ~= sign_slope_after && normalized_fft_data(idx_desired - 1) < normalized_fft_data(idx_desired) && normalized_fft_data(idx_desired + 1) < normalized_fft_data(idx_desired)
                    % fprintf('Peak found around the driving frequency. Storing data\n');
                    amplitude_vector = [amplitude_vector, amplitude]; % Pulls amplitude from fft calculation
                    initial_position_vector = [initial_position_vector, probe_data(1)];
                    phase_vector = [phase_vector, angle(normalized_fft_data(idx_desired))];
                    valid_probe_numbers = [valid_probe_numbers, nn];
                else
                    % fprintf('*** Alert: No peak found around the driving frequency. ***\n');
                end
            else
                % fprintf('*** Alert: No peak found around the driving frequency. ***\n');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Semi log plot (because exponential)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform linear fit
coefficients = polyfit(initial_position_vector, log(abs(amplitude_vector)), 1);

% Extract slope and intercept
slope = coefficients(1);
intercept = coefficients(2);

% Create a linear fit line
fit_line = exp(intercept) * exp(initial_position_vector.*slope);

% Convert coefficients to string
equation_str = sprintf('y = %.4f * exp(%.4f)', exp(intercept), slope);

% Plot original data and linear fit
figure;
semilogy(initial_position_vector, abs(amplitude_vector), 'bo', 'DisplayName', 'Data');
hold on;
semilogy(initial_position_vector, fit_line, 'r-', 'DisplayName', 'Linear Fit');
xlabel('Distance');
ylabel('Particle Oscillation Amplitude');
% title('Linear Fit of Attenuation of Oscillation in Probes', 'FontSize', 16);
    % Set the title with variables
title(sprintf('f=%.2f, k_n=%.2f, gamma_n=%.2f, P=%.2f, alpha=%.2f', driving_frequency, kn, gamma_n, dimensionless_p, slope), 'FontSize', 12);
legend('show');
grid on; 

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Wavenumber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unwrapped_phase_vector = unwrap(phase_vector);

% Plot initial position vs. phase as dots
figure;
scatter(initial_position_vector, unwrapped_phase_vector, 'o');
grid on;
hold on;  % Keep the plot for adding the fitted line

% Fit a line to the data
p = polyfit(initial_position_vector, unwrapped_phase_vector, 1);
fitted_line = polyval(p, initial_position_vector);

% Plot the fitted line
plot(initial_position_vector, fitted_line, '-r');

% Store the slope of the line as wavenumber
wavenumber = p(1);
wavespeed = driving_frequency/wavenumber;

% Label the axes
xlabel('z(t=0)');
ylabel('\Delta\phi');

% Customizing y-axis to show multiples of pi
y_max = max(unwrapped_phase_vector);  % Get the maximum y value
y_min = min(unwrapped_phase_vector);  % Get the minimum y value
yticks = [ceil(y_min/pi)*pi:pi:floor(y_max/pi)*pi];  % Define y-ticks in steps of pi
yticklabels = arrayfun(@(x) sprintf('%.2f\\pi', x/pi), yticks, 'UniformOutput', false);  % Create custom y-tick labels
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);  % Apply custom ticks and labels

% Set the title with variables
title(sprintf('f=%.2f, k_n=%.2f, gamma_n=%.2f, P=%.2f, k=%.2f', driving_frequency, kn, gamma_n, dimensionless_p, wavenumber), 'FontSize', 12);

% Hold off to finish the plotting
hold off;
