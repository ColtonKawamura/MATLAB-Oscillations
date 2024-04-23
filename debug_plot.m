time = tvec;
driving_frequency = omega_D/6.2832;
kn = K;
gamma_n = Bv;
dimensionless_p=P;
driving_amplitude=A;

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
valid_probe_numbers = [];
initial_phase_offset = 0;

    if(~left_wall_list(nn)) # Executes the following block only if the particle indexed by nn is not on the left wall.
        x_temp = x_all(nn,:); # extracts and stores the time-series positions of the particle indexed by nn from the array x_all AKA looks at and picks out 
       
        if length(unique(x_temp))>100 # Checks if x_temp has more than 100 unique values AKA only process data that moves
            %  fprintf('*** Doing fft fit  ***\n');
            probe_data = x_temp;
            average_dt = mean(diff(time));
            sampling_freq = 1/average_dt;
            Fn = sampling_freq / 2; % Nyquist frequency

            % Determine the index to start from the second half of the data / Round down to lowest number of half elements in probe_data ,then +1 because MATLAB is one-based indexed, so shifts start_index to start of 2nd half
            start_index = floor(numel(probe_data) / 2) + 1;

            % Select only the last half of the probe_data
            last_half_data = probe_data(start_index:end);

            % Adjust the number of elements to reflect the 2nd half of data
            half_number_elements_time = numel(last_half_data);

            % Center the data on zero for mean
            centered_data = last_half_data - mean(last_half_data);

            % Perform FFT normalized by the number of elements
            normalized_fft_data = fft(centered_data) / half_number_elements_time;

            % Frequency vector calculation adjusted for the new data length
            freq_vector = linspace(0, 1, fix(half_number_elements_time / 2) + 1) * Fn;

            % Index vector for the frequencies
            index_vector = 1:numel(freq_vector);

            % Find the dominant frequency and its amplitude
            %   Need to double it because when signal is centered and normalized, original energy at each f represented equaly of both sides of 0.  Looking at positive values of f, so need to double for correct amplitude.
            %   distributed in both positive and negative. double the abs accounts for this
            [amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
            dominant_frequency = freq_vector(idx_max);

            % Find the index of the frequency closest to driving frequency
            desired_frequency = driving_frequency;

            % Find the index of the closest frequency to the desired frequency
            [~, idx_desired] = min(abs(freq_vector - desired_frequency));

            % Check if there is a peak around the desired frequency and amplitude is greater than 
                if idx_desired > 1 && idx_desired < numel(freq_vector) && amplitude > cutoff_amplitude
                    %  fprintf('*** Checking for Slope  ***\n');
                    % Calculate the sign of the slope before and after the desired frequency
                    sign_slope_before = sign(normalized_fft_data(idx_desired) - normalized_fft_data(idx_desired - 1));
                    sign_slope_after = sign(normalized_fft_data(idx_desired + 1) - normalized_fft_data(idx_desired));
                    
                    % Check if the signs of the slopes are different and if the values on both sides are greater than the value at the desired frequency
                    if sign_slope_before ~= sign_slope_after && normalized_fft_data(idx_desired - 1) < normalized_fft_data(idx_desired) && normalized_fft_data(idx_desired + 1) < normalized_fft_data(idx_desired)
                        %  fprintf('Peak found around the driving frequency. Storing data\n');

                        amplitude_vector = [amplitude_vector, amplitude]; % Pulls amplitude from fft calculation
                        initial_position_vector = [initial_position_vector, probe_data(1)];
                        phase_vector = [phase_vector, angle(normalized_fft_data(idx_desired))];
                        % phase_vector = [phase_vector, s(2)];
                        valid_probe_numbers = [valid_probe_numbers, nn];
                    else
                        %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');
                    end
                else
                    %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');mail
                end
        end
    end

% Plot Frequency Spectrum of Particle nn
figure; % Create a new figure window
stem(freq_vector, abs(normalized_fft_data(index_vector)) * 2); % Plot the frequency spectrum
xlabel('Frequency (Hz)'); % Label the x-axis
ylabel('Amplitude'); % Label the y-axis
title(sprintf('Frequency Spectrum of Particle %f', nn), 'FontSize', 12); % Title for the plot with font size
grid on; % Turn on the grid for easier visualization

% Plot Oscillation of Particle nn
figure;
plot(time, x_all(nn,:) - x0(nn));
xlabel('Time'); % Label the x-axis
ylabel('Position'); % Label the y-axis
title(sprintf('Full Time Oscillation of Particle %f', nn), 'FontSize', 12); % Set title with font size
grid on; % Turn on the grid for easier visualization

% Plot 2nd Half Oscillation of Particle nn
figure;
plot(time(start_index:end), x_all(nn, start_index:end) - x0(nn));
xlabel('Time'); % Label the x-axis
ylabel('Position'); % Label the y-axis
title(sprintf('Second Half of Oscillation of Particle %f', nn), 'FontSize', 12); % Set title with font size
grid on; % Turn on the grid for easier visualization

