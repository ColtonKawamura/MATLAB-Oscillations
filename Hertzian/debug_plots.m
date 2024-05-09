figure; % Create a new figure window
stem(freq_vector, abs(normalized_fft_data) * 2); % Plot the frequency spectrum
xlabel('Frequency (Hz)'); % Label the x-axis
ylabel('Amplitude'); % Label the y-axis
title('Frequency Spectrum'); % Title for the plot
grid on; % Turn on the grid for easier visualization

close all
figure;
plot3(x_all(nn,:), y_all(nn,:), z_all(nn,:), 'LineWidth', 2);
grid on;
title('3D Line Plot');
xlabel('X');
ylabel('Y');
zlabel('Z');

nn = 30
figure; % Create a new figure window
stem(frequency_spectra{nn}.Frequencies, frequency_spectra{nn}.Amplitudes * 2); % Plot the frequency spectrum
xlabel('Frequency (Hz)'); % Label the x-axis
ylabel('Amplitude'); % Label the y-axis
title(sprintf('Frequency Spectrum for Particle %d', frequency_spectra{nn}.ParticleIndex), 'FontSize', 12); % Title for the plot
grid on; % Turn on the grid for easier visualization

% Collect all frequencies into one array
all_frequencies = [];
for k = 1:numel(frequency_spectra)
    all_frequencies = [all_frequencies, frequency_spectra{k}.Frequencies];
end

% Plot histogram using hist
figure;
[n, x] = hist(all_frequencies, 40);  % 50 bins
bar(x, n);  % Plot as bar chart
title('Histogram of Frequencies','FontSize', 12);
xlabel('Frequency (Hz)');
ylabel('Count');
grid on;

% Find the index closest to x_0=target
target = 1.55; [~, nn]= min(abs(x0-target))
figure; plot(time_vector, y_all(nn,:))

% Calculate the number of elements and sampling frequency
number_elements_time = numel(time_vector);
average_dt = mean(diff(time_vector));
sampling_freq = 1 / average_dt;
nyquist_freq = sampling_freq / 2;
% FFT and Frequency Vector Setup
freq_vector = linspace(0, nyquist_freq, floor(number_elements_time/2)+1);
centered_data = position_nn - mean(position_nn);
normalized_fft_data = fft(centered_data) / number_elements_time;
magnitude_spectrum = abs(normalized_fft_data(1:floor(end/2)+1)) * 2;
% Plot the frequency spectrum
figure;
plot(freq_vector, magnitude_spectrum);
title(['Frequency Spectrum for Particle at Index ' num2str(nn)]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
grid on;