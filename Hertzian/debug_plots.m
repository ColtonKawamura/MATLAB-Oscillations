figure; % Create a new figure window
stem(freq_vector, abs(normalized_fft_data(index_vector)) * 2); % Plot the frequency spectrum
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
