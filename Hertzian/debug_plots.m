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