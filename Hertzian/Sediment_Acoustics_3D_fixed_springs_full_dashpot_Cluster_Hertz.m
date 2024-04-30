% function Sediment_Acoustics_3D_fixed_springs_full_dashpot_Cluster_Hertz(K, M, Bv, w_D, Nt, N, P, W, seed)
%% Molecular Dynamics Simulator (Adapted from Mark D. Shattuck, CCNY)
%%% For copy pasting: [~, index_of_closest]=min(abs(x0-53.53))

% Set up initial conditions and visualization
% Add random initial velocities
% Replace periodic boundaries with fixed walls
% Replace Euler with Velocity Verlet
% Add "Nplotskip" so every frame is not plotted
% Add gravity
%%%%% CHANGE FORCE-DETECTION TO SPEED UP %%%%%
% Add dissipation during collisions
% Add dt calculation based on sqrt(m/k)
% Add Ek(nt) storage inside loop

% % % Manual variables for troubleshooting
K = 100;
M = 1
Bv = 1;
w_D = 6.28;
Nt =500;
N = 10000;
P=0.05;
W = 5;
seed = 1;
tolerance = 0.5;



% close all
K_in = K;
PackingName = ['Hertz_' 'N' num2str(N) '_P' num2str(P) '_Width' num2str(W) '_Seed' num2str(seed)];

Filename = strcat('outputs3D_fulldash_Hertz/', PackingName, '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');

if exist(Filename)
    return
end

load(['Packings3D_newP_Hertz/' PackingName '.mat']);

K = K_in;
% N = length(Dn);
flag = true;
Lx0 = Lx;
D = min(Dn);
Lz = W*D;

B=0;



A = P_target/100;


dt = pi*sqrt(M/K)*0.05;
ax_old = 0*x;
ay_old = 0*y;
az_old = 0*z;
vx = 0*x;
vy = 0*y;
vz = 0*z;

Ek = zeros(1,Nt);
Ep = zeros(1,Nt);
g = 0;

x_all = zeros(length(x),Nt);
y_all = x_all;
z_all = x_all;

%% initial positions

x0 = x;
y0 = y;
z0 = z;

%% Make neighbor lists with initial spring lengths

skin = 0; % Add extra spac to the calculated sitance between particles, prob influencing collision and interaction calcs
Zn_list = []; %
neighbor_list_all = [];
Dnm_list_all = [];

for nn = 1:N
    percent_done = nn/N*100;
    fprintf('*** Currently %d percent done with building neighbor list. ***\n', percent_done);

    neighbor_list_nn = [];
    Dnm_list_nn = [];

    for mm = [1:nn-1,nn+1:N] % Runs over all other particles (mm) to check if they are neighbors of the current particle (nn)

        dy = y(mm)-y(nn); % Calcs the difference in the y-coords
        dy = dy - round(dy/Ly)*Ly; % also adjusts for periodic boundary conditions
        Dnm = (1+skin)*(Dn(nn) + Dn(mm))/2; % this is the interaction distance, can be adjust based on "skin"

        if(abs(dy) <= Dnm) % checks if the abs difference in y-coords is less than the interaction distance (Dnm)

            dz = z(mm)-z(nn); % if good, calculates distance in z and x
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2; % calcs euclidean distance dnm between particles

            if(dnm < Dnm^2) % checks if euclidean distance is within interaction range

                neighbor_list_nn = [neighbor_list_nn, mm]; % if good, adds to neighbor list for this particle
                Dnm_list_nn = [Dnm_list_nn, Dnm];

            end
        end
    end

    neighbor_list_all{nn} = neighbor_list_nn;
    Dnm_list_all{nn} = Dnm_list_nn;
    Zn_list = [Zn_list;length(Dnm_list_nn)];
end


% identify wall particles
left_wall_list = (x<Dn/2);
right_wall_list = (x>Lx-Dn/2);
bulk_list = ~(left_wall_list | right_wall_list);

%% Main Loop
% P = 0;
for nt = 1:Nt
    percent_done = nt/Nt*100;
    fprintf('*** Currently %d percent done with simulation. ***\n', percent_done);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% First step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_all(:,nt) = x;
    y_all(:,nt) = y;
    z_all(:,nt) = z;

    x  =  x+vx*dt+ax_old.*dt.^2/2;
    y  =  y+vy*dt+ay_old.*dt.^2/2;
    z  =  z+vz*dt+az_old.*dt.^2/2;

    x(left_wall_list) = x0(left_wall_list)+A*sin(w_D*dt*nt);
    y(left_wall_list) = y0(left_wall_list);
    z(left_wall_list) = z0(left_wall_list);
    x(right_wall_list) = x0(right_wall_list);
    y(right_wall_list) = y0(right_wall_list);
    z(right_wall_list) = z0(right_wall_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Interaction detector and Force Law %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fx = zeros(1,N);
    Fy = zeros(1,N);
    Fz = zeros(1,N);
    Zn = zeros(1,N);

    for nn = 1:N
        Dnm_list = Dnm_list_all{nn};
        neighbor_list = neighbor_list_all{nn};
        for mm_counter = 1:length(neighbor_list)
            mm = neighbor_list(mm_counter);
            dy = y(mm)-y(nn);
            dy = dy - round(dy/Ly)*Ly;
            Dnm = Dnm_list(mm_counter);
            %             if(abs(dy) <= Dnm)
            dz = z(mm)-z(nn);
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2;
            %                 if(dnm < Dnm^2)
            dnm = sqrt(dnm);
            if dnm< Dnm

                F = -K*(Dnm-dnm).^(1.5);

                % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
                %                     m_red = M*M/(M+M);
                dvx = vx(nn)-vx(mm);
                dvy = vy(nn)-vy(mm);
                dvz = vz(nn)-vz(mm);

                Fx(nn) = Fx(nn)+F.*dx/dnm-Bv*dvx;  % particle-particle Force Law
                %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
                Fy(nn) = Fy(nn)+F.*dy/dnm-Bv*dvy;
                %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
                Fz(nn) = Fz(nn)+F.*dz/dnm-Bv*dvz;
                %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
                Zn(nn) = Zn(nn) + 1;
                %                     Zn(mm) = Zn(mm) + 1;
                Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
                %                 end
                %             end
            end
        end
    end
    %
    % Fx = Fx - B.*vx;
    % Fy = Fy - B.*vy;
    % Fz = Fz - B.*vz;

    Fx(left_wall_list) = 0;
    Fx(right_wall_list) = 0;


    %     Fx = Fx-K*(x-Dn/2).*(x<Dn/2);  % Left wall
    % Fy = Fy-K*(y-D/2).*(y<D/2);  % Bottom wall

    %     Fx = Fx-K*(x-(Lx-Dn/2)).*(x>Lx-Dn/2);  % Right wall
    % Fy = Fy-K*(y-(Ly-D/2)).*(y>Ly-D/2);  % Top wall

    %     y=mod(y,Ly); %periodic boundaries for top and bottom

    Ek(nt) = 1/2*M*sum((vx).^2+(vy).^2);
    Ek(nt) = Ek(nt)/N;
    Ep(nt) = Ep(nt)/N;

    ax = Fx./M;
    ay = Fy./M;
    az = Fz./M;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Second step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vx = vx+(ax_old+ax).*dt/2;
    vy = vy+(ay_old+ay).*dt/2;
    vz = vz+(az_old+az).*dt/2;

    ax_old = ax;
    ay_old = ay;
    az_old = az;
end

%%%%% Post Processing %%%%%

tvec = (1:Nt)*dt;
omega_D = w_D;
[~,isort] = sort(x0);
iskip = 1;
list = [];
b_start = 0;
offset_guess = 0;

problem_children_index = [];
time_vector = tvec;
driving_frequency = omega_D/6.2832;
kn = K;
gamma_n = Bv;
dimensionless_p=P;
driving_amplitude=A;

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];
initial_phase_offset = 0;


for nn = isort(1:iskip:end)
    fprintf('*** Analyzing particle %d . ***\n', nn)
    if(~left_wall_list(nn))
        % fprintf('*** Particle %d not a wall particle. ***\n', nn);
        x_temp = x_all(nn,:);

        if length(unique(x_temp))>10
            % fprintf('*** Particle %d has more than 10 unique values. ***\n', nn);

            i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');

            if ~isempty(i0) % If there's not anything wrong with the inital oscilation index, move on

                initial_oscilation_time_vector = tvec(i0:end);
                index_initial_oscillation = i0;
                fit_x = initial_oscilation_time_vector;
                fit_y = x_temp(i0:end)-x0(nn);

                % Center the data on y = 0
                mean_y = mean(fit_y);
                centered_fit_y = fit_y - mean_y;

                % Detect zero crossings
                shifted_y = circshift(centered_fit_y, [0, 1]); % % Shift centered_fit_y circularly by one position to the right
                zero_crossing_indices = find(centered_fit_y .* shifted_y <= 0); % Find indices where the product of adjacent elements is less than or equal to zero

                if numel(zero_crossing_indices)>=3 % If particle nn has more that three zero-crossings with it's mean, move on
                    % Define function to fit
                    fit_function = @(b, X) b(1) .* sin(2 * pi * driving_frequency * X - b(2));

                    % Define least-squares cost function
                    cost_function = @(b) sum((fit_function(b, fit_x) - fit_y).^2);

                    % Initial guess
                    initial_amplitude =  driving_amplitude;
                    initial_guess = [initial_amplitude; initial_phase_offset];

                    % Initial guess            
                    initial_amplitude =  driving_amplitude;
                    initial_guess = [initial_amplitude; initial_phase_offset];

                    % Perform fitting
                    [s, ~, ~] = fminsearch(cost_function, initial_guess);

                    %Update the initial phase offset for the next iteration
                    initial_phase_offset = s(2);

                    % Calculate correlation coefficient
                    [R, ~] = corrcoef(fit_function(s, fit_x), fit_y);
                    R2 = R(1, 2)^2;

                    % Check if R2 is less than 0.5
                    if R2 < 0.7
                        % fprintf('*** Particle %d did not make the cut. ***\n', nn); 
                    else
                        % If good, store the vector and probe number
                        initial_position_vector = [initial_position_vector, x_temp(1)];
                        amplitude_vector = [amplitude_vector, s(1)];
                        phase_vector = [phase_vector, s(2)];
                        cleaned_particle_index = [cleaned_particle_index, nn];
                        % For Troubleshooting, plot the oscilations to make sure the simulation is working
                        plot(time_vector,x_all(nn,:)-mean(x_all(nn,:)))
                        legend([' nn = ' num2str(nn)]);
                        hold on
                    end
                else % This let's use know the particle index of the particles that have something wrong.
                    % fprintf('*** Empty Vector Found. Storing Index for Debugging. ***\n');
                    problem_children_index = [problem_children_index, nn];
                end
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
simulation_attenuation = coefficients(1);
intercept = coefficients(2);

% Create a linear fit line
fit_line = exp(intercept) * exp(initial_position_vector.*simulation_attenuation);

% Plot original data and linear fit
figure;
semilogy(initial_position_vector, abs(amplitude_vector), 'bo', 'DisplayName', 'Data');
hold on;
semilogy(initial_position_vector, fit_line, 'r-', 'DisplayName', 'Linear Fit');
xlabel('Distance');
ylabel('Particle Oscillation Amplitude');
title_str = sprintf('f=%d, k_n=%d, gamma_n=%d, P=%d, alpha=%d', driving_frequency, kn, gamma_n, P, simulation_attenuation);
title(title_str, 'FontSize', 12);
% title(sprintf('f=%d, k_n=%d, gamma_n=%d, P=%d, alpha=%d', driving_frequency, kn, gamma_n, P, simulation_attenuation), 'FontSize', 12);
legend('show');
grid on;
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fft fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tvec = (1:Nt)*dt;
omega_D = w_D;
[~,isort] = sort(x0);
iskip = 50;
list = [];
b_start = 0;
offset_guess = 0;

problem_children_index = [];
time_vector = tvec;
driving_frequency = omega_D/6.2832;
kn = K;
gamma_n = Bv;
dimensionless_p=P;
driving_amplitude=A;
freq_match_tolerance = 0.05;

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];
initial_phase_offset = 0;


for nn = isort(1:iskip:end) % Sorts them by increments of iskip...for iskip>1, speeds things up
    if(~left_wall_list(nn)) % Executes the following block only if the particle indexed by nn is not on the left wall.
        x_temp = x_all(nn,:); % extracts and stores the time_vector-series positions of the particle indexed by nn from the array x_all AKA looks at and picks out 
       
        if length(unique(x_temp))>10 % Checks if x_temp has more than 100 unique values AKA only process data that moves
            %  fprintf('*** Doing fft fit  ***\n');
            % FFT "engine"
            particle_position = x_temp;
            average_dt = mean(diff(time_vector));
            sampling_freq = 1/average_dt;
            Fn = sampling_freq/2; % Nyquist frequency 
            number_elements_time = numel(time_vector);
            centered_data = particle_position-mean(particle_position); %Center the data on zero for mean
            normalized_fft_data = fft(centered_data)/number_elements_time; 
            freq_vector = linspace(0, 1, fix(number_elements_time/2)+1)*Fn;
            index_vector = 1:numel(freq_vector);

            % Find the dominant frequency and its max amplitude
            %   Need to double it because when signal is centered, power is
            %   distributed in both positive and negative. double the abs accounts for this
            [max_particle_amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
            dominant_frequency = freq_vector(idx_max);

            % Find the index of the frequency closest to driving frequency
            desired_frequency = driving_frequency;

            % Find the index of the closest frequency to the desired frequency
            [~, idx_desired] = min(abs(freq_vector - desired_frequency));

            % Check if there is a peak around the desired frequency and amplitude is greater than 
                if idx_desired > 1 && idx_desired < numel(freq_vector) 
                    %  fprintf('*** Checking for Slope  ***\n');
                    % Calculate the sign of the slope before and after the desired frequency
                    sign_slope_before = sign(normalized_fft_data(idx_desired) - normalized_fft_data(idx_desired - 1));
                    sign_slope_after = sign(normalized_fft_data(idx_desired + 1) - normalized_fft_data(idx_desired));
                    
                    % Check if the signs of the slopes are different and if the values on both sides are greater than the value at the desired frequency
                    if sign_slope_before ~= sign_slope_after && abs(normalized_fft_data(idx_desired - 1)) < abs(normalized_fft_data(idx_desired)) && abs(normalized_fft_data(idx_desired + 1)) < abs(normalized_fft_data(idx_desired)) && abs(dominant_frequency - desired_frequency) < freq_match_tolerance
                        %  fprintf('Peak found around the driving frequency. Storing data\n');

                        amplitude_vector = [amplitude_vector, max_particle_amplitude]; % Pulls amplitude from fft calculation
                        initial_position_vector = [initial_position_vector, particle_position(1)];
                        phase_vector = [phase_vector, angle(normalized_fft_data(idx_desired))];
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
% Semi log plot (because exponential)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform linear fit
coefficients = polyfit(initial_position_vector, log(abs(amplitude_vector)), 1);

% Extract fitted_attenuation and intercept
fitted_attenuation = coefficients(1);
intercept = coefficients(2);

% Create a linear fit line
fit_line = exp(intercept) * exp(initial_position_vector.*fitted_attenuation);

% Convert coefficients to string
equation_str = sprintf('y = %.4f * exp(%.4f)', exp(intercept), fitted_attenuation);

% Plot original data and linear fit
figure;
semilogy(initial_position_vector, abs(amplitude_vector), 'bo', 'DisplayName', 'Data');
hold on;
semilogy(initial_position_vector, fit_line, 'r-', 'DisplayName', 'Linear Fit');
xlabel('Distance');
ylabel('Particle Oscillation Amplitude');
% title('Linear Fit of Attenuation of Oscillation in Probes', 'FontSize', 16);
    % Set the title with variables
title(sprintf('f=%.2f, k_n=%.2f, gamma_n=%.2f, P=%.2f, alpha=%.2f', driving_frequency, kn, gamma_n, dimensionless_p, fitted_attenuation), 'FontSize', 12);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harmonic analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tvec = (1:Nt)*dt;
omega_D = w_D;
[~,isort] = sort(x0);
iskip = 10;
list = [];
b_start = 0;
offset_guess = 0;

problem_children_index = [];
time_vector = tvec;
driving_frequency = omega_D/6.2832;
kn = K;
gamma_n = Bv;
dimensionless_p=P;
driving_amplitude=A;
freq_match_tolerance = 0.05;

% Initialize output vectors
initial_position_vector = [];
amplitude_vector = [];
phase_vector = [];
cleaned_particle_index = [];
initial_phase_offset = 0;
frequency_spectra = {};  % Cell array to store frequency spectra of particles with more than two peaks

for nn = isort(1:iskip:end) % Sorts them by increments of iskip...for iskip>1, speeds things up
    if(~left_wall_list(nn)) % Executes the following block only if the particle indexed by nn is not on the left wall.
        x_temp = x_all(nn,:); % extracts and stores the time_vector-series positions of the particle indexed by nn from the array x_all AKA looks at and picks out 
       
        if length(unique(x_temp))>10 % Checks if x_temp has more than 100 unique values AKA only process data that moves
            %  fprintf('*** Doing fft fit  ***\n');
            % FFT "engine"
            particle_position = x_temp;
            average_dt = mean(diff(time_vector));
            sampling_freq = 1/average_dt;
            Fn = sampling_freq/2; % Nyquist frequency 
            number_elements_time = numel(time_vector);
            centered_data = particle_position-mean(particle_position); %Center the data on zero for mean
            normalized_fft_data = fft(centered_data)/number_elements_time; 
            freq_vector = linspace(0, 1, fix(number_elements_time/2)+1)*Fn;
            index_vector = 1:numel(freq_vector);

            % Find the dominant frequency and its max amplitude
            %   Need to double it because when signal is centered, power is
            %   distributed in both positive and negative. double the abs accounts for this
            [max_particle_amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
            dominant_frequency = freq_vector(idx_max);

            % Find the index of the frequency closest to driving frequency
            desired_frequency = driving_frequency;

            % Find the index of the closest frequency to the desired frequency
            [~, idx_desired] = min(abs(freq_vector - desired_frequency));

            % Check if there is a peak around the desired frequency 
            if idx_desired > 1 && idx_desired < numel(freq_vector) 
                %  fprintf('*** Checking for Slope  ***\n');
                % Calculate the sign of the slope before and after the desired frequency
                sign_slope_before = sign(normalized_fft_data(idx_desired) - normalized_fft_data(idx_desired - 1));
                sign_slope_after = sign(normalized_fft_data(idx_desired + 1) - normalized_fft_data(idx_desired));
                
                % Check if the signs of the slopes are different and if the values on both sides are greater than the value at the desired frequency
                if sign_slope_before ~= sign_slope_after && abs(normalized_fft_data(idx_desired - 1)) < abs(normalized_fft_data(idx_desired)) && abs(normalized_fft_data(idx_desired + 1)) < abs(normalized_fft_data(idx_desired)) && abs(dominant_frequency - desired_frequency) < freq_match_tolerance
                    %  fprintf('Peak found around the driving frequency. Storing data\n');

                    amplitude_vector = [amplitude_vector, max_particle_amplitude]; % Pulls amplitude from fft calculation
                    initial_position_vector = [initial_position_vector, particle_position(1)];
                    phase_vector = [phase_vector, angle(normalized_fft_data(idx_desired))];
                    cleaned_particle_index = [cleaned_particle_index, nn];
                else
                    %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');
                end
            else
            end
            number_bins = length(normalized_fft_data);
            found_peaks = [];
            for i = 2:number_bins-1 % starts at 2 and ends at end-1 because can't check slopes around the last and first points
                if abs(normalized_fft_data(i)) > abs(normalized_fft_data(i-1)) && abs(normalized_fft_data(i)) > abs(normalized_fft_data(i+1))
                    found_peaks = [found_peaks, i];
                end
            end
            if length(found_peaks) > 2
                frequency_spectra{end+1} = struct('ParticleIndex', nn, 'Amplitudes', abs(normalized_fft_data(index_vector)), 'Frequencies', freq_vector);
                fprintf('*** Alert: More Frequency Peaks Found. ***\n');
            end
        end
    end
end