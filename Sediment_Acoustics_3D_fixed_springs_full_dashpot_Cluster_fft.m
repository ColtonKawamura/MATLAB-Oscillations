function Sediment_Acoustics_3D_fixed_springs_full_dashpot_Cluster_fft(K, M, Bv, w_D, Nt, N, P, W, seed, tolerance)
%% Sediment_Acoustics_3D_fixed_springs_full_dashpot_Cluster_fft(100, 1, 7.59, 6.28, 300, 5000, 0.05, 5, 5, 0.05)

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
w_D = 1.28;
Nt =3000;
N = 5000;
P=0.05;
W = 5;
seed = 5;
tolerance = 5;


% close all
K_in = K;
PackingName = ['N' num2str(N) '_P' num2str(P) '_Width' num2str(W) '_Seed' num2str(seed)];

Filename = strcat('outputs3D_fulldash_v2/', PackingName, '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');

if exist(Filename)
    fprintf('*** Alert: Output for this already exists. ***\n')
    return
end

load(['Packings3D/' PackingName '.mat']);

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

skin = 0;
Zn_list = [];
neighbor_list_all = [];
spring_list_all = [];
for nn = 1:N
    neighbor_list_nn = [];
    spring_list_nn = [];
    for mm = [1:nn-1,nn+1:N]
        dy = y(mm)-y(nn);
        dy = dy - round(dy/Ly)*Ly;
        Dnm = (1+skin)*(Dn(nn) + Dn(mm))/2;
        if(abs(dy) <= Dnm)
            dz = z(mm)-z(nn);
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2;
            if(dnm < Dnm^2)

                neighbor_list_nn = [neighbor_list_nn, mm];
                spring_list_nn = [spring_list_nn, sqrt(dnm)];

            end
        end
    end
    neighbor_list_all{nn} = neighbor_list_nn;
    spring_list_all{nn} = spring_list_nn;
    Zn_list = [Zn_list;length(spring_list_nn)];
end


% identify wall particles
left_wall_list = (x<Dn/2);
right_wall_list = (x>Lx-Dn/2);
bulk_list = ~(left_wall_list | right_wall_list);

%% Main Loop
% P = 0;
for nt = 1:Nt


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
        spring_list = spring_list_all{nn};
        neighbor_list = neighbor_list_all{nn};
        for mm_counter = 1:length(neighbor_list)
            mm = neighbor_list(mm_counter);
            dy = y(mm)-y(nn);
            dy = dy - round(dy/Ly)*Ly;
            Dnm = spring_list(mm_counter);
            %             if(abs(dy) <= Dnm)
            dz = z(mm)-z(nn);
            dz = dz - round(dz/Lz)*Lz;

            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2+dz.^2;
            %                 if(dnm < Dnm^2)
            dnm = sqrt(dnm);

            F = -K*(Dnm/dnm-1);

            % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
            %                     m_red = M*M/(M+M);
            dvx = vx(nn)-vx(mm);
            dvy = vy(nn)-vy(mm);
            dvz = vz(nn)-vz(mm);

            Fx(nn) = Fx(nn)+F.*dx-Bv*dvx;  % particle-particle Force Law
            %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
            Fy(nn) = Fy(nn)+F.*dy-Bv*dvy;
            %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
            Fz(nn) = Fz(nn)+F.*dz-Bv*dvz;
            %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
            Zn(nn) = Zn(nn) + 1;
            %                     Zn(mm) = Zn(mm) + 1;
            Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
            %                 end
            %             end
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
    ay = Fy./M-g;
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

% Colton's FFT "engine" adaptation of existing variables
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

for nn = isort(1:iskip:end) # Sorts them by increments of iskip...for iskip>1, speeds things up
    if(~left_wall_list(nn)) # Executes the following block only if the particle indexed by nn is not on the left wall.
        x_temp = x_all(nn,:); # extracts and stores the time-series positions of the particle indexed by nn from the array x_all AKA looks at and picks out 
       
        if length(unique(x_temp))>10 # Checks if x_temp has more than 100 unique values AKA only process data that moves
            %  fprintf('*** Doing fft fit  ***\n');
            % FFT "engine"
            probe_data = x_temp;
            average_dt = mean(diff(time));
            sampling_freq = 1/average_dt;
            Fn = sampling_freq/2; % Nyquist frequency 
            number_elements_time = numel(time);
            centered_data = probe_data-mean(probe_data); %Center the data on zero for mean
            normalized_fft_data = fft(centered_data)/number_elements_time; 
            freq_vector = linspace(0, 1, fix(number_elements_time/2)+1)*Fn;
            index_vector = 1:numel(freq_vector);

            % Find the dominant frequency and its amplitude
            %   Need to double it because when signal is centered, power is
            %   distributed in both positive and negative. double the abs accounts for this
            [amplitude, idx_max] = max(abs(normalized_fft_data(index_vector)) * 2);
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
                    if sign_slope_before ~= sign_slope_after && normalized_fft_data(idx_desired - 1) < normalized_fft_data(idx_desired) && normalized_fft_data(idx_desired + 1) < normalized_fft_data(idx_desired) && abs(dominant_frequency - desired_frequency) < tolerance
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
                    %  fprintf('*** Alert: No peak found around the driving frequency. ***\n');
                end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Damping Warning for Data Anaylsis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the damping ratio
zeta = Bv / (2 * sqrt(M * K)); % Damping ratio formula

%% Check for overdamping
if zeta > 1
    warning('The system is overdamped! Damping ratio (zeta) is %f, which is greater than 1.', zeta);
elseif zeta == 1
    fprintf('The system is critically damped.\n');
else
    fprintf('The system is underdamped.\n');
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
ylabel('Probe Oscillation Amplitude');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change variable names for more verbosity 
attenuation = slope;
cleaned_particle_index = valid_probe_numbers; % Index of particles that made it through fft frequency “threshold” 

% Save the file
filename = sprintf('output_N%d_P%d_W%s_seed%d_K%d_Bv%d_wd%d_M%d.mat', N, P, W, seed, K, Bv, w_D, M);
save(filename, 'cleaned_particle_index', 'initial_position_vector', 'amplitude_vector', 'attenuation', 'unwrapped_phase_vector', 'wavenumber', 'wavespeed', 'N', 'P', 'W', 'seed', 'K', 'Bv', 'w_D', 'M');