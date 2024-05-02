% function [asp_rat_counts,asp_rat_bins,rot_ang_counts,rot_ang_bins]=Sediment_Acoustics_2D_fixed_springs_nonaffine(K, M, Bv, w_D, Nt, N, P, W, seed, nnmax, plotdebug)
%% Molecular Dynamics Simulator (Adapted from Mark D. Shattuck, CCNY)

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

K = 100; % Spring Constant
M = 1; % Mass
Bv = 1; % gamma_n or dampening constant
w_D = 6.28; % oscillation frequency
Nt =1000; % Time steps
N = 1000; % Number of Particles
P=0.1; % Pressure
W = 5; % Width of channel
seed = 2;
plotdebug = 1;
nnmax = N;

% close all

K_in = K; % Why is this done? There is no K_in, maybe for loading
PackingName = ['N' num2str(N) '_P' num2str(P) '_Width' num2str(W) '_Seed' num2str(seed)];

load(['Packings/' PackingName '.mat']);

K = K_in; % Reassign the initial spring constant back to K after potentially modified by loaded settings - still not sure why this is done?
% N = length(Dn);

x_all = zeros([N,Nt]); % Pre allocating outside of loops because we know what size it should be
y_all = zeros([N,Nt]); 

flag = true;

Lx0 = Lx; % Store the initial value of the horizontal length of the simulation area, useful for boundary conditions or scaling.
B=0;
A = P_target/100;
dt = pi*sqrt(M/K)*0.05;

% Initialize old acceleration and velocity arrays to zero for all particles, preparing for the first computation step.
ax_old = 0*x;
ay_old = 0*y;
vx = 0*x;
vy = 0*y;
Ek = zeros(1,Nt);
Ep = zeros(1,Nt);
g = 0;

% initial positions
x0 = x;
y0 = y;

%% Make neighbor lists with initial spring lengths
skin = 0;
Zn_list = [];
neighbor_list_all = [];
spring_list_all = [];

for nn = 1:N
    % For each particle 'nn', initialize temporary lists to collect its neighbors and the lengths of springs connecting them.
    neighbor_list_nn = [];
    spring_list_nn = [];

    for mm = [1:nn-1,nn+1:N] % Check against all other particles 'mm' except itself.

        % Calculate the y-distance considering periodic boundary conditions in the vertical direction (Ly is the system height).
        dy = y(mm)-y(nn);
        dy = dy - round(dy/Ly)*Ly;

        % Calculate effective diameter for interaction, adjusted by 'skin'.
        Dnm = (1+skin)*(Dn(nn) + Dn(mm))/2;

        if(abs(dy) <= Dnm) % Only consider particle pairs within the vertical interaction range.
                        
            % Calculate horizontal distance and the square of the direct distance.
            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2;

            if(dnm < Dnm^2) % If within the circular interaction range,
                
                % Store the neighbor and the spring length (sqrt of the direct distance).
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
P = 0; % Why change the pressure? To "relax" the packing before 
for nt = 1:Nt


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% First step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_all(:,nt) = x;
    y_all(:,nt) = y;

    x  =  x+vx*dt+ax_old.*dt.^2/2;
    y  =  y+vy*dt+ay_old.*dt.^2/2;

    x(left_wall_list) = x0(left_wall_list)+A*sin(w_D*dt*nt);
    y(left_wall_list) = y0(left_wall_list);
    x(right_wall_list) = x0(right_wall_list);
    y(right_wall_list) = y0(right_wall_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Interaction detector and Force Law %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fx = zeros(1,N);
    Fy = zeros(1,N);
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
            dx = x(mm)-x(nn);
            dnm = dx.^2+dy.^2;
            %                 if(dnm < Dnm^2)
            dnm = sqrt(dnm);

            F = -K*(Dnm/dnm-1);

            % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
            %                     m_red = M*M/(M+M);

            dvx = vx(nn)-vx(mm);
            dvy = vy(nn)-vy(mm);

            Fx(nn) = Fx(nn)+F.*dx-Bv*dvx;  % particle-particle Force Law
            %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
            Fy(nn) = Fy(nn)+F.*dy-Bv*dvy;
            %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;


            %                     v_dot_r=((vx(nn)-vx(mm))*dx + (vy(nn)-vy(mm))*dy);
            %                     Fdiss = Bv * v_dot_r;
            %
            %                     Fx(nn) = Fx(nn)+F.*dx-Fdiss.*dx/dnm;  % particle-particle Force Law
            % %                     Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
            %                     Fy(nn) = Fy(nn)+F.*dy-Fdiss.*dy/dnm;
            % %                     Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
            Zn(nn) = Zn(nn) + 1;
            %                     Zn(mm) = Zn(mm) + 1;
            Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
            %                 end
            %             end
        end
    end
    %
    Fx = Fx - B.*vx;
    Fy = Fy - B.*vy;

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Second step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    vx = vx+(ax_old+ax).*dt/2;
    vy = vy+(ay_old+ay).*dt/2;

    ax_old = ax;
    ay_old = ay;
end

%%%%% Post Processing %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Snapshot of xy displacement for particles at ilist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tvec = (1:Nt)*dt;
omega_D = w_D;
[~,isort] = sort(x0);
ilist = [1,50,105,248,500,1000];
figure(1044)
tiledlayout(2,3,'TileSpacing','tight')
iskip = 10;
list = [];
b_start = 0;

% What is intent of this loop? Just debugging for particles in ilist?
for nn = isort(ilist)
    % if(~left_wall_list(nn))
    x_temp = x_all(nn,:);
    y_temp = y_all(nn,:);
    % if length(unique(x_temp))>100

    % Find the first index i0 where x_temp exceeds the midpoint between x0(nn) and its maximum value.
    i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');

    % Calculate the amplitude as the distance from the initial position.
    amp = sqrt((x_temp-x0(nn)).^2+(y_temp-y0(nn)).^2);

    % Calculate the angle from the initial position.
    ang = atan2((y_temp-y0(nn)),(x_temp-x0(nn)));

    % Shorten the time vector from the found index i0 to the end.
    t_temp = tvec(i0:end);

    % Compute the median amplitude.
    medamp = median(amp);

    % Create a time array for plotting.
    t = 1:length(amp);
    % figure(1040), plot(t(amp>medamp),atan(tan(ang(amp>medamp))),'o')
    % figure(1041), plot(t,amp,'o')
    % figure(1043), plot(x_temp,y_temp,'.'), axis equal
    
    % Check if plotting is enabled for debugging.
    if plotdebug

        % Focus on figure 1044 and plot the next tile.
        figure(1044),
        nexttile

        % Normalize and plot the positions by some scale A with equal axes.
        plot((x_temp-x0(nn))/A,(y_temp-y0(nn))/A,'.'), axis equal
        % nn
        title(['$x(t=0)=' num2str(round(x0(nn))) '$'],'Interpreter','latex')
        axis([-2 2 -2 2])
        xlabel('$x/A$','Interpreter','latex')
        ylabel('$y/A$','Interpreter','latex')
    end
    % pause
    % fo = fitoptions('Method','NonlinearLeastSquares',...
    %     'Lower',[0,-inf,omega_D*0.999,x0(nn)*0.999],...
    %     'Upper',[2*A,inf,omega_D*1.001,x0(nn)*1.001],...
    %     'StartPoint',[(max(x_temp(end-100:end))-x0(nn)) b_start omega_D x0(nn)]);
    % ft = fittype('a*sin(b-c*x)+d','options',fo);
    %
    % if i0+4<length(x_temp)
    %     [curve2,gof2] = fit(t_temp(i0:end)',amp(i0:end)',ft);
    %     if gof2.rsquare > 0.8 && curve2.a > A/100
    %         b_start = curve2.b;
    %         list = [list;[x0(nn),curve2.a,curve2.b]];
    %     end
    % end
    % end
    % end

end
% pause

tvec = (1:Nt)*dt;
omega_D = w_D;
[~,isort] = sort(x0);
iskip = 1;
list = [];
b_start = 0;

ellipse_stats = zeros(length(isort(1:iskip:end)),3);
nn_list = isort(1:iskip:nnmax);

for nn_counter = 1:length(nn_list)
    waitbar(nn_counter/length(nn_list));
    nn = nn_list(nn_counter);
    if(~left_wall_list(nn))
        x_temp = x_all(nn,:);
        y_temp = y_all(nn,:);
        if length(unique(x_temp))>100

            % find index when particles rises to .5 its max amplitude
            i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');

            % Calculate the amplitude as the Euclidean distance from the initial position (x0(nn), y0(nn)).
            amp = sqrt((x_temp-x0(nn)).^2+(y_temp-y0(nn)).^2);

            % Calculate the angle from the initial position to each point in Cartesian coordinates.
            ang = atan2((y_temp-y0(nn)),(x_temp-x0(nn)));
            t_temp = tvec(i0:end);

            medamp = median(amp);
            t = 1:length(amp);
            % figure(1040), plot(t(amp>medamp),atan(tan(ang(amp>medamp))),'o')
            % figure(1041), plot(t,amp,'o')
            % figure(1043), plot(x_temp,y_temp,'.'), axis equal
            % title(['$x(t=0)=' num2str(round(x0(nn))) '$'],'Interpreter','latex')

            % Normalize the x and y coordinates by subtracting the initial position and scaling by A.
            X = (x_temp(100:end)-x0(nn))/A;
            Y = (y_temp(100:end)-y0(nn))/A;

            % Why are these set to their same values? 
            X = X; 
            Y = Y;

            % Compute the radius squared for each point relative to the origin.
            R = X.^2+Y.^2;

            % Find the maximum radius squared from the array R.
            R_max = max(R);

            % Find the index of the first occurrence where R equals R_max.
            max_ind = find(R == R_max,1,"first");

            % Extract the X and Y coordinates corresponding to the maximum radius.
            X_max = X(max_ind);
            Y_max = Y(max_ind);

            % a0 = [R_max R_max atan(Y_max/X_max)];
            % a0 = [1 0.1 0];
            % options = optimset('Display','iter');
            % c = [0 0];
            % f = @(a) ((((X-c(1))*cos(a(3))+(Y-c(2))*sin(a(3))).^2)/a(1).^2 + (((X-c(1))*sin(a(3))-(Y-c(2))*cos(a(3))).^2)/a(2).^2 -1);
            % %[af fval] = fminsearch(f,a0);
            % options.FunctionTolerance = 1e-8;
            % [af,resnorm,residual,exitflag,output] = lsqnonlin(f, a0, [], [], options);
            % af = fminsearch(f,a0)
            af = [];

            % *** NEEDS fit_ellipse FUNCTION *** Fit an ellipse to the latter half of the data points, returning a structure with ellipse parameters.
            ellipse_t = fit_ellipse(X(round(end/2):end)',Y(round(end/2):end)');
            % nn_counter

            % Check if the ellipse fitting structure 'ellipse_t' is not empty and contains  parameter 'a' that comes from the fit_ellipse function
            if ~isempty(ellipse_t) && ~isempty(ellipse_t.a)
                % if nn_counter==336
                
                % end

                % Assign the semi-major axis, semi-minor axis, and negative rotation angle from 'ellipse_t' to 'af'
                af(1) = ellipse_t.a; % sub axis (radius) of the X axis of the non-tilt ellipse
                af(2) = ellipse_t.b; % sub axis (radius) of the Y axis of the non-tilt ellipse
                af(3) = -ellipse_t.phi; % orientation in radians of the ellipse (tilt)
            else
                % Define initial parameter estimates for fitting based on maximum radius and the angle at maximum point
                a0 = [R_max R_max atan(Y_max/X_max)];
                c = [0 0]; % Center of transformation set to origin

                % Define the function 'f' to minimize an ellipse fit error function parameterized by 'a'
                f = @(a) ((((X-c(1))*cos(a(3))+(Y-c(2))*sin(a(3))).^2)/a(1).^2 + (((X-c(1))*sin(a(3))-(Y-c(2))*cos(a(3))).^2)/a(2).^2 -1);
                %[af fval] = fminsearch(f,a0);

                % Use nonlinear least squares to optimize parameters 'a' starting from 'a0'
                af = lsqnonlin(f, a0);%, [], [], options);
            end

            if ~isempty(af)
                % af(1:2) = abs(af(1:2));

                % If the semi-major axis (af(1)) is less than the semi-minor axis (af(2))
                if af(1)<af(2)

                    % Rotate the ellipse by -pi/2 to correct the orientation
                    af(3) = af(3)-pi/2;

                    % Swap the semi-major and semi-minor axes to maintain convention
                    temp = af(2);
                    af(2) = af(1);
                    af(1) = temp;
                end

                % Normalize the rotation angle af(3) within the range -pi/2 to pi/2
                af(3) = (mod(af(3)+pi/2,pi)-pi/2);

                % plot(X,Y,'*'), hold on
                if plotdebug
                    figure(1044), clf, plot((x_temp(100:20:end)-x0(nn))/A,(y_temp(100:20:end)-y0(nn))/A,'ro'), axis equal, hold on
                    % plot(X_max,Y_max,'ro','markerfacecolor','r','markersize',20)

                    % plot(c(1), c(2), 'r*')
                    axis equal
                    theta_rot = af(3);
                    v = linspace(0, 1);
                    rx = af(1);
                    ry = af(2);
                    x = rx*cos(2*pi*v);
                    y = ry*sin(2*pi*v);
                    xx = x*cos(theta_rot) - y*sin(theta_rot);
                    yy = x*sin(theta_rot) + y*cos(theta_rot);
                    % figure
                    plot(xx, yy,'k-','linewidth',2)
                    plot([0 af(1)*cos(af(3))],[0 af(1)*sin(af(3))],'k-')
                    plot([0 -af(2)*sin(af(3))],[0 af(2)*cos(af(3))],'k-')
                    plot([0 af(1)*cos(af(3))],[0 0],'k--')
                    text(0.4, 0.1, '$a$','Interpreter','latex','fontsize',14)
                    text(-0.1, 0.02, '$b$','Interpreter','latex','fontsize',14)
                    text(0.8,0, '$\theta$','Interpreter','latex','fontsize',14)
                    
                    xlim([1.05*min((x_temp(100:20:end)-x0(nn))/A), 1.05*max((x_temp(100:20:end)-x0(nn))/A)])
                    ylim([1.05*min((y_temp(100:20:end)-y0(nn))/A), 1.05*max((y_temp(100:20:end)-y0(nn))/A)])
                    xlabel('$\Delta x / A$','Interpreter','latex','fontsize',18)
                    ylabel('$\Delta y / A$','Interpreter','latex','fontsize',18)
                    % grid
                    % axis([-1  1    -1  1])
                    % axis('equal')
                    % t=0:pi/10:2*pi;
                    % plot(c(1) + af(1)*cos(t), c(2) + af(2)*sin(t), 'r')
                    % axis([-2 2 -2 2])
                    drawnow
                    pause
                end

                % Store the fitted ellipse parameters (semi-major axis, semi-minor axis, rotation angle) in the ellipse_stats matrix for the current particle.
                ellipse_stats(nn_counter,:) = af;
            end
            % pause
            % fo = fitoptions('Method','NonlinearLeastSquares',...
            %     'Lower',[0,-inf,omega_D*0.999,x0(nn)*0.999],...
            %     'Upper',[2*A,inf,omega_D*1.001,x0(nn)*1.001],...
            %     'StartPoint',[(max(x_temp(end-100:end))-x0(nn)) b_start omega_D x0(nn)]);
            % ft = fittype('a*sin(b-c*x)+d','options',fo);
            %
            % if i0+4<length(x_temp)
            %     [curve2,gof2] = fit(t_temp(i0:end)',amp(i0:end)',ft);
            %     if gof2.rsquare > 0.8 && curve2.a > A/100
            %         b_start = curve2.b;
            %         list = [list;[x0(nn),curve2.a,curve2.b]];
            %     end
            % end
        end
    end
end

% Copy the ellipse_stats matrix to a new variable for processing.
ellipse_stats_nonzero = ellipse_stats;

% Convert the rotation angles from radians to degrees for easier interpretation.
ellipse_stats_nonzero(:,3) = ellipse_stats_nonzero(:,3)*180/pi;

% Ensure the semi-major and semi-minor axes values are positive.
ellipse_stats_nonzero(:,1:2) = abs(ellipse_stats_nonzero(:,1:2));

% Remove any rows where the semi-major axis is zero.
ellipse_stats_nonzero = ellipse_stats_nonzero(ellipse_stats_nonzero(:,1)~=0,:);

% Calculate the histogram of aspect ratios (semi-minor axis divided by semi-major axis) using bins of 0.05.
[asp_rat_counts,asp_rat_bins] = histcounts((ellipse_stats_nonzero(:,2))./(ellipse_stats_nonzero(:,1)),0:.05:1);

% Calculate the histogram of absolute rotation angles using bins of 5 degrees up to 90 degrees.
[rot_ang_counts,rot_ang_bins] = histcounts(abs(ellipse_stats_nonzero(:,3)),0:5:90);

if plotdebug
    figure
    histogram((ellipse_stats_nonzero(:,2))./(ellipse_stats_nonzero(:,1)),0:0.05:1);
    set(gca,'YScale','log')
    figure
    % h2 = histogram(mod(ellipse_stats_nonzero(:,3)/pi,floor(ellipse_stats_nonzero(:,3)/pi)),-1/2:0.05:1/2);
    histogram(abs(ellipse_stats_nonzero(:,3)),0:5:90);
    set(gca,'YScale','log')
end

% set(gca,'yscale','log')

% %%%% Output Data %%%%%
% uscore = (PackingName == '_');
% Filename = strcat('outputs2D_normdash/', PackingName, '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
% % Filename = strcat('outputs2D', PackingName(find(uscore,1,'first')+9:end-4), '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
% %Filename = strcat(PackingName(1:end-4), '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
% dlmwrite(Filename, list, 'precision', '%.8f');
%  %dlmwrite(strcat('Data2D/', Filename), x_all, 'precision', '%.8f');