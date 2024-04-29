function Sediment_Acoustics_3D_fixed_springs_full_dashpot_Cluster(K, M, Bv, w_D, Nt, N, P, W, seed)
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


% close all
K_in = K;
PackingName = ['N' num2str(N) '_P' num2str(P) '_Width' num2str(W) '_Seed' num2str(seed)];

Filename = strcat('outputs3D_fulldash/', PackingName, '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');

if exist(Filename)
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
for nn = isort(1:iskip:end)
    if(~left_wall_list(nn))
        x_temp = x_all(nn,:);

        if length(unique(x_temp))>100

            i0 = find(x_temp>x0(nn)+0.5*(max(x_temp)-x0(nn)),1,'first');

            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,-inf,omega_D*0.999,x0(nn)*0.999],...
                'Upper',[2*A,inf,omega_D*1.001,x0(nn)*1.001],...
                'StartPoint',[(max(x_temp(end-100:end))-x0(nn)) b_start omega_D x0(nn)]);
            ft = fittype('a*sin(b-c*x)+d','options',fo);

            if i0+4<length(x_temp)
                [curve2,gof2] = fit(tvec(i0:end)',x_temp(i0:end)',ft);
                if gof2.rsquare > 0.8 && curve2.a > A/100
                    b_start = curve2.b;
                    list = [list;[x0(nn),curve2.a,curve2.b]];
                end
            end
        end
    end

end

%%%% Output Data %%%%%
uscore = (PackingName == '_');
Filename = strcat('outputs3D_fulldash/', PackingName, '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
% Filename = strcat('outputs2D', PackingName(find(uscore,1,'first')+9:end-4), '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
%Filename = strcat(PackingName(1:end-4), '_K', num2str(K), '_Bv', num2str(Bv), '_wD', num2str(w_D), '_M', num2str(M), '.dat');
dlmwrite(Filename, list, 'precision', '%.8f');
%dlmwrite(strcat('Data2D/', Filename), x_all, 'precision', '%.8f');