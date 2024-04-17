function Sediment_Acoustics_3D_Packing_cell_list_newP(N, K, D, G, M, P_old, P_target_new, W_factor, seed, plotit)
%Function to create 2D packing with following input parameters:
% N, Number of Particles
% K, spring constant
% D, Average Diameter
% G, Ratio of large to small particles (typically 1.4)
% M, mass of particles
% P_thres, targeted threshold pressure
% W_factor, Factor of the width vs number of particles

%% Set up section
K_in = K;
PackingName = ['N' num2str(N) '_P' num2str(P_old) '_Width' num2str(W_factor) '_Seed' num2str(seed)];

load(['Packings3D/' PackingName '.mat']);

PackingNameNew = ['N' num2str(N) '_P' num2str(P_target_new) '_Width' num2str(W_factor) '_Seed' num2str(seed)];
filename = ['Packings3D_newP/' PackingNameNew '.mat'];

K = K_in;
% N = length(Dn);
flag = true;
Lx0 = Lx;
B=0.1;
Bv = 0.1; % dissipation factor

%% initial decompress to zero pressure
r_init = P;
Ly = Ly * (1+r_init);
y = y.*(1+r_init);

%% Physical parameters
g = 0;
P_target = P_target_new;
r = P_target;
flag = true;
flag2 = true;
fast_compress_flag = false;

%% Display Parameters
% plotit = 1;  % plot ?
plot_KE = 0;
Nplotskip = 200;  % number of timesteps to skip before plotting

%% Simulation Parmeters
dt = pi*sqrt(M/K)*0.05;
Nt = 1e8; % Number of steps

%% Initial Conditions
T = 0.001;
vx = sqrt(T)*randn(1,N);
vx = vx - mean(vx);
vy = sqrt(T)*randn(1,N);
vy = vy - mean(vy);
vz = sqrt(T)*randn(1,N);
vz = vz - mean(vz);

ax_old = 0*x;
ay_old = 0*y;
az_old = 0*z;
Ek = zeros(1,Nt);
Ep = zeros(1,Nt);

%% Verlet cell parameters
cell_width = 2*G*D;
jj_max = ceil(Lx/cell_width);

cell_num = ceil(x/cell_width);
for jj = 1:jj_max
    cell_list{jj} = find(cell_num == jj);
end

%% Setup Plotting
if plotit && 0
    figure(1), clf;
    h=zeros(1,2*N);
    for np = 1:N
        h(np) = rectangle('Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
        h(np+N)=rectangle('Position',[Lx Ly Dn(np) Dn(np)],'Curvature',[1 1],'edgecolor','b');
    end
    axis('equal');
    axis([0 Lx 0 Ly]);
    pause

    figure(2), clf;
%     figure(3), clf;
end
%% Main Loop
for nt = 1:Nt

    % plot particles
    if(plotit && mod(nt,Nplotskip) == 0)
        % if flag
        %     figure(1);
        %     for np = 1:N
        %         set(h(np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np) Dn(np) Dn(np)]);
        %     end
        %     Np=N;
        %     ii=find(y<Dn/2);
        %     for nn=1:length(ii)
        %         np=ii(nn);
        %         set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)+Ly Dn(np) Dn(np)]);
        %     end
        %     Np=Np+length(ii);
        %     %Top wall
        %     ii=find(y>Ly-Dn/2);
        %     for nn=1:length(ii)
        %         np=ii(nn);
        %         set(h(nn+Np),'Position',[x(np)-.5*Dn(np) y(np)-.5*Dn(np)-Ly Dn(np) Dn(np)]);
        %     end
        %     Np=Np+length(ii);
        %     figure(1), ylim([0, Ly]);
        %     title(num2str(Ly));
        % 
        % end
        figure(2), semilogy(nt, Ek(nt-1), 'ro');
        hold on, semilogy(nt, Ep(nt-1), 'bs');
        hold on, plot(nt,P,'kx')

%         figure(3), plot(nt,Ly,'ro',nt,Ly_min,'bx',nt,Ly_max,'kp'), hold on
        drawnow;

    elseif (plot_KE && mod(nt,Nplotskip) == 0)
        figure(3), plot3(x,y,z,'k.')
        axis equal
        drawnow
        pause
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% First step in Verlet integration %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x  =  x+vx*dt+ax_old.*dt.^2/2;
    y  =  y+vy*dt+ay_old.*dt.^2/2;
    z  =  z+vz*dt+az_old.*dt.^2/2;

    y = mod(y,Ly);
    z = mod(z,Lz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Re-assign particles to cells %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if flag2 == true || mod(nt,10000) == 0
        flag2 = false;
        cell_num = ceil(x/cell_width);
        for jj = 1:jj_max
            cell_list{jj} = find(cell_num == jj);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Interaction detector and Force Law %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Fx = zeros(1,N);
    Fy = zeros(1,N);
    Fz = zeros(1,N);
    Zn = zeros(1,N);

    for jj = 1:jj_max

        if jj == 1
            mm_list = [cell_list{1} cell_list{2}];
        elseif jj == jj_max
            mm_list = [cell_list{jj_max-1} cell_list{jj_max}];
        else
            mm_list = [cell_list{jj-1} cell_list{jj} cell_list{jj+1}];
        end

        for nn = cell_list{jj}

            mm_list_nn = mm_list(mm_list~=nn);

            for mm = mm_list_nn
                dy = y(mm)-y(nn);
                dy = dy - round(dy/Ly)*Ly;
                Dnm = (Dn(nn) + Dn(mm))/2;
                if(abs(dy) <= Dnm)
                    dz = z(mm)-z(nn);
                    dz = dz - round(dz/Lz)*Lz;
                
                    dx = x(mm)-x(nn);
                    dnm = dx.^2+dy.^2+dz.^2;
                    if(dnm < Dnm^2)
                        dnm = sqrt(dnm);

                        F = -K*(Dnm/dnm-1);

                        % dissipation force = B * m_red * v_N, v_N is normal component of velocity diff
                        m_red = M*M/(M+M);
                        v_dot_r=(vx(nn)-vx(mm))*dx + (vy(nn)-vy(mm))*dy + (vz(nn)-vz(mm))*dz;
                        Fdiss = Bv * m_red * v_dot_r;

                        Fx(nn) = Fx(nn)+F.*dx-Fdiss.*dx/dnm;  % particle-particle Force Law
                        %                         Fx(mm) = Fx(mm)-F.*dx+Fdiss.*dx/dnm;
                        Fy(nn) = Fy(nn)+F.*dy-Fdiss.*dy/dnm;
                        %                         Fy(mm) = Fy(mm)-F.*dy+Fdiss.*dy/dnm;
                        Fz(nn) = Fz(nn)+F.*dz-Fdiss.*dz/dnm;
                        Zn(nn) = Zn(nn) + 1;
                        %                         Zn(mm) = Zn(mm) + 1;
                        Ep(nt) = Ep(nt) + 0.5*K*(Dnm-dnm)^2;
                    end
                end
            end
        end

    end

    Fx = Fx - B.*vx;
    Fy = Fy - B.*vy;
    Fz = Fz - B.*vz;

    LW_contacts = x<Dn/2;
    Fx = Fx-K*(x-Dn/2).*(LW_contacts);  % Left wall
    % Fy = Fy-K*(y-D/2).*(y<D/2);  % Bottom wall

    RW_contacts = x>Lx-Dn/2;
    Fx = Fx-K*(x-(Lx-Dn/2)).*(RW_contacts);  % Right wall
    % Fy = Fy-K*(y-(Ly-D/2)).*(y>Ly-D/2);  % Top wall
    
    y=mod(y,Ly); %periodic boundaries for top and bottom

    Ek(nt) = 1/2*M*sum((vx).^2+(vy).^2+(vz).^2);
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
    
%     no_cont_list = (Zn == 0 & ~LW_contacts & ~RW_contacts);
%     vx(no_cont_list) = 0;
%     vy(no_cont_list) = 0;
%     ax(no_cont_list) = 0;
%     ay(no_cont_list) = 0;
% 

    ax_old = ax;
    ay_old = ay;
    az_old = az;

    tot_contacts = sum(Zn)/2;
    wall_contacts = sum(LW_contacts) + sum(RW_contacts);
    num_rattlers = sum(Zn==0);
    excess_contacts = tot_contacts + wall_contacts - 2*(N-num_rattlers)+1-1;

    P = sqrt(2*Ep(nt)/(K));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% COMPRESSION DECISIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if (P<P_target && Ek(nt)<1e-8)
        Ly = Ly * (1-r);
        y = y.*(1-r);
        flag = true;
        flag2 = true;
        nt_compress = nt;
    elseif(P>P_target && Ek(nt)<1e-12)%sum(Cn)/2>(3*sum(Cn>0)-2))
        break;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

% if nt == Nt
%     'hello, world!!!'
% end
% 
disp(['number of excess contacts = ' num2str(sum(Zn)/2 + sum(LW_contacts) + sum(RW_contacts) - 2*N)])

save(filename, 'x', 'y', 'z', 'Dn', 'Lx', 'Ly', 'Lz', 'K', 'P_target', 'P');