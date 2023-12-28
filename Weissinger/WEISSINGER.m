
close all
clear
clc


flpath = pwd;
addpath(append(flpath,'/main/'));

tic

% AERODYNAMIC properties
alpha     = 0;
beta      = 0;

% GEOMETRY
delta     = 0;
lambda    = 0;
root      = 8;
L         = 15;
taper     = 1;
AOA       = 0;

% DISCRETIZATION 
M         = 10;
N         = 5;
alpha_vec = linspace(-4,8,30);

flag = "plot";

% Panel creation function
[PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

% System matrix generation

toll        = 1e-4;
[MATRIX]    = BIOTSAVART(PANELwing,AOA,M,N,L,toll);

% Known term vector
[b]         = VECTORcomp(PANELwing,alpha,beta,M,N);

% System solving
GAMMA       = MATRIX\b;

% Plotting spanwise GAMMA distribution
GAMMAplot(GAMMA,M,N);

% LIFT computation
U                = 1;
rho              = 1;
S                = (root + root/taper) * L*cos(lambda/180*pi);
[~,L_vec,Cl,~,~] = LIFTcomp(GAMMA,PANELwing,lambda,M,N,rho,U,S,"yes");

% Induced velocity computation
[v_ind,alpha_ind] = INDUCEDVELcomp(GAMMA,PANELwing,alpha,M,N,U,"yes");

% DRAG computation
[D,D_vec,Cd]      = DRAGcomp(L_vec,-alpha_ind,alpha,rho,U,S,M);

% Cl and Cd - alpha computation
[Cl_vec,Cd_vec]   = CLAplot(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,"yes");

toc

% Removing path
flpath = pwd;
rmpath(append(flpath,'/main/'));

%% Coefficients computation varying wing geometry 
close all
clear 
clc

% Initializing path
flpath = pwd;
addpath(append(flpath,'/main/'));

tic

% TAPER RATIO study
% GEOMETRY
    delta     = 0;
    lambda    = 0;
    root      = 8;
    L         = 30;
    AOA       = 0;

% AERODYNAMIC properties
    beta      = 0;

% DISCRETIZATION
    M         = 10;
    N         = 5;
    alpha_vec = linspace(0,5,30);

% TAPER RATIO discretization properties
    TAPERvec  = 1:1:5;

% Plotting options
flag = "noplot";

figure 

subplot(3,1,1)
hold on

for taper = TAPERvec
    % Panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % System matrix generation
    % Setting tolerance to avoid singular MATRIX 
    toll        = 1e-4;
    [MATRIX]    = BIOTSAVART(PANELwing,AOA,M,N,L,toll);
    
    % Surfaces computation
    S               = (root + root/taper) * L*cos(lambda/180*pi);
    
    % Cl and Cd - alpha computation
    [Cl_vec,Cd_vec] = CLAplot(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
    % Plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
end

grid on
grid minor

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\lambda \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\lambda = " + string(TAPERvec) + "$";
legend(TEXT,'Interpreter','latex','location','best');

% DIHEDRAL study
% GEOMETRY 
    lambda    = 0;
    root      = 8;
    L         = 30;
    taper     = 1;
    AOA       = 0;

% DIHEDRAL discretization 
    DELTAvec  = 0:1:5;

% Plotting options
flag = "noplot";

subplot(3,1,2)
hold on

for delta = DELTAvec
    % Panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % System matrix generation
    
    toll        = 1e-4;
    [MATRIX]    = BIOTSAVART(PANELwing,AOA,M,N,L,toll);
    
    % Surfaces computation
    S = (root + root/taper) * L*cos(lambda/180*pi);
    
    % Cl-alpha and Cd-alpha computation
    [Cl_vec,Cd_vec] = CLAplot(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
    % Plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
    drawnow
    
end

grid on
grid minor

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\Delta \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\Delta = " + string(DELTAvec) + "$";
legend(TEXT,'Interpreter','latex','location','best');

% SWEEP ANGLE study
    delta     = 0;
    root      = 8;
    L         = 30;
    taper     = 1;
    AOA       = 0;

% DISCRETIZATION 
    LAMBDAvec = 0:10:30;

% Plotting options
flag = "noplot";

subplot(3,1,3)
hold on

for lambda = LAMBDAvec
    
    % Panel creation function 
    [PANELwing] = PANELING(delta,lambda,AOA,root,taper,L,M,N,flag,[0,0,0]);

    % System matrix generation
    
    toll        = 1e-4;
    [MATRIX]    = BIOTSAVART(PANELwing,AOA,M,N,L,toll);
    
    % Surfaces computation
    S = (root + root/taper) * L*cos(lambda/180*pi);
    
    % Cl-alpha and Cd-alpha computation
    [Cl_vec,Cd_vec] = CLAplot(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag);
    
    % Plotting results
    plot(Cd_vec,Cl_vec,'LineWidth',3);
    drawnow

end

grid on

xlabel('$C_{D}$','Interpreter','latex');
ylabel('$C_{L}$','Interpreter','latex');
TEXT = "$\Lambda \ , \ \beta = " + string(beta) + "$"; 
title(TEXT,'Interpreter','latex');
TEXT = "$\Lambda = " + string(LAMBDAvec) + "$";
legend(TEXT,'Interpreter','latex','location','best');

toc

% Removing path
flpath = pwd;
rmpath(append(flpath,'/main/'));

%% Wing interactions
close all
clear
clc


flpath = pwd;
addpath(append(flpath,'/main/'));

tic

% AERODYNAMIC ANGLES
    alpha     = 0;
    beta      = 0;
    alpha_vec = linspace(-10,10,30);
    
% 1ST WING GEOMETRY     - A320
    delta1    = 5;
    lambda1   = 25;
    root1     = 6;
    L1        = 34;
    taper1    = 2;
    AOA1      = 5;
    transl1   = [0,0,0]; 
    M1        = 7;
    N1        = 3;
    
% 2ND WING GEOMETRY (tail)
    delta2    = 15;
    lambda2   = 5;
    root2     = 4;
    L2        = 10;
    taper2    = 1;
    AOA2      = 5;
    transl2   = [35,0,0]; 
    M2        = 7;
    N2        = 3;

% Plotting options
flag = "plot";

% Panel creation function 
[PANELwing1] = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,transl1);
[PANELwing2] = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,transl2);

% Assemblying MATRIX
tol          = 1e-4;
PANELwing    = [PANELwing1,PANELwing2];
[MATRIX]     = BIOTSAVART_multi(PANELwing,[AOA1,AOA2],[M1,M2],[N1,N2],(L1+L2)/2,tol);    

% Assembling vector 
[b]          = VECTORcomp_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);

% System solving
GAMMA        = MATRIX\b; 

% Plotting gamma distribution for the 2 wings
GAMMAplot(GAMMA(1:N1*2*M1)    ,M1,N1);
GAMMAplot(GAMMA(N1*2*M1+1:end),M2,N2);

% Surfaces computation
S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

% Cl-alpha and Cl-Cd computation
CLAplot_multi(MATRIX,PANELwing,beta,[AOA1,AOA2],[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],alpha_vec,"yes",4);

% AERODYNAMIC PROPERTIES FOR UNDISTURBED WINGS
% WING STUDY %
    flag = "noplot";

    % Panel creation function 
    [PANELwing]       = PANELING(delta1,lambda1,AOA1,root1,taper1,L1,M1,N1,flag,[0,0,0]);

    % System matrix generation
    % Setting tolerance to avoid singular MATRIX 
    toll              = 1e-4;
    [MATRIX]          = BIOTSAVART(PANELwing,AOA1,M1,N1,L1,toll);

    % Cl-alpha and Cd-alpha computation --2nd wing
    [Cl_vec1,Cd_vec1] = CLAplot(MATRIX,PANELwing,beta,lambda1,AOA1,M1,N1,S1,alpha_vec,"no");

% TAIL STUDY %
    flag = "noplot";

    % Panel creation function 
    [PANELwing]       = PANELING(delta2,lambda2,AOA2,root2,taper2,L2,M2,N2,flag,[0,0,0]);

    % System matrix generation
    % Setting tolerance to avoid singular MATRIX 
    toll              = 1e-4;
    [MATRIX]          = BIOTSAVART(PANELwing,AOA2,M2,N2,L2,toll);

    % Cl-alpha and Cd-alpha computation --2nd wing
    [Cl_vec2,Cd_vec2] = CLAplot(MATRIX,PANELwing,beta,lambda2,AOA2,M2,N2,S2,alpha_vec,"no");

% PLOTTING DATA for comparing the different plots -- [WING+TAIL; WING; TAIL]
figure (4)
hold on

subplot(2,2,1)
hold on
plot(AOA1+alpha_vec,Cl_vec1,'or','LineWidth',5);
legend("WING + TAIL","WING ALONE",'location','best')

subplot(2,2,2)
hold on
plot(Cd_vec1,Cl_vec1,'or','LineWidth',5);
legend("WING + TAIL","WING ALONE",'location','best')

subplot(2,2,3)
hold on
plot(AOA2+alpha_vec,Cl_vec2,'or','LineWidth',5);
legend("WING + TAIL","TAIL ALONE",'location','best')

subplot(2,2,4)
hold on
plot(Cd_vec2,Cl_vec2,'or','LineWidth',5);
legend("WING + TAIL","TAIN ALONE",'location','best')

exportgraphics(figure (3),"figures\Geometry.png",'Resolution',500);

toc

% Removing path
flpath = pwd;
rmpath(append(flpath,'/main/'));

%% ground effect
close all
clear
clc

% Initializing path
flpath = pwd;
addpath(append(flpath,'/main/'));

% plot font settings
fSize = 9;
fSizeLeg = fSize * 0.8;
fSizeTit = fSize * 1.1;

tic

% AERODYNAMIC ANGLES
alpha   = 0;
beta    = 5;
AOA_vec = -1:0.2:1;
% AOA_vec = [-2,-1,0,1,2,3,4];
unit_chord = 1;

% 1ST WING GEOMETRY
delta1    = 5;
lambda1   = 25;
root1     = 6;
L1        = 34;
taper1    = 1;
%AOA1    = 2.5;
h_s=[0.4,0.6,0.8,1];        % Quota h dal suolo [m]
transl1 = zeros(length(h_s),3);
transl1(:,3) = h_s';
M1      = 7;
N1      = 3;

% 2ND WING GEOMETRY (metodo dell'immagine)
delta2  = -delta1;
lambda2 = lambda1;
root2   = root1;
L2      = L1;
taper2  = taper1;
transl2 = - transl1;
M2      = M1;
N2      = N1;

[r,c]=size(transl1);
Cl_vec1=zeros(length(AOA_vec),length(h_s));
Cd_vec1=zeros(length(AOA_vec),length(h_s));
Cl_vec2=zeros(length(AOA_vec),length(h_s));
Cd_vec2=zeros(length(AOA_vec),length(h_s));

for j=1:r
    flag = "noplot";
    i=0;
    for AOA = AOA_vec
        i=i+1;
        [PANELwing1] = PANELING(delta1,lambda1, AOA,root1,taper1,L1,M1,N1,flag,transl1(j,:));
        [PANELwing2] = PANELING(delta2,lambda2,-AOA,root2,taper2,L2,M2,N2,flag,transl2(j,:));

        % Assemblying MATRIX
        tol          = 1e-4;
        PANELwing    = [PANELwing1,PANELwing2];
        [MATRIX]     = BIOTSAVART_multi(PANELwing,[AOA,-AOA],[M1,M2],[N1,N2],transl1(j,3)*3/4,tol);

        % Assembling vector
        [b]          = VECTORcomp_multi(PANELwing,alpha,beta,[M1,M2],[N1,N2]);

        % Surfaces computations
        S1           = (root1 + root1/taper1) * L1*cos(lambda1/180*pi);
        S2           = (root2 + root2/taper2) * L2*cos(lambda2/180*pi);

        % Cl-alpha computation
        [Cl_vec1(i,j),Cd_vec1(i,j),Cl_vec2(i,j),Cd_vec2(i,j)] = CLAplot_multi(MATRIX,PANELwing,beta,[AOA,-AOA],[lambda1,lambda2],[M1,M2],[N1,N2],[S1,S2],0,"noplot");
    end
end

% Plot Cl-alpha (parametro altezza)
figure_GroundEffectAlpha = figure (1);  hold on;
cmap_hfloor = winter(length(h_s));
for i=1:length(h_s)
    plot(AOA_vec,Cl_vec1(:,i),'o-','LineWidth',1,'Color',cmap_hfloor(i,:));
    axis padded
    grid on
    hold on
end

% Colorbar definition
colormap('winter');
clim([h_s(1),h_s(end)])
cbar_hloor = colorbar('Ticks',h_s);
cbar_hloor.Label.String = '$h/c$ [-]';
cbar_hloor.Label.Interpreter = 'latex';

% Title and legend definition
fontsize(7,"points");
title('$C_l$ vs $\alpha$', 'fontsize', fSizeTit,'interpreter', 'latex')

% Axis definition
xticks(AOA_vec)
xlabel('$\alpha$ [$^{\circ}$]','fontsize', fSize,'interpreter', 'latex')
ylabel('$C_{l}$','fontsize', fSize,'interpreter', 'latex')

% Plot Cl-altezza (parametro alpha)
figure_GroundEffectHfloor = figure (2);  hold on;
cmap_alpha = winter(length(AOA_vec));     % colormap definition

for i=1:length(AOA_vec)
    plot(h_s/unit_chord,Cl_vec1(i,:), 'o-', 'LineWidth', 1, 'Color', cmap_alpha(i,:))
    axis padded
    grid on
    hold on
end

% Colorbar definition
colormap('winter');
clim([AOA_vec(1),AOA_vec(end)])
cbar_alpha = colorbar('Ticks',AOA_vec);
cbar_alpha.Label.String = '$\alpha$ [$^{\circ}$]';
cbar_alpha.Label.Interpreter = 'latex';

% Title definition
fontsize(7,"points");
title('$C_l$ vs h/c','fontsize', fSizeTit, 'interpreter', 'latex')

% Axis definition
xticks(h_s/unit_chord)
xlabel('$h/c$ [-]','interpreter', 'latex')
ylabel('$C_{l}$','interpreter', 'latex')
fontname(figure_GroundEffectHfloor,"Palatino Linotype")
set(figure_GroundEffectHfloor,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_GroundEffectHfloor,"figures\GroundEffectHfloor.png",'Resolution',500);

% Plot polare di resistenza
figure_DragPolar = figure (3);  hold on;
cmap_hfloor = winter(length(h_s));
for i=1:length(h_s)
    plot(Cd_vec1(:,i),Cl_vec1(:,i),'o-','LineWidth',1,'Color',cmap_hfloor(i,:));
    axis padded
    grid on
    hold on
end

% Colorbar definition
colormap('winter');
clim([h_s(1),h_s(end)])
cbar_hloor = colorbar('Ticks',h_s);
cbar_hloor.Label.String = '$h/c$ [-]';
cbar_hloor.Label.Interpreter = 'latex';

% Title and legend definition
fontsize(7,"points");
title('$Drag$ $Polar$','fontsize', fSizeTit,'interpreter', 'latex')

% Axis definition
%xticks(Cd_vec1)
xlabel('$C_{d}$','fontsize', fSize,'interpreter', 'latex')
ylabel('$C_{l}$','fontsize', fSize,'interpreter', 'latex')

% Plot polare di resistenza con zoom
figure_DragPolarZoom = figure (4);  hold on;
cmap_hfloor = winter(length(h_s));
for i=1:length(h_s)
    plot(Cd_vec1(:,i),Cl_vec1(:,i),'o-','LineWidth',1,'Color',cmap_hfloor(i,:));
    axis padded
    grid on
    hold on
end

% Colorbar definition
colormap('winter');
clim([h_s(1),h_s(end)])
cbar_hloor = colorbar('Ticks',h_s);
cbar_hloor.Label.String = '$h/c$ [-]';
cbar_hloor.Label.Interpreter = 'latex';

% Title and legend definition
fontsize(7,"points");
title('$Drag$ $Polar$ $zoom$', 'fontsize', fSizeTit,'interpreter', 'latex')

% Axis definition
%xticks(Cd_vec1)
xlabel('$C_{d}$','fontsize', fSize,'interpreter', 'latex')
xlim([0 2])
ylabel('$C_{l}$','fontsize', fSize,'interpreter', 'latex')
ylim([0.05 0.1])

% AERODYNAMIC PROPERTIES FOR UNDISTURBED WINGS
% WING STUDY %
flag = "noplot";

% Panel creation function
[PANELwing]       = PANELING(delta1,lambda1,0,root1,taper1,L1,M1,N1,flag,[0,0,0]);

% System matrix generation
% Setting tolerance to avoid singular MATRIX
toll              = 1e-4;
[MATRIX]          = BIOTSAVART(PANELwing,0,M1,N1,L1,toll);

% Cl-alpha and Cd-alpha computation - 1st wing
[Cl_vec1u,Cd_vec1u] = CLAplot(MATRIX,PANELwing,beta,lambda1,0,M1,N1,S1,AOA_vec,"no");

% Adding UNDISTURBED WING data
figure (1)
hold on
und_plot1 = plot(AOA_vec,Cl_vec1u,'-r','LineWidth',1.5);
legend(und_plot1,'No ground','fontsize', fSizeLeg, 'interpreter', 'latex', 'Location','northwest')
grid on
axis padded
fontname(figure_GroundEffectAlpha,"Palatino Linotype")
box on
set(figure_GroundEffectAlpha,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_GroundEffectAlpha,"figures\GroundEffectAlpha.png",'Resolution',500);


figure (3)
hold on
und_plot3 = plot(Cd_vec1u,Cl_vec1u,'-r','LineWidth',1.5);
legend(und_plot3,'No ground','fontsize', fSizeLeg, 'interpreter', 'latex', 'Location','northeast')
grid on
axis padded
fontname(figure_DragPolar,"Palatino Linotype")
box on
set(figure_DragPolar,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_DragPolar,"figures\DragPolar.png",'Resolution',500);


figure (4)
hold on
und_plot4 = plot(Cd_vec1u,Cl_vec1u,'-r','LineWidth',1.5);
legend(und_plot4,'No ground','fontsize', fSizeLeg, 'interpreter', 'latex', 'Location','southeast')
grid on
axis padded
xlim([-0.5e-4 2.5e-4])    
ylim([0 0.1])
fontname(figure_DragPolarZoom,"Palatino Linotype")
box on
set(figure_DragPolarZoom,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_DragPolarZoom,"figures\DragPolarZoom.png",'Resolution',500);

toc

% Removing path
flpath = pwd;
rmpath(append(flpath,'/main/'));
