clc
close all
clear 

addpath mat_functions
run inputAirfoil.m

% plot font settings
fSize = 9;
fSizeLeg = fSize * 0.8;
fSizeTit = fSize * 1.1;


%% Plot andamento Cp-corda su un singolo profilo

U_inf = 1;
alpha_Cp = 2 ;
calettAng = 0;
EffettoSuolo = false;

dorsoInterv_end = floor(NPannelli/2);
ventreInterv_start = ceil(NPannelli/2);

% Xfoil output
[x,Cp_xfoil] = createCp(CodiceProfilo{1},NPannelli(1),Chord(1),alpha_Cp);

% Hess-Smith output
[~,Cp_pAlpha,Centro] = ParamHessSmith(U_inf, alpha_Cp, NCorpi, calettAng, EffettoSuolo);

% Comparative plot
figure_CpCompare = figure;  hold on;
    plot(Centro{1}((1:dorsoInterv_end),1), -Cp_pAlpha{1}((1:dorsoInterv_end),1), '-', 'LineWidth', 1 , 'Color', 'b')
    plot(Centro{1}((ventreInterv_start:NPannelli), 1), -Cp_pAlpha{1}((ventreInterv_start:NPannelli), 1), '-', 'LineWidth', 1 , 'Color', 'r')
    plot(Centro{1}(1:NPannelli,1), -Cp_xfoil(1:NPannelli,1), 'o', 'MarkerSize', 4 , 'Color', 'k')
    
    % title and legend definition
    fontsize(7,"points");
    legend('HS method dorso', 'HS method ventre','Xfoil', 'fontsize', fSizeLeg, 'interpreter', 'latex')
    t_string = ['$C_P$ at $\alpha=$ ',num2str(alpha_Cp), '$^{\circ}$'];
    title(t_string,'fontsize', fSizeTit, 'fontweight', 'bold', 'interpreter', 'latex')
    
    % axis definition
    ylabel('$-C_{P}$ [-]','fontsize', fSize, 'interpreter', 'latex')
    xlabel('$x/c$ [-]','fontsize', fSize, 'interpreter', 'latex')
    
    % layout definition
    axis padded
    grid on
    box on
    fontname(figure_CpCompare,"Palatino Linotype")
    set(gcf,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout    
exportgraphics(figure_CpCompare,'figures\plot_CpCompare.png','Resolution',500);
   

%% Plot andamento Cl-alpha su un singolo profilo

U_inf = 1;
alpha_Cl = (-2:1:2)';
calettAng = 0;
EffettoSuolo = false;

% Xfoil output
alphaMin = min(alpha_Cl);
alphaMax = max(alpha_Cl);
step = norm(alpha_Cl(1)-alpha_Cl(2));
[alpha_xfoil,ClAlpha_xfoil] = createClAlpha(CodiceProfilo{1},NPannelli(1),alphaMin,alphaMax,step);

% Hess-Smith output
Cl_pAlpha = cell(length(alpha_Cl),1);
Cl_pAlpha_vector = zeros(length(alpha_Cl),1);
for i=1:length(alpha_Cl)
    [Cl_pAlpha{i},~,~] = ParamHessSmith(U_inf, alpha_Cl(i), NCorpi, calettAng, EffettoSuolo);
    Cl_local = Cl_pAlpha{i};
    Cl_pAlpha_vector(i,1) = Cl_local{1};
end

% Comparative Plot
figure_ClAlphaCompare = figure;  hold on;
    plot(alpha_Cl(:,1), Cl_pAlpha_vector(:,1), '-', 'LineWidth', 1 , 'Color', 'b')
    plot(alpha_xfoil,ClAlpha_xfoil, 'o', 'MarkerSize', 5 , 'Color', 'r')
    
    % title and legend definition
    fontsize(7,"points");
    legend('HS method ', 'Xfoil', 'fontsize', fSizeLeg, 'interpreter', 'latex', 'Location','northwest')
    t_string = '$C_{l/\alpha}$';
    title(t_string,'fontsize', fSizeTit, 'fontweight', 'bold', 'interpreter', 'latex')

    % axis definition
    xlabel('$\alpha$ [$^{\circ}$]','fontsize', fSize, 'interpreter', 'latex')
    ylabel('$C_{l}$ [-]','fontsize', fSize, 'interpreter', 'latex')

    % layout definition
    axis square
    axis padded
    grid on
    box on
    fontname(figure_ClAlphaCompare,"Palatino Linotype")
    set(gcf,'units','centimeters','position',[0,0,7,7]);   % setted for the report layout
exportgraphics(figure_ClAlphaCompare,'figures\plot_ClAlphaCompare.png','Resolution',500);
   

%% Plot andamento Cl-altezza-alpha su un singolo profilo

% Condizione 
U_inf = 1;
alpha = 0;
calettAng = (-2:2)';

% Caso senza effetto suolo
EffettoSuolo = false;
Cl_matrix_NoGround = zeros(length(calettAng),1);

for i=1:length(calettAng)
    [Cl,Cp,Centro] = ParamHessSmith(U_inf, alpha, NCorpi, calettAng(i), EffettoSuolo);
    Cl_matrix_NoGround(i) = Cl{1,1};
end

% Caso effetto suolo
EffettoSuolo = true;
altezzaSuolo = (0.5:0.2:1.5)';
Cl_matrix_Ground = zeros(length(calettAng),length(altezzaSuolo));

for i=1:length(calettAng)
    for j=1:length(altezzaSuolo)
        [Cl,Cp,Centro] = ParamHessSmith(U_inf, alpha, NCorpi, calettAng(i), EffettoSuolo, altezzaSuolo(j));
        Cl_matrix_Ground(i,j) = Cl{1,1};
    end
end


% Plot Cl-altezza (parametro alpha)
figure_GroundEffectHfloor = figure;  hold on;
    cmap_alpha = winter(length(calettAng));     % colormap definition

    for i=1:length(calettAng)
        plot(altezzaSuolo(:,1),Cl_matrix_Ground(i,:), 'o-', 'LineWidth', 1, 'Color', cmap_alpha(i,:) )
        hold on
    end

    % colorbar definition
    colormap('winter');
    clim([calettAng(1),calettAng(end)])
    cbar_alpha = colorbar('Ticks',calettAng);
    cbar_alpha.Label.String = '$\alpha$ [$^{\circ}$]';
    cbar_alpha.Label.Interpreter = 'latex';
    
    % title definition
    fontsize(7,"points");
    t_string = '$C_l$ vs h';
    title(t_string,'fontsize', fSizeTit, 'fontweight', 'bold', 'interpreter', 'latex')

    % axis definition
    xticks(altezzaSuolo(:,1))
    xlabel('$h$ [m]','fontsize', fSize, 'interpreter', 'latex')
    ylabel('$C_{l}$','fontsize', fSize, 'interpreter', 'latex')

    % layout settings
    axis padded
    grid on
    box on
    fontname(figure_GroundEffectHfloor,"Palatino Linotype")
    set(gcf,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_GroundEffectHfloor,'figures\plot_GroundEffectHfloor.png','Resolution',500);


% Plot Cl-alpha (parametro altezza)
figure_GroundEffectAlpha = figure;  hold on;
    cmap_hfloor = winter(length(altezzaSuolo));     % colormap definition

    xf_plot = plot(calettAng(:,1), Cl_matrix_NoGround(:,1), '--');

    for i=1:length(altezzaSuolo)
        plot(calettAng(:,1),Cl_matrix_Ground(:,i), 'o-', 'LineWidth', 1, 'Color', cmap_hfloor(i,:) )
        hold on
    end
    
    % colorbar definition
    colormap('winter');
    clim([altezzaSuolo(1),altezzaSuolo(end)])
    cbar_hloor = colorbar('Ticks',altezzaSuolo);
    cbar_hloor.Label.String = '$h$ [m]';
    cbar_hloor.Label.Interpreter = 'latex';
    
    % title and legend definition
    fontsize(7,"points");
    legend(xf_plot,'xfoil', 'fontsize', fSizeLeg, 'interpreter', 'latex', 'Location','northwest')
    t_string = '$C_l$ vs $\alpha$';
    title(t_string,'fontsize', fSizeTit, 'fontweight', 'bold', 'interpreter', 'latex')
    
    % axis definition
    xticks(calettAng(:,1))
    xlabel('$h$ [m]','fontsize', fSize, 'interpreter', 'latex')
    ylabel('$C_{l}$','fontsize', fSize, 'interpreter', 'latex')
    
    % layout settings
    axis padded
    grid on
    box on
    fontname(figure_GroundEffectAlpha,"Palatino Linotype")
    set(figure_GroundEffectAlpha,'units','centimeters','position',[0,0,10,7]);   % setted for the report layout
exportgraphics(figure_GroundEffectAlpha,'figures\plot_GroundEffectAlpha.png','Resolution',500);


%% Plot andamento Cp-altezza-alpha su un singolo profilo
 
U_inf = 1;
alpha = 0;
NCorpi = 1;
calettAng = [-4 4 8]';
altezzaSuolo = [0.5 1]';

cmap_cp = hsv(length(altezzaSuolo)+1);     % colormap definition
clString = cell(length(altezzaSuolo)+1,1);      % cl string cell definition

tiledlayout(1,length(calettAng));        % multiplot initialization

for i=1:length(calettAng)
    
    nexttile
    hold on;
    
    for j=1:length(altezzaSuolo)
        % Caso effetto suolo
        [Cl_parGroud,Cp_paraGround,Centro] = ParamHessSmith(U_inf, alpha, NCorpi, calettAng(i), true, altezzaSuolo(j));
        plot(Centro{1}(:,1), -Cp_paraGround{1}(:,1), '-', 'LineWidth', 1, 'Color', cmap_cp(j,:))
        clString{j} = ['$C_l$ = ', num2str(Cl_parGroud{1})];
    end

    % Caso senza effetto suolo
    [Cl_parNoGroud,Cp_paraNoGround,Centro] = ParamHessSmith(U_inf, alpha, NCorpi, calettAng(i), false);
    plot(Centro{1}(:,1), -Cp_paraNoGround{1}(:,1), '-', 'LineWidth', 1, 'Color', cmap_cp(end,:))
    clString{end} = ['$C_l$ = ', num2str(Cl_parNoGroud{1})]; 

    % title and legend definition
    legend(clString, 'fontsize', fSizeLeg, 'interpreter', 'latex')
    t_string = ['$\alpha=$ ',num2str(calettAng(i)), '$^{\circ}$'];
    title(t_string,'fontsize', fSizeTit, 'fontweight', 'bold', 'interpreter', 'latex')

    % axis definition
    ylabel('$-C_{P}$ [-]','fontsize', fSize, 'interpreter', 'latex')
    ylim([-1.2 4.5])                                                    % To set manually for the moment
    xlabel('$x/c$ [-]','fontsize', fSize, 'interpreter', 'latex')
    grid on
    box on
    fontname("Palatino Linotype")
end

    % colorbar definition
    alphaColorBar = [altezzaSuolo; 1.5];
    colormap(cmap_cp);
    clim([alphaColorBar(1), alphaColorBar(end)])
    cbar_CpFloor = colorbar('Ticks',alphaColorBar,'TickLabelsMode', 'manual');
    cbar_CpFloor.TickLabels = {'0.5','1','$\infty$'};
    cbar_CpFloor.TickLabelInterpreter = 'latex'; 
    cbar_CpFloor.Label.String = '$h$ [m]';
    cbar_CpFloor.Label.Interpreter = 'latex';

set(gcf,'units','centimeters','position',[0,0,20,7]);   % setted for the report layout    
exportgraphics(gcf,'figures\plot_CpCompare.png','Resolution',500);
