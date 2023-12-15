clc
close all
clear 

addpath mat_functions
run Airfoil_input.m

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
figure
plot(Centro{1}((1:dorsoInterv_end),1), -Cp_pAlpha{1}((1:dorsoInterv_end),1), '-')
hold on
plot(Centro{1}((ventreInterv_start:NPannelli), 1), -Cp_pAlpha{1}((ventreInterv_start:NPannelli), 1), '-')
plot(Centro{1}(1:NPannelli,1), -Cp_xfoil(1:NPannelli,1), '*')
title('$C_P$', 'interpreter', 'latex')
%legend(legend_string, 'interpreter', 'latex')
ylabel('$C_{P}$', 'interpreter', 'latex')
xlabel('$x/c$', 'interpreter', 'latex')
axis padded
grid on



%% Plot andamento Cl-alpha su un singolo profilo

U_inf = 1;
alpha_Cl = (-2:1:2)' ;
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
figure
plot(alpha_Cl(:,1), Cl_pAlpha_vector(:,1), '-')
hold on
plot(alpha_xfoil,ClAlpha_xfoil, '*')
title('$C_{l/\alpha}$', 'interpreter', 'latex')
xlabel('\alpha')
ylabel('$C_{l}$', 'interpreter', 'latex')
axis square
grid on


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
altezzaSuolo = (0:1)';
Cl_matrix_Ground = zeros(length(calettAng),length(altezzaSuolo));

for i=1:length(calettAng)
    for j=1:length(altezzaSuolo)
        [Cl,Cp,Centro] = ParamHessSmith(U_inf, alpha, NCorpi, calettAng(i), EffettoSuolo, altezzaSuolo(j));
        Cl_matrix_Ground(i,j) = Cl{1,1};
    end
end


% Plot Cl-altezza (parametro alpha)



% Plot Cl-altezza (parametro alpha)



% Plot Cl-altezza-alpha 3D



