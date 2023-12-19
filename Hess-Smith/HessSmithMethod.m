%% Hess Smith Method
% Implemented for multi airfoil configuration and gound effect
% Current case:     NACA xxxx

clc
close all
clear 

addpath mat_functions

%% Input

U_inf0 = 1;   % Velocità all'infinito [m/s]
alpha = 0;   % Angolo di incidenza flusso [°]
U_inf_x = U_inf0 * cos(deg2rad(alpha));
U_inf_y = U_inf0 * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

NCorpi = 1;             % Numero di corpi da analizzare
EffettoSuolo = true;    % Effetto suolo presente o meno

CodiceProfilo = cell(NCorpi, 1);
Chord = zeros(NCorpi, 1);
NPannelli = zeros(NCorpi, 1);
LE_X_Position = zeros(NCorpi, 1);
LE_Y_Position = zeros(NCorpi, 1);
CalettAng_deg = zeros(NCorpi, 1);

% Corpo 1
CodiceProfilo{1} = '0012';
Chord(1) = 1;
NPannelli(1) = 101;
LE_X_Position(1) = 0;
LE_Y_Position(1) = 2;   % Aggiungi qui altezza da suolo
CalettAng_deg(1) = 10;

% % Corpo 2
% CodiceProfilo{2} = '23012';
% Chord(2) = 1;
% NPannelli(2) = 101;
% LE_X_Position(2) = 2;
% LE_Y_Position(2) = 3;
% CalettAng_deg(2) = 25;


%% Creazione profilo

Corpi = cell(NCorpi, 1);
% Numero profilo:
for i=1:NCorpi
    [x,y]=createProfile(CodiceProfilo{i},NPannelli(i),Chord(i));
    Corpi{i}.x = x;
    Corpi{i}.y = y; 
end


%% Calettamento profilo

% Rotazione profili e plot profili ruotati
figure
for i=1:NCorpi
    
    % Coordinate del centro di rotazione
    xc_rotC = 0.25 * Chord(i) ;  % 1/4*chord
    yc_rotC = 0 ;

    % Shifto a centro di rotazione
    x_shift = Corpi{i}.x - xc_rotC;  % shifted data
    y_shift = Corpi{i}.y - yc_rotC;

    % Ruoto le coordinate dell'angolo di calettamentoo
    CalettAng_rad = deg2rad(CalettAng_deg(i));
    x_shiftRot =  x_shift .* cos(CalettAng_rad) + y_shift .* sin(CalettAng_rad); 
    y_shiftRot = -x_shift .* sin(CalettAng_rad) + y_shift .* cos(CalettAng_rad);
        
    % Ritorno alle coordinate reali sottraendo lo "shift" fatto
    Corpi{i}.x = x_shiftRot + xc_rotC;
    Corpi{i}.y = y_shiftRot + yc_rotC;
    
    % Plot della configurazione finale dei profili
    plot(Corpi{i}.x(:, 1), Corpi{i}.y(:, 1))
    hold on
    axis equal
end
title('Configurazione finale ruotata', 'interpreter', 'latex')
axis equal
hold off


%% Traslazione profilo

% Traslazione profili e plot profili traslati
figure
for i=1:NCorpi
    
    Corpi{i}.x = Corpi{i}.x + LE_X_Position(i);
    Corpi{i}.y = Corpi{i}.y + LE_Y_Position(i); 

    plot(Corpi{i}.x(:, 1), Corpi{i}.y(:, 1))
    hold on
end   
title('Configurazione intermedia traslata', 'interpreter', 'latex')
axis equal
hold off


%% Metodo dell'immagine effetto suolo

if EffettoSuolo
 
    NCorpi_mirror = NCorpi;
    NPannelli_mirror = NPannelli;
    Chord_mirror = Chord;
    Corpi_mirror = Corpi;
    for i=1:NCorpi_mirror       
        Corpi_mirror{i}.y = -Corpi{i}.y;
    end
end


%% Creazione di una struttura di pannelli

Centro = cell(NCorpi, 1);
Normale = cell(NCorpi, 1);
Tangente = cell(NCorpi, 1);
Estremo_1 = cell(NCorpi, 1);
Estremo_2 = cell(NCorpi, 1);
alpha = cell(NCorpi, 1);
lunghezza = cell(NCorpi, 1);
L2G_TransfMatrix = cell(NCorpi, 1);
G2L_TransfMatrix = cell(NCorpi, 1);

% Creazione pannellizzazione dei corpi
for i = 1:NCorpi
    [Centro{i}, Normale{i}, Tangente{i}, Estremo_1{i}, Estremo_2{i}, alpha{i}, lunghezza{i}, L2G_TransfMatrix{i}, G2L_TransfMatrix{i}] = CreaStrutturaPannelli(Corpi{i});
end        

% Creazione pannelli dei corpi caso effetto suolo
if EffettoSuolo

    Centro_mirror = cell(NCorpi, 1);
    Normale_mirror = cell(NCorpi, 1);
    Tangente_mirror = cell(NCorpi, 1);
    Estremo_1_mirror = cell(NCorpi, 1);
    Estremo_2_mirror = cell(NCorpi, 1);
    alpha_mirror = cell(NCorpi, 1);
    lunghezza_mirror = cell(NCorpi, 1);
    L2G_TransfMatrix_mirror = cell(NCorpi, 1);
    G2L_TransfMatrix_mirror = cell(NCorpi, 1);

    for i = 1:NCorpi_mirror
        [Centro_mirror{i}, Normale_mirror{i}, Tangente_mirror{i}, Estremo_1_mirror{i}, Estremo_2_mirror{i}, alpha_mirror{i}, lunghezza_mirror{i}, L2G_TransfMatrix_mirror{i}, G2L_TransfMatrix_mirror{i}] = CreaStrutturaPannelli(Corpi_mirror{i});
    end
end

% Plot di tutte le singole pannellizzazioni
title_string = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi    
    figure;
    plot(Centro{Corpo_i}(:, 1), Centro{Corpo_i}(:, 2), 'o-')
    hold on   
    title_string{Corpo_i} = strcat("Pannellizzazione corpo ", num2str(Corpo_i));
    title(title_string{Corpo_i}, 'interpreter', 'latex')
    axis equal
    hold off
end


%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori
NCols = sum(NPannelli) + NCorpi;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);


%% Creazione della matrice A

% Definizione Matrice A11
indexStart_riga = 0;
for Corpo_i = 1:NCorpi
    for i = 1:NPannelli(Corpo_i)
        
        index_i = indexStart_riga + i; % riga
        indexStart_colonna = 0;
    
        Centro_qui = Centro{Corpo_i}(i, :)';
        Normale_qui = Normale{Corpo_i}(i, :)';
    
        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)
    
                index_j = indexStart_colonna + j;  % colonna
    
                Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';
                Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';    
                L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));
    
                Ujs = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                Ujv = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                
                % Caso effetto suolo
                if EffettoSuolo

                    Estremo_1_qui_mirror = Estremo_1_mirror{Corpo_j}(j, :)';
                    Estremo_2_qui_mirror = Estremo_2_mirror{Corpo_j}(j, :)';
    
                    L2G_TransfMatrix_qui_mirror = squeeze(L2G_TransfMatrix_mirror{Corpo_j}(j, :, :));
                    G2L_TransfMatrix_qui_mirror = squeeze(G2L_TransfMatrix_mirror{Corpo_j}(j, :, :));

                    Ujs_mirror = ViSorgente(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                    Ujv_mirror = ViVortice(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                else
                    Ujs_mirror = 0;
                    Ujv_mirror = 0;
                end

                matriceA(index_i, index_j) = dot( (Ujs + Ujs_mirror), Normale_qui);
                matriceA(index_i, sum(NPannelli)+Corpo_j) = matriceA(index_i, sum(NPannelli)+Corpo_j) + dot( (Ujv - Ujv_mirror), Normale_qui);
    
            end
    
            indexStart_colonna = indexStart_colonna + NPannelli(Corpo_j);            
        end
    end

    indexStart_riga = indexStart_riga + NPannelli(Corpo_i);
end

% Restante parte della Matrice A
for Corpo_i = 1:NCorpi
    
    Centro_Start = Centro{Corpo_i}(1, :)';
    Tangente_Start = Tangente{Corpo_i}(1, :)';     
    Centro_End = Centro{Corpo_i}(end, :)';
    Tangente_End = Tangente{Corpo_i}(end, :)'; 
    
    indexStart_colonna = 0;
        
    for Corpo_j = 1:NCorpi
        b = 0;
        for j = 1:NPannelli(Corpo_j)
    
            index_j = indexStart_colonna + j;

            Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';                
            Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

            Ujs_start = ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
            Ujv_start = ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

            Ujs_end = ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
            Ujv_end = ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

            % Caso effetto suolo
            if EffettoSuolo

                Estremo_1_qui_mirror = Estremo_1_mirror{Corpo_j}(j, :)';
                Estremo_2_qui_mirror = Estremo_2_mirror{Corpo_j}(j, :)';
                L2G_TransfMatrix_qui_mirror = squeeze(L2G_TransfMatrix_mirror{Corpo_j}(j, :, :));
                G2L_TransfMatrix_qui_mirror = squeeze(G2L_TransfMatrix_mirror{Corpo_j}(j, :, :));

                Ujs_mirror_start = ViSorgente(Centro_Start, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                Ujv_mirror_start = ViVortice(Centro_Start, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
            
                Ujs_mirror_end = ViSorgente(Centro_End, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                Ujv_mirror_end = ViVortice(Centro_End, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);     
            else
                Ujs_mirror_start = 0;
                Ujv_mirror_start = 0;

                Ujs_mirror_end = 0;
                Ujv_mirror_end = 0;
            end

            a = dot((Ujs_start + Ujs_mirror_start), Tangente_Start);
            b = b + dot((Ujv_start - Ujv_mirror_start), Tangente_Start);

            a = a + dot((Ujs_end + Ujs_mirror_end), Tangente_End);
            b = b + dot((Ujv_end - Ujv_mirror_end), Tangente_End);

            matriceA(sum(NPannelli) + Corpo_i, index_j) = a;

        end
        
        matriceA(sum(NPannelli) + Corpo_i, sum(NPannelli) + Corpo_j) = b;
        indexStart_colonna = indexStart_colonna + NPannelli(Corpo_j);
    end
end


%% Creazione del termine noto

% Termine noto
indexStart = 0;
for Corpo_i = 1:NCorpi
    for j = 1:NPannelli(Corpo_i)

        Normale_qui = Normale{Corpo_i}(j, :)'; 
        index = indexStart + j;
        TermineNoto(index) = - dot(U_inf, Normale_qui);
    end
    
    Tangente_1 = Tangente{Corpo_i}(1, :)';  
    Tangente_end = Tangente{Corpo_i}(end, :)';
    TermineNoto(sum(NPannelli) + Corpo_i) = - dot(U_inf, (Tangente_1 + Tangente_end));
    
    indexStart = indexStart + NPannelli(Corpo_i);

end


%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);

sigma_mia = cell(NCorpi,1);
gamma_mia = zeros(NCorpi,1);

sigma_start = 0;
sigma_end = 0;

% Definizione vettore sorgenti e vorticità
for Corpo_i=1:NCorpi
    % Sigma
    sigma_start = sigma_end + 1;
    sigma_end = sigma_start + NPannelli(Corpo_i) - 1;
    sigma_mia{Corpo_i} = Soluzione(sigma_start:sigma_end,1);
    
    % Gamma
    gamma_selector = (sum(NPannelli) + Corpo_i);
    gamma_mia(Corpo_i) = Soluzione(gamma_selector,1);
end


%% Calcolo del cp e della velocità sui pannelli

U_Pannelli = cell(NCorpi, 1);
Ut_Pannelli = cell(NCorpi, 1);
Un_Pannelli = cell(NCorpi, 1);
Cp = cell(NCorpi, 1);

% Inizializzazione vettore velocità
for Corpo_i = 1:NCorpi    
    U_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),2);
    Ut_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),1);
    Un_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),1);    
end

% Calcolo del campo di velocità
for Corpo_i = 1:NCorpi
    for i = 1:NPannelli(Corpo_i)

        U_Pannelli{Corpo_i}(i, :) = U_inf'; 
        Centro_qui = Centro{Corpo_i}(i, :)';
        Tangente_qui = Tangente{Corpo_i}(i, :)'; 
        Normale_qui = Normale{Corpo_i}(i, :)'; 
    
        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)

                Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                % Caso effetto suolo
                if EffettoSuolo
    
                    Estremo_1_qui_mirror = Estremo_1_mirror{Corpo_j}(j, :)';
                    Estremo_2_qui_mirror = Estremo_2_mirror{Corpo_j}(j, :)';
                    L2G_TransfMatrix_qui_mirror = squeeze(L2G_TransfMatrix_mirror{Corpo_j}(j, :, :));
                    G2L_TransfMatrix_qui_mirror = squeeze(G2L_TransfMatrix_mirror{Corpo_j}(j, :, :));
    
                    U_Sorgente_mirror = ViSorgente(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                    U_Vortice_mirror = ViVortice(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                
                else
                    U_Sorgente_mirror = 0;
                    U_Vortice_mirror = 0;
                end

                U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                U_Pannelli{Corpo_i}(i, :) = U_Pannelli{Corpo_i}(i, :) + ...
                                            sigma_mia{Corpo_j}(j) .* U_Sorgente' + gamma_mia(Corpo_j) .* U_Vortice' + ...
                                            sigma_mia{Corpo_j}(j) .* U_Sorgente_mirror' - gamma_mia(Corpo_j) .* U_Vortice_mirror' ;
            end
        end
        
        Ut_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Tangente_qui);
        Un_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Normale_qui);
    end
    
    Cp{Corpo_i} = 1-Ut_Pannelli{Corpo_i}.^2/norm(U_inf)^2;
end

% Calcolo del CL
Cl = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
    Cl_qui = 0;
    for i = 1:NPannelli(Corpo_i)

        Normale_qui = Normale{Corpo_i}(i, :)';
        lunghezza_qui = lunghezza{Corpo_i}(i);

        Cl_qui = Cl_qui + (-Cp{Corpo_i}(i) * (lunghezza_qui .* dot(Normale_qui, U_inf_normal)));
    end
    
    Cl{Corpo_i} = Cl_qui / Chord(Corpo_i);
end


%% Plot Cp e profili

% Plot andamento Cp su profilo
figure
legend_string = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
%   plot(Centro{Corpo_i}(:, 1) - min(Centro{Corpo_i}(:, 1)), -Cp{Corpo_i}, '-')
    plot(Centro{Corpo_i}(:, 1), -Cp{Corpo_i}, '-')
    hold on
    legend_string{Corpo_i} = strcat("Corpo ", num2str(Corpo_i));
end
title('$C_P$', 'interpreter', 'latex')
legend(legend_string, 'interpreter', 'latex')

% Plot profili analizzati
figure
legend_string = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
%   plot(Centro{Corpo_i}(:, 1) - min(Centro{Corpo_i}(:, 1)), -Cp{Corpo_i}, '-')
    plot(Corpi{Corpo_i}.x, Corpi{Corpo_i}.y, '-')
    hold on
    legend_string{Corpo_i} = strcat("Corpo ", num2str(Corpo_i));
end
axis equal
title('Profili analizzati', 'interpreter', 'latex')
legend(legend_string, 'interpreter', 'latex')


%% Render della configurazione studiata

% Attiva/disattiva render a colori della corrente
ifSaveFigures = true;

if ifSaveFigures

    xMin = 2000;
    xMax = -2000;
    yMin = 2000;
    yMax = -2000;
    
    for Corpo_i = 1:NCorpi      
        xMin = min(xMin, min(Centro{Corpo_i}(:, 1)));
        xMax = max(xMax, max(Centro{Corpo_i}(:, 1)));
        
        yMin = min(yMin, min(Centro{Corpo_i}(:, 2)));
        yMax = max(yMax, max(Centro{Corpo_i}(:, 2)));
    end
    
    xMin = xMin - 1;
    xMax = xMax + 1;
    yMin = yMin - 1;
    yMax = yMax + 1;
        
    Nx = 400;
    Ny = 400;
        
    x = linspace(xMin,xMax, Nx);
    y = linspace(yMin, yMax, Ny);

    [X, Y] = meshgrid(x, y);

    isIn = zeros(Nx, Ny);
    
    parfor i = 1:Nx
%       i
        for j = 1:Ny
            for Corpo_i = 1:NCorpi
                Boundary = [Corpi{Corpo_i}.x Corpi{Corpo_i}.y];
                if (inpolygon(X(i, j), Y(i, j), Boundary(:, 1), Boundary(:, 2)))
                    isIn(i, j) = 1;
                end
            end
        end
    end

    U_Mesh = zeros(Nx, Ny);
    V_Mesh = zeros(Nx, Ny);
    U_Mesh_Mag = zeros(Nx, Ny);
    Cp_Mesh = zeros(Nx, Ny);

    t = cputime;
    parfor PointIndex_i = 1:Nx
%       PointIndex_i
        for PointIndex_j = 1:Ny

            if(~isIn(PointIndex_i, PointIndex_j))

                U = U_inf'; 
                Centro_qui = [X(PointIndex_i, PointIndex_j); Y(PointIndex_i, PointIndex_j)];

                for Corpo_j = 1:NCorpi
                    for j = 1:NPannelli(Corpo_j)

                        Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                        Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                        L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                        G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                        % Caso effetto suolo
                        if EffettoSuolo
            
                            Estremo_1_qui_mirror = Estremo_1_mirror{Corpo_j}(j, :)';
                            Estremo_2_qui_mirror = Estremo_2_mirror{Corpo_j}(j, :)';
                            L2G_TransfMatrix_qui_mirror = squeeze(L2G_TransfMatrix_mirror{Corpo_j}(j, :, :));
                            G2L_TransfMatrix_qui_mirror = squeeze(G2L_TransfMatrix_mirror{Corpo_j}(j, :, :));
            
                            U_Sorgente_mirror = ViSorgente(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                            U_Vortice_mirror = ViVortice(Centro_qui, Estremo_1_qui_mirror, Estremo_2_qui_mirror, L2G_TransfMatrix_qui_mirror, G2L_TransfMatrix_qui_mirror);
                        
                        else
                            U_Sorgente_mirror = 0;
                            U_Vortice_mirror = 0;
                        end

                        U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                        U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                        U = U + sigma_mia{Corpo_j}(j) .* U_Sorgente' + gamma_mia(Corpo_j) .* U_Vortice' + ...
                            sigma_mia{Corpo_j}(j) .* U_Sorgente_mirror' - gamma_mia(Corpo_j) .* U_Vortice_mirror'

                    end
                end

                U_Mesh(PointIndex_i, PointIndex_j) = U(1);
                V_Mesh(PointIndex_i, PointIndex_j) = U(2);
                U_Mesh_Mag(PointIndex_i, PointIndex_j) = norm(U);
                Cp_Mesh(PointIndex_i, PointIndex_j) = 1 - norm(U) / norm(U_inf);

            else
                U_Mesh(PointIndex_i, PointIndex_j) = NaN;
                V_Mesh(PointIndex_i, PointIndex_j) = NaN;
                U_Mesh_Mag(PointIndex_i, PointIndex_j) = NaN;
                Cp_Mesh(PointIndex_i, PointIndex_j) = NaN;
            end
        end
    end

    SavingNameStart = "./figures/test_";

    UMag_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
%   legend("$U_{MAG}$", 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_UMag.eps');
    saveas(UMag_fig, SavingName);
    
    
    U_fig = figure;
    contourf(X, Y, U_Mesh, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_U.eps');
    saveas(U_fig, SavingName);
    
    
    V_fig = figure;
    contourf(X, Y, V_Mesh, 1000,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $V$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_V.eps');
    saveas(V_fig, SavingName);


    Cp_fig = figure;
    contourf(X, Y, Cp_Mesh, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Cp.eps');
    saveas(Cp_fig, SavingName);
    
    
    Streamlines_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    streamslice(X, Y, U_Mesh, V_Mesh, 10);
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x, Corpi{Corpo_j}.y, '-k')
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$ e linee di corrente', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Streamlines.eps');
    saveas(Streamlines_fig, SavingName);
    
end
