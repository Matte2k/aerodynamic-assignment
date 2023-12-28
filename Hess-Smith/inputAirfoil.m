% Script che inizializza e definisce le airfoil oggetto dello studio
% parametrico attraverso Hess-Smith
%
%   script impiegato in:
%   -  "mat_functions\ParamHessSmith.m
%   -  "main.m

NCorpi=1;   % Numero di corpi considerati

%% Inizializzazione delle variabili

CodiceProfilo = cell(NCorpi, 1);
Chord = zeros(NCorpi, 1);
NPannelli = zeros(NCorpi, 1);
LE_X_Position = zeros(NCorpi, 1);
LE_Y_Position = zeros(NCorpi, 1);
CalettAng_deg = zeros(NCorpi, 1);


%% Corpo 1

CodiceProfilo{1} = '0012';
Chord(1) = 1;
NPannelli(1) = 101;
LE_X_Position(1) = 0;
LE_Y_Position(1) = 0;
CalettAng_deg(1) = 0;


%% Corpo 2 - NOT IMPLEMENTED YET
% 
% if NCorpi==2
%     CodiceProfilo{2} = '23012';
%     Chord(2) = 1;
%     NPannelli(2) = 101;
%     LE_X_Position(2) = 2;
%     LE_Y_Position(2) = 3;
%     CalettAng_deg(2) = 0;
% end

