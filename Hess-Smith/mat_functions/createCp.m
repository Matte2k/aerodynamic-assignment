function [x,y] = createCp(Profilo,NPannelli,Chord,alpha)
% Sfruttiamo XFoil per la creazione del profilo, sperando che sia
% all'interno del database

    fileID = fopen('XFoilCp.txt','w');
    fprintf(fileID, ['naca ' ' '  Profilo, '\n\n']);
    fprintf(fileID,'pane\n\n');

    fprintf(fileID,'gdes\n');
    fprintf(fileID,'tgap 0 0 \n');
    
    fprintf(fileID,'exec \n\n\n');

    fprintf(fileID,'ppar\n');
    fprintf(fileID, ['n ' ' ' num2str(NPannelli+1)  '\n\n\n']);
    % NPannelli+1 perchè Xfoil da coordinate bordi pannelli iniziando e
    % terminando su TE

    filename = strcat('NACA', Profilo, '_Cp.dat');

    fprintf(fileID,'oper\n');
    fprintf(fileID, ['alfa ' ' ' num2str(alpha) ' \n']);
    fprintf(fileID, ['cpwr ' ' ' filename '\n\n\n']);
    
    
    fprintf(fileID,'quit \n\n');
    fclose(fileID);

    % Str2Exec = strcat("xfoil < XFoilInput.txt ");
    % Riscrittura che permette uscita da xfoil automatica in CommandWindow
    Str2Exec = strcat("xfoil < XFoilCp.txt > /dev/null 2>&1");

    system(Str2Exec);

    Cp_xfoil = importXfoilProfile(filename);
    
    % Prima flippa i vettori
    x = flipud(Cp_xfoil.x);
    y = flipud(Cp_xfoil.y);
    
    x = x.*Chord;

end