function [x,y] = createClAlpha(Profilo,NPannelli,alphaMin,alphaMax,step)
% Sfruttiamo XFoil per la creazione del profilo, sperando che sia
% all'interno del database


    fileID = fopen('XFoilClAlpha.txt','w');
    fprintf(fileID, ['naca ' ' '  Profilo, '\n\n']);
    fprintf(fileID,'pane\n\n');

    fprintf(fileID,'gdes\n');
    fprintf(fileID,'tgap 0 0 \n');
    
    fprintf(fileID,'exec \n\n\n');

    fprintf(fileID,'ppar\n');
    fprintf(fileID, ['n ' ' ' num2str(NPannelli+1)  '\n\n\n']);
    % NPannelli+1 perchè Xfoil da coordinate bordi pannelli iniziando e
    % terminando su TE

    filename = strcat('NACA', Profilo, '_ClAlpha.dat');

    % Elimino precedente file contente polari se esiste
    if exist(filename,'file')==2
        delete(filename);
    end

    fprintf(fileID,'oper\n');
    fprintf(fileID,'pacc\n');
    fprintf(fileID, [' ' filename '\n\n']);

    fprintf(fileID, 'aseq\n');
    fprintf(fileID, [' ' num2str(alphaMin) '\n']);
    fprintf(fileID, [' ' num2str(alphaMax) '\n']);
    fprintf(fileID, [' ' num2str(step) '\n']);
    fprintf(fileID,'pacc\n\n\n');

    
    fprintf(fileID,'quit \n\n');
    fclose(fileID);

    % Str2Exec = strcat("xfoil < XFoilInput.txt ");
    % Riscrittura che permette uscita da xfoil automatica in CommandWindow
    Str2Exec = strcat("xfoil < XFoilClAlpha.txt > /dev/null 2>&1");

    system(Str2Exec);

    ClAplha_xfoil = importXfoilClAlpha(filename,13);
    x = ClAplha_xfoil.alpha;
    y = ClAplha_xfoil.CL;
    

end