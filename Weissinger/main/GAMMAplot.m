function GAMMAplot(GAMMA,M,N)
% This function plots GAMMA ditribution spanwise
%
% INPUT: 
%   GAMMA : circulation distribution
%   M     : spanwise # of discretization points for the semiwing
%   N     : chordwise # of discretization points

GAMMA1   = zeros(N*M,1);
GAMMA2   = zeros(N*M,1);
GAMMAvec = zeros(N,2*M); 
    
for i=0:N-1
    for j=1:M
       GAMMA1(j) = GAMMA(i*M + j);
    end

    for j=1:M
       GAMMA2(j) = GAMMA(i*M + j + N*M); 
    end
    
    for j=1:M
        GAMMAvec(i+1,j) = GAMMA1(M+1-j);
    end
    
    for j=1:M
        GAMMAvec(i+1,j+M) = GAMMA2(j);
    end
end

% spanwise discretization
x = linspace(-1,1,2*M);
X = zeros(N,2*M);
for i=1:N
   X(i,:) = x; 
end

figure
hold on
for i=1:N
    plot(X(i,:),GAMMAvec(i,:),'-o','LineWidth',1.5);
end
title('$\Gamma  \ distribution$','Interpreter','latex');
xlabel('$\%SPAN$','Interpreter','latex');
ylabel('$\Gamma$','Interpreter','latex');
TEXT = "$chord \ position \ N " + string(1:N) + "$";
legend(TEXT,'Interpreter','latex');
grid on
axis padded
end