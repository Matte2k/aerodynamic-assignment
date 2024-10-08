function [Cl_vec1,Cd_vec1,Cl_vec2,Cd_vec2] = CLAplot_multi(MATRIX,PANELwing,beta,AOA,lambda,M,N,S,alpha_vec,flag,num)
% This function computes the Cl alpha plot of a 3D wing given geometry and
% flow conditions
% 
% INPUT: 
%   MATRIX    : system matrix -- describes the non penetration condition
%   PANELwing : PANEL class array
%   beta      : sideslip angle
%   M         : spanwise # of discretization points
%   N         : chordwise # of discretization points
%   S         : wing surface
%   alpha_vec : set AOA
%   flag      : toggling plot procedure

% initializing variables
Cl_vec1 = zeros(size(alpha_vec));
Cd_vec1 = zeros(size(alpha_vec));
Cl_vec2 = zeros(size(alpha_vec));
Cd_vec2 = zeros(size(alpha_vec));
U       = 1;
rho     = 1;

for i=1:length(alpha_vec)
        
    alpha = alpha_vec(i);
    
    % system known vector 
    [b]   = VECTORcomp_multi(PANELwing,alpha,beta,[M(1),M(2)],[N(1),N(2)]);

    % solve system
    GAMMA = MATRIX\b;

    % computing LIFT
    [~,L_vec1,Cl_vec1(i),~,~] = LIFTcomp(GAMMA(1:N(1)*2*M(1)),    PANELwing(1:N(1)*2*M(1)),    lambda(1),M(1),N(1),rho,U,S(1),"no");
    [~,L_vec2,Cl_vec2(i),~,~] = LIFTcomp(GAMMA(N(1)*2*M(1)+1:end),PANELwing(N(1)*2*M(1)+1:end),lambda(2),M(2),N(2),rho,U,S(2),"no");
    
    % computing induced velocity 
    [~,alpha_ind1]   = INDUCEDVELcomp(GAMMA(1:N(1)*2*M(1)),    PANELwing(1:N(1)*2*M(1)),    alpha,M(1),N(1),U,"no"); 
    [~,alpha_ind2]   = INDUCEDVELcomp(GAMMA(N(1)*2*M(1)+1:end),PANELwing(N(1)*2*M(1)+1:end),alpha,M(2),N(2),U,"no"); 

    % computing DRAG
    [~,~,Cd_vec1(i)] = DRAGcomp(L_vec1,-alpha_ind1,alpha,rho,U,S(1),M(1));
    [~,~,Cd_vec2(i)] = DRAGcomp(L_vec2,-alpha_ind2,alpha,rho,U,S(2),M(2));
end 

% plotting procedure
if(flag == "yes")
    if(nargin == 11)
        figure(num)
        hold on
    else
        figure
        hold on
    end

    subplot(2,2,1)
    plot(AOA(1) + alpha_vec,Cl_vec1,'k-','LineWidth',2);
    xlabel('$\alpha$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    title('WING','Interpreter','latex');
    hold on
    grid on 
    axis padded

    subplot(2,2,2)
    plot(Cd_vec1,Cl_vec1,'k-','LineWidth',2);
    xlabel('$C_D$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    title('WING','Interpreter','latex');
    hold on
    grid on 
    axis padded

    subplot(2,2,3)
    plot(AOA(2) + alpha_vec,Cl_vec2,'k-','LineWidth',2);
    xlabel('$\alpha$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    title('TAIL','Interpreter','latex');
    hold on
    grid on 
    axis padded

    subplot(2,2,4)
    plot(Cd_vec2,Cl_vec2,'k-','LineWidth',2);
    xlabel('$C_D$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    title('TAIL','Interpreter','latex');
    hold on
    grid on 
    axis padded    
end 
end