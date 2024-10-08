function [Cl_vec,Cd_vec] = CLAplot(MATRIX,PANELwing,beta,lambda,AOA,M,N,S,alpha_vec,flag,num)
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

% initializing variables
Cl_vec = zeros(size(alpha_vec));
Cd_vec = zeros(size(alpha_vec));

for i=1:length(alpha_vec)
        
    alpha = alpha_vec(i);
    
    % system known vector 
    [b]         = VECTORcomp(PANELwing,alpha,beta,M,N);

    % solve system
    GAMMA       = MATRIX\b;

    % computing LIFT
    U                       = 1;
    rho                     = 1;
    [~,L_vec,Cl_vec(i),~,~] = LIFTcomp(GAMMA,PANELwing,lambda,M,N,rho,U,S,"no");
    
    % computing induced velocity 
    [~,alpha_ind]   = INDUCEDVELcomp(GAMMA,PANELwing,alpha,M,N,U,"no"); 

    % computing DRAG
    [~,~,Cd_vec(i)] = DRAGcomp(L_vec,-alpha_ind,alpha,rho,U,S,M);
end 

if(flag == "yes")
    % plotting procedure
    if(nargin == 11)
        figure(num)
        hold on
    else
        figure
        hold on
    end
       
    figure (6)
    plot(AOA + alpha_vec,Cl_vec,'k-','LineWidth',2);
    xlabel('$\alpha$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    grid on
    axis padded

    figure (7)
    plot(Cd_vec,Cl_vec,'k-','LineWidth',2);
    xlabel('$C_D$','Interpreter','latex');
    ylabel('$C_L$','Interpreter','latex');
    grid on 
    axis padded
end
end