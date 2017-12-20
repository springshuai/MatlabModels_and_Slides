function [T_nash,R_nash,A_nash,Reg_nash,Ubar_nash,T_coop,R_coop,A_coop,Reg_coop,Ubar_coop] = NashEquRev(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
   
global x rd ru t c_b rho_d rho_u rho_n gamma theta_bar
% regulation renumeration
rd = 0.7;
ru = 2.1;

% whole sale prices $/kWh
t = 0.03;

% energy demand per car kWh
c_b = 50;

% probabilities:
rho_d=0.48;
rho_u=0.48;
rho_n = 1-rho_d-rho_u;

% user sensitiveity and preferences
gamma = 0.05;
theta_bar = 0.3;

for i=1:1:11
    x(i) = 0+0.1*(i-1);
    [T_nash(i,:),R_nash(i,:),A_nash(i,:),Reg_nash(i,:),Ubar_nash(i,:),...
     T_coop(i,:),R_coop(i,:),A_coop(i,:),Reg_coop(i,:),Ubar_coop(i,:)]...
    = optTrTs(x(i),rd,ru);
end

% data into tex file


end

