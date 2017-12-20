function [R] = Revenue_toSolve (x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
syms  theta_bar guamma rho_d rho_u r_u r_d p_d t t_r  t_c C_B Delta 

    rho_n = 1-rho_u-rho_d;
    P_n = x*p_d;
    
    E_r = Delta*t*p_d*(rho_u*r_u*x - (1-r_d)*(1-x)*rho_d - x);
    
    P_bar = rho_d*p_d + rho_n*P_n;
    
    P_A = P_bar - guamma*sqrt(rho_u*(P_bar)^2+rho_d*(p_d-P_bar)^2+rho_n*(P_n-P_bar)^2);
    
    A_c = exp(-(t_c-t_r)*C_B/((p_d-P_A)*theta_bar));
    A_r = exp(-t_r*C_B/(P_A*theta_bar)) - A_c;
    
    
    R = A_r*C_B*(t_r+E_r/(P_bar*Delta)) + A_c*C_B*(t_c - t);



end

