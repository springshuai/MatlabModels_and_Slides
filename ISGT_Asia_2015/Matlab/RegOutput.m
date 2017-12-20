function [Output, Sign] = RegOutput(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

P_d = 20;
% Probabilities;
rho_u = 0.5;% the probability of receiving an regulation "up" demand 
rho_d = 0.5;%..."down"...
rho_n = 1-rho_u-rho_d; %..."null"..., 

gamma = 0;

t=0.03;
r_u_min=2-rho_u+gamma*rho_u^(-0.5)*(1-rho_u)^1.5;
r_d_min=1-(rho_d-gamma*sqrt(rho_d-(rho_d^2)));
% rd_Step = 0.1*r_d_min;
% ru_Step = 0.1*r_u_min;

rd_Step = 0.025;
ru_Step = 0.1;

rd_Start = 0.4;
rd_End =0.8;
rd_Num = length(rd_Start:rd_Step:rd_End);

ru_Start = 1.5;
ru_End = 2.5;
ru_Num = length(ru_Start:ru_Step:ru_End);

Output = zeros;
    for d=1:1:rd_Num
    r_d(d)=rd_Start+(d-1)*rd_Step;
    
        for u=1:1:ru_Num
        r_u(u)=ru_Start+(u-1)*ru_Step;
        [ R_compare, slope, fun1, fun2, U_bar, w_P_vec, Ar_vec, Ac_vec] = R_Tc_Tr_Correction( t, gamma, r_u(u), r_d(d), rho_u, rho_d);
        Num_x = length (w_P_vec);
            if R_compare(1,3,1)>R_compare(1,3,Num_x)
                Sign(d,u) = -1;
            elseif R_compare(1,3,1)<R_compare(1,3,Num_x)
                Sign(d,u) = 1;
            elseif R_compare(1,3,1) == 0
                Sign(d,u) = 0;
            else
                Sign(d,u) = 2;
            end
            
        Output(d,u) = max(P_d*Ar_vec(1),P_d*Ar_vec(Num_x));
        end
        
    end

end

