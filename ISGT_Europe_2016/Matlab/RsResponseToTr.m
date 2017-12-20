function [ Rs_matrix,Ts_theo,Rs_theo] = RsResponseToTr( tr,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global t c_b rho_d rho_u rho_n gamma p_d p_n p_bar p_a const const_r const_s

% whole sale prices $/kWh
t = 0.03;

% energy demand per car kWh
c_b = 50;

% probabilities:
rho_d=0.48;
rho_u=0.48;
rho_n = 1-rho_d-rho_u;

gamma = 0.05; % user sensitiveity and preference
p_d = 20; % maximum power

theta_bar = 0.3;
 
Ts_vec = linspace(0.03,0.16,131);
Rs_matrix = zeros(length(Ts_vec),1);


% \simplecharging station revenue funtion
    function [r] = R_simple(tr,ts)
        if tr/p_a < ts/p_d
            r = exp(-(ts-tr)*const_s)*(ts-t)*c_b;
        else
            r = exp(-ts*const)*(ts-t)*c_b;
        end
    end

for ts_id = 1:1:length(Ts_vec)
        ts = Ts_vec(ts_id);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----powers and depending parameters -----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_n = x*p_d; % regulation-null power
        p_bar = rho_d*p_d + rho_n*p_n; % average power obtained 
        
        % user perceived recharging power:
        p_a = p_bar - gamma*sqrt(rho_u*(p_bar)^2+rho_d*(p_d-p_bar)^2+rho_n*(p_n-p_bar)^2);
        
        const = c_b/(p_d*theta_bar);
        const_r = c_b/(p_a*theta_bar);
        const_s = c_b/((p_d-p_a)*theta_bar);
        
        Rs_matrix(ts_id)=R_simple(tr,ts);
end
    
if tr < (t+(p_d-p_a)*theta_bar/c_b)*p_a/p_d
    Ts_theo = t+(p_d-p_a)*theta_bar/c_b;
elseif tr > (t+p_d*theta_bar/c_d)*p_a/p_d
    Ts_thero = t+p_d*theta_bar/c_d;
else
    Ts_theo = tr*p_d/p_a;
end
Rs_theo = R_simple(tr,Ts_theo);
            
end
        

