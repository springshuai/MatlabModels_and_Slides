function [ Er, social_theory, social_nash, social_mono, Reg_nash, Reg_mono] = Social( x, ru, rd )
%SOCIAL_OPT Summary of this function goes here
%   Detailed explanation goes here

global t c_b rho_d rho_u rho_n gamma theta_bar

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----powers and depending parameters -----%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  p_d p_n p_bar p_a const const_r const_s

p_d = 20; % maximum power
p_n = x*p_d; % regulation-null power
p_bar = rho_d*p_d + rho_n*p_n; % average power obtained 
p_a = p_bar - gamma*sqrt(rho_u*(p_bar)^2+rho_d*(p_d-p_bar)^2+rho_n*(p_n-p_bar)^2);
const = c_b/(p_d*theta_bar);
const_r = c_b/(p_a*theta_bar);
const_s = c_b/((p_d-p_a)*theta_bar);
Er = t*p_d*(rho_u*ru*x-rho_d*(1-rd)*(1-x)-x)/p_bar; % unit is $/kWh
ts_nash = t+(p_d-p_a)*theta_bar/c_b;
thr_Er_nash = (1-exp(-const_s*ts_nash))/(const_r+const_s*exp(-const_s*ts_nash));
Rr_pos = ts_nash*p_a/p_d+Er;

% \regcharging station revenue function
    function [r] = R_reg(ts, tr)
        if tr/p_a<ts/p_d
            r = (exp(-tr*const_r)-exp(-(ts-tr)*const_s))*(tr+Er)*c_b;
        else
            r = 0;
        end
    end

% \simplecharging station revenue funtion
function [r] = R_simple(ts,tr)
    if tr/p_a < ts/p_d
        r = exp(-(ts-tr)*const_s)*(ts-t)*c_b;
    else
        r = exp(-ts*const)*(ts-t)*c_b;
    end
end

% ts_idx = linspace(0,0.15,16);
% tr_idx = linspace(0,0.05,51);
% social_matrix = zeros(length(ts_idx), length(tr_idx));
% for i = 1:1:16
% 
%     for j = 1:1:51
%         social_matrix(i,j) = R_reg(ts_idx(i), tr_idx(j)) + R_simple(ts_idx(i), tr_idx(j))...
%                 + User_welfare( theta_bar, c_b, p_d, p_a, ts_idx(i), tr_idx(j));
%     end
% end

social_theory = R_reg(t, -Er) + R_simple(t, -Er)...
            + User_welfare( theta_bar, c_b, p_d, p_a, t, -Er);
        
[t_nash,R_nash,A_nash,Reg_nash,uber_nash,t_mono,R_mono,A_mono,Reg_mono,uber_mono] = optTrTs(x, rd, ru);
social_nash = R_nash(1)+R_nash(2)+uber_nash;
social_mono = R_mono+uber_mono;
end

