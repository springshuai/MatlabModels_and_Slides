function [ x_opt, tr_matrix, Rr_matrix, ts_nash, ts_idx, userfare] = BestRespondingTr(ru, rd ,ts_input)
%BESTRESPONDINGTR Summary of this function goes here
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

global  p_d p_n p_bar p_a const const_r const_s

% \regcharging station revenue function
    function [r] = R_reg(tr, ts)
        if tr/p_a<ts/p_d
            r = (exp(-tr*const_r)-exp(-(ts-tr)*const_s))*(tr+Er)*c_b;
        else
            r = 0;
        end
    end


%
ts_idx = linspace(ts_input,ts_input,1); % simplecharging price higher than wholesale price 't'
x = linspace(0,1,101);
tr_vec = zeros(length(ts_idx),1);
Rr_vec = zeros(length(ts_idx),1);

x_opt = zeros(length(ts_idx),1);
tr_opt = zeros(length(ts_idx),1);
tr_matrix = zeros(length(ts_idx),length(x));
Rr_opt = zeros(length(ts_idx),1);
Rr_matrix = zeros(length(ts_idx),length(x));
Rr_temp = 0;
userfare = zeros(length(ts_idx),1);
for i = 1:1:length(ts_idx)
    for j=1:1:length(x)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----powers and depending parameters -----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        p_d = 20; % maximum power
        p_n = x(j)*p_d; % regulation-null power
        p_bar = rho_d*p_d + rho_n*p_n; % average power obtained 
        p_a = p_bar - gamma*sqrt(rho_u*(p_bar)^2+rho_d*(p_d-p_bar)^2+rho_n*(p_n-p_bar)^2);
        const = c_b/(p_d*theta_bar);
        const_r = c_b/(p_a*theta_bar);
        const_s = c_b/((p_d-p_a)*theta_bar);
        Er = t*p_d*(rho_u*ru*x(j)-rho_d*(1-rd)*(1-x(j))-x(j))/p_bar; % unit is $/kWh
        ts_nash = t+(p_d-p_a)*theta_bar/c_b;
        
        fun_tr = @(tr) exp(-tr*const_r)-exp(-(ts_idx(i)-tr)*const_s)...
        -(exp(-tr*const_r)*const_r+exp(-(ts_idx(i)-tr)*const_s)*const_s)*(tr+Er);
    
% 'thr_Er'is the threshold above which FREE regcharging is the optimal choice;
% this threshold increases as 'ts' increase.
%      thr_Er(i) = (1-exp(-const_s*ts_idx(i)))/(const_r+const_s*exp(-const_s*ts_idx(i)));
      
% the solutions in 'tr_vec' should give the theoretical best-responding tr
      tr_vec(i) = fzero(fun_tr,[0,1]);
      if  tr_vec(i)/p_a > ts_idx(i)/p_d
          tr_vec(i) = NaN;
      end
% 'Rr_vec' gives the revenue at the theoretical best-responding 'tr_vec'
      Rr_vec(i) = R_reg(tr_vec(i), ts_idx(i));
      if Rr_vec(i) > Rr_temp
          x_opt(i) = x(j);
          tr_opt(i) = tr_vec(i);
          Rr_opt(i) = Rr_vec(i);
          Rr_temp = Rr_vec(i);
      end
      
      tr_matrix(i,j) = tr_vec(i);
      Rr_matrix(i,j) = Rr_vec(i);
    end
    
end
end

