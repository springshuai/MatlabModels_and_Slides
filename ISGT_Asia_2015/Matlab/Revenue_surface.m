function [R_matrix, P_A] = Revenue_surface()
% Record the revenue on T_c-T_r plane
% Define the search region of T_c times T_r
 Start = 0;
 Step = 0.01;
 End = 1;
 Num = (End-Start)/Step+1;
 R_matrix = zeros(Num,Num,3);
gamma=0.1;
t=0.03;

Delta = 0.05; % is the time duration of one regulation, unit is 'hour'

theta_bar = 0.3; %is user's valuation over recharging power

C_B = 50; % is the energy demand of each EV, in 'kWh'

rho_d=0.5;
rho_u=0.5;
rho_n = 1-rho_u-rho_d; %..."null"..., 
r_u=2;
r_d=0.5;
P_d = 20; % namly power of regulation down, i.e. the maximum recharging power per EV, unit is 'kW'
w_P=0.8;
P_n = w_P*P_d;
E_r = Delta*t*P_d*(rho_u*r_u*w_P - (1-r_d)*(1-w_P)*rho_d - w_P); % average net gain after doing one regulaiton, unit is '$/slot'    
P_bar = rho_d*P_d + rho_n*P_n; % Average recharging power while regulating, unit is 'kW'
P_A = P_bar - gamma*sqrt(rho_u*(P_bar)^2+rho_d*(P_d-P_bar)^2+rho_n*(P_n-P_bar)^2);% the acknowledged '85 perce

function [R] = Revenue (a_r, t_r, a_c, t_c)
    R = a_r*C_B*(t_r+E_r/(P_bar*Delta)) + a_c*C_B*(t_c - t);
end

% Probability function 
function [A_r, A_c] = Prob (t_r, p_a, t_c, p_d)
    % T_r/P_A < T_c/P_d, so user chose to 'quit', 'regulate' or 'charge'
    if t_r <= t_c*p_a/p_d
            A_c = exp(-(t_c-t_r)*C_B/((p_d-p_a)*theta_bar));
            A_r = exp(-t_r*C_B/(p_a*theta_bar)) - A_c;
    else 
            A_c = exp(-t_c*C_B/(p_d*theta_bar));
            A_r = 0;
%             A_r = exp(-t_r*C_B/(p_a*theta_bar));
%             A_c = exp(-t_c*C_B/(p_d*theta_bar)) - exp(-t_r*C_B/(p_a*theta_bar));
    end
end

for j=1:1:Num
    T_c = Start+(j-1)*Step;


    for i = 1:1:Num
        T_r = Start+(i-1)*Step;

        [alpha_r,alpha_c] = Prob (T_r, P_A, T_c, P_d);
        
        R_matrix (i,j,3) = Revenue (alpha_r, T_r, alpha_c, T_c);
        R_matrix (i,j,1) = T_c;       
        R_matrix (i,j,2) = T_r;

    end
end

% Get the max revenue on the plane
Best_trial = zeros(1,3);
for m = 1:1:Num
    for n = 1:1:Num
        if R_matrix (m,n,3) > Best_trial
            Best_trial = R_matrix (m,n,:);
        end
    end
end
end