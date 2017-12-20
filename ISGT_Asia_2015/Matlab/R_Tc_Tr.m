function [R_compare, Best_trial, P_check] = R_Tc_Tr( t, r_u, r_d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Prices
% t = 0.03; % wholesale price of energy unit is '$/kWh'
% r_u = 3; % renumeration on regulation up, normalized over 't'
% r_d = 0.3; % price discound when regulation down, normalized over 't'

% Time
Delta = 0.05; % is the time duration of one regulation, unit is 'hour'

% Utility
theta_bar = 0.3; %is user's valuation over recharging power

% Capacity
C_B = 50; % is the energy demand of each EV, in 'kWh'

% Probabilities;
rho_u = 0.4;% the probability of receiving an regulation "up" demand 
rho_d = 0.4;%..."down"...
rho_n = 1-rho_u-rho_d; %..."null"..., 

% Powers
P_d = 20; % namly power of regulation down, i.e. the maximum recharging power per EV, unit is 'kW'


% Revenue function
function [R] = Revenue (a_r, t_r, a_c, t_c)
    R = a_r*C_B*(t_r+E_r/(P_bar*Delta)) + a_c*C_B*(t_c - t);
end

% Probability function 
function [A_r, A_c] = Prob (t_r, p_a, t_c, p_d)
    % T_r/P_A < T_c/P_d, so user chose to 'quit', 'regulate' or 'charge'
    if t_r <= t_c*p_a/p_d
            A_r = exp(-t_r*C_B/(p_a*theta_bar)) - exp(-t_c*C_B/(p_d*theta_bar));
            A_c = exp(-t_c*C_B/(p_d*theta_bar));
    % T_c/P_d < T_r/P_A, so user chose to 'quit', 'charge' or 'regulate'
    else
            A_r = exp(-t_r*C_B/(p_a*theta_bar));
            A_c = exp(-t_c*C_B/(p_d*theta_bar)) - exp(-t_r*C_B/(p_a*theta_bar));
    end
end

Start_w_P = 0.01;
End_w_P = 1;
Num_w_P = (End_w_P-Start_w_P)/0.01+1;
R_compare = zeros(Num_w_P,2,4);

% Compute the benchmark case revenue
T_0 = t+P_d*theta_bar/C_B;
R_0 = C_B*(T_0-t)*exp(-T_0*C_B/(P_d*theta_bar));

indx = 1;
for w_P = Start_w_P:0.01:End_w_P % named as the waight of P_n over P_d, w_p <= 1 without unit
    P_n = w_P*P_d; % the null state recharging power, unit is 'kW'
%     indx = 100*w_P;
    
    % Deducted parameters 
    P_bar = rho_d*P_d + rho_n*P_n; % Average recharging power while regulating, unit is 'kW'
    P_A = P_bar - sqrt(rho_u*(P_bar)^2+rho_d*(P_d-P_bar)^2+rho_n*(P_n-P_bar)^2)/sqrt(10);% the acknowledged '85 percentage' power in'kW'
    E_r = Delta*t*P_d*(rho_u*r_u*w_P - rho_d*(1-r_d)*(1-w_P) - w_P); % average net gain after doing one regulaiton, unit is '$/slot'
    P_check(indx)=P_A;

    % Find the theoretical optimal T_c and T_r in region 'case 1', i.e. t_r <= t_c*p_a/p_d
    fun_case1 = @(x) 1-C_B*(x+E_r/(P_bar*Delta))/(P_A*theta_bar)...
                -exp(C_B*(x/P_A - (x+t+E_r/(P_bar*Delta)+P_d*theta_bar/C_B)/P_d)/theta_bar);

    Theo_case1_Tr = fzero(fun_case1,[0 1]);
    Theo_case1_Tc = Theo_case1_Tr+t+E_r/(P_bar*Delta)+P_d*theta_bar/C_B;

    Best_theo_case1(1,2) = Theo_case1_Tr;
    Best_theo_case1(1,1) = Theo_case1_Tc;
    Best_theo_case1(1,3) = theta_bar*(P_A*exp(-Theo_case1_Tr*C_B/(P_A*theta_bar))+(P_d-P_A)*exp(-Theo_case1_Tc*C_B/(P_d*theta_bar)));
    [Theo_case1_alpha_r,Theo_case1_alpha_c] = Prob (Theo_case1_Tr, P_A, Theo_case1_Tc, P_d);
    Best_theo_case1(1,4) = Revenue (Theo_case1_alpha_r, Theo_case1_Tr, Theo_case1_alpha_c, Theo_case1_Tc);
 
    R_compare(indx,1,1)=Best_theo_case1(1,1);
    R_compare(indx,1,2)=Best_theo_case1(1,2);
    R_compare(indx,1,3)=Best_theo_case1(1,3);
    R_compare(indx,1,4)=Best_theo_case1(1,4);
    % Find the theoretical optimal T_c and T_r in region 'case 2', i.e. t_r > t_c*p_a/p_d
    fun_case2 = @(y) 1-C_B*(y-t)/(P_d*theta_bar)...
                -exp(C_B*(y/P_d - (y-t-E_r/(P_bar*Delta)+P_A*theta_bar/C_B)/P_A)/theta_bar);

    Theo_case2_Tc = fzero(fun_case2,[0 1]);
    Theo_case2_Tr = Theo_case2_Tc-t-E_r/(P_bar*Delta)+P_A*theta_bar/C_B;

    Best_theo_case2(1,1) = Theo_case2_Tc;
    Best_theo_case2(1,2) = Theo_case2_Tr;
    Best_theo_case2(1,3) = theta_bar*(P_d*exp(-Theo_case2_Tc*C_B/(P_d*theta_bar))-(P_d-P_A)*exp(-Theo_case2_Tr*C_B/(P_A*theta_bar)));
    [Theo_case2_alpha_r,Theo_case2_alpha_c] = Prob (Theo_case2_Tr, P_A, Theo_case2_Tc, P_d);
    Best_theo_case2(1,4) = Revenue (Theo_case2_alpha_r, Theo_case2_Tr, Theo_case2_alpha_c, Theo_case2_Tc);

    R_compare(indx,2,:)=Best_theo_case2;
    
    indx = indx+1;
end


%         This function records the revenue on T_c-T_r plane
%         Define the search region of T_c times T_r
        Start = 0;
        Step = 0.01;
        End = 1;
        Num = (End-Start)/Step+1;
        R_matrix = zeros(Num,Num,3);
        i = 1;
        for T_c = Start:Step:End
            j = 1;

            for T_r = Start:Step:End

                [alpha_r,alpha_c] = Prob (T_r, P_A, T_c, P_d);

                R_matrix (i,j,3) = Revenue (alpha_r, T_r, alpha_c, T_c);
                R_matrix (i,j,1) = T_c;       
                R_matrix (i,j,2) = T_r;

                j = j+1;
            end
            i = i+1;

        end

%         Get the max revenue on the plane
        Best_trial = zeros(1,3);
        for m = 1:1:Num
            for n = 1:1:Num
                if R_matrix (m,n,3) > Best_trial
                    Best_trial = R_matrix (m,n,:);
                end
            end
        end
    
end

