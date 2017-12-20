function [ R_compare, U_bar, w_P_vec, Ar_vec, Ac_vec, R_matrix, P_A] = R_Tc_Tr_Correction( gamma, r_u, r_d, rho_u, rho_d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t=0.03;
r_u_min=2-rho_u+gamma*rho_u^(-0.5)*(1-rho_u)^1.5;
r_d_min=1-(rho_d-gamma*sqrt(rho_d-(rho_d^2)));
% Time
Delta = 0.05; % is the time duration of one regulation, unit is 'hour'

% Utility
theta_bar = 0.3; %is user's valuation over recharging power

% Capacity
C_B = 50; % is the energy demand of each EV, in 'kWh'

% Probabilities;
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
            A_c = exp(-(t_c-t_r)*C_B/((p_d-p_a)*theta_bar));
            A_r = exp(-t_r*C_B/(p_a*theta_bar)) - A_c;
    else 
            A_c = exp(-t_c*C_B/(p_d*theta_bar));
            A_r = 0;
%             A_r = exp(-t_r*C_B/(p_a*theta_bar));
%             A_c = exp(-t_c*C_B/(p_d*theta_bar)) - exp(-t_r*C_B/(p_a*theta_bar));
    end
end

w_P_Start = 0;
w_P_End = 1;
w_P_Step = 0.05;
w_P_Num = length(w_P_Start:w_P_Step:w_P_End);
R_compare = zeros(3,3,w_P_Num);
slope=rho_u*r_u+rho_d*(1-r_d)-1;
% r_u_min=2-rho_u+gamma*rho_u^(-0.5)*(1-rho_u)^1.5;
% r_d_min=rho_d-gamma*sqrt(rho_d-(rho_d^2));
Ar_vec = zeros;
Ac_vec = zeros;
for id = 1:1:w_P_Num

    w_P = w_P_Start+(id-1)*w_P_Step; % named as the waight of P_n over P_d, w_p <= 1 without unit
    w_P_vec(id) = w_P;    
    P_n = w_P*P_d; % the null state recharging power, unit is 'kW'

    E_r = Delta*t*P_d*(rho_u*r_u*w_P - (1-r_d)*(1-w_P)*rho_d - w_P); % average net gain after doing one regulaiton, unit is '$/slot'    

    P_bar = rho_d*P_d + rho_n*P_n; % Average recharging power while regulating, unit is 'kW'
    P_A = P_bar - gamma*sqrt(rho_u*(P_bar)^2+rho_d*(P_d-P_bar)^2+rho_n*(P_n-P_bar)^2);% the acknowledged '85 percentage' power in'kW'
    
    T_c_theo = t + P_d*theta_bar/C_B;
    T_r_theo = P_A*theta_bar/C_B-E_r/(P_bar*Delta);
    fun1(id) = T_c_theo/P_d - T_r_theo/P_A;
    fun2(id) = rho_u*r_u*w_P-rho_d*(1-r_d)*(1-w_P)-w_P+(rho_d+rho_n*w_P)*(rho_d+rho_n*w_P-gamma*sqrt(rho_n*w_P^2+rho_d-(rho_d+rho_n*w_P)^2));  
    if fun1(id)*fun2(id)<0
        Error(id)=1;
    end
        
    [A_r_theo, A_c_theo] = Prob (T_r_theo, P_A, T_c_theo, P_d);
    Ar_vec(id) = A_r_theo;
    Ac_vec(id) = A_c_theo;
    R_compare(1,1,id) = T_c_theo;
    if fun1(id)>0
        R_compare(1,2,id) = T_r_theo;
    else
        R_compare(1,2,id) = inf;
    end
    R_compare(1,3,id) = Revenue (A_r_theo, T_r_theo, A_c_theo, T_c_theo); 
 
    % Compute the benchmark case revenue
    T_0 = t+P_d*theta_bar/C_B;
    R_0 = C_B*(T_0-t)*exp(-T_0*C_B/(P_d*theta_bar));
    R_compare(2,1,id)=T_0;
    R_compare(2,2,id)=NaN;
    R_compare(2,3,id)=R_0;
    
% The following part records the revenue on T_c-T_r plane and find the max
% Define the search region of T_c times T_r on a 0.5*0.5 plane
    Start = 0;
    End = 0.5;
    Step = 0.01;
    Num = length(Start:Step:End);
    for c = 1:1:Num
        tc = Start+(c-1)*Step;
        for r = 1:1:Num
            tr = Start+(r-1)*Step;
            [alpha_r,alpha_c] = Prob (tr, P_A, tc, P_d);
            R_exh (c,r,1) = tc;       
            R_exh (c,r,2) = tr;
            R_exh (c,r,3) = Revenue (alpha_r, tr, alpha_c, tc);            
        end
    end
    
% Define the search region of T_c times T_r around the theoretical optimal
% 'T_c_theo' and 'T_r_theo'    
    tc_Start = 0;
    tc_End = 0.2;
    tr_Start = 0;
    tr_End = 0.2;
    Step = 0.01;

    tc_Num = length(tc_Start:Step:tc_End);
    tr_Num = length(tr_Start:Step:tr_End);

    for i = 1:1:tc_Num
        T_c = tc_Start+(i-1)*Step;

        for j = 1:1:tr_Num
            T_r = tr_Start+(j-1)*Step;

            [alpha_r,alpha_c] = Prob (T_r, P_A, T_c, P_d);

            R_matrix (i,j,1) = T_c;       
            R_matrix (i,j,2) = T_r;
            R_matrix (i,j,3) = Revenue (alpha_r, T_r, alpha_c, T_c);
            
            % Get the max revenue on the plane
            if R_matrix (i,j,3) > R_compare(3,:,id)
                R_compare(3,1,id) = R_matrix (i,j,1);
                R_compare(3,2,id) = R_matrix (i,j,2);
                R_compare(3,3,id) = R_matrix (i,j,3);
            end            

        end

    end
    U_bar(1,id)=User_welfare( theta_bar, C_B, P_d, P_A, T_c_theo, T_r_theo );
    U_bar(2,id)=User_welfare( theta_bar, C_B, P_d, P_A, T_c_theo, inf );
    U_bar(3,id)=User_welfare( theta_bar, C_B, P_d, P_A, R_compare(3,1,id), R_compare(3,2,id) );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%------write the result into a Tex file named 3D_Tc_Tr.tex----%%%%%%
% Parameters = ['$P_n/P_d = $' num2str(w_P) ', $\gamma = $' num2str(gamma) ', $t = $' num2str(t) ...
%     ', $\rho_d = $' num2str(rho_d) ', $\rho_u =$' num2str(rho_u) ', $r_d = $' ...
%     num2str(r_d) '%r_u = ' num2str(r_u)];
% x_label = '$T_c$';
% y_label = '$T_r$';
% z_label = 'Revenue';
% 
% name_file = '3D_Tc_Tr.tex';
% fid = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' name_file], 'w');
% % fprintf(fid, '%s\n', ['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
% %     x_label ',ylabel=' y_label ']']);
% 
% fprintf(fid, '%s\n', Parameters);
% 
% fprintf(fid, '%s\n', '%this is the revenue on the big 0.5*0.5 plane, step is 0.01, maximized (1.7969) at tc = 0.15, tr = 0.02')
% for c=1:1:Num
%     for r=1:1:Num
%         fprintf(fid, '%s','(',num2str(R_exh(c,r,1)),',',num2str(R_exh(c,r,2)),',',num2str(R_exh(c,r,3)),')');
%     end
% end
% 
% fprintf(fid, '%s\n', '%this is the revenue on the small 0.02*0.02 plane, step is 0.001, maximized (1.7999) at tc = 0.15, tr = 0.0225')
% for i = 1:1:tc_Num
%     for j = 1:1:tr_Num
%         fprintf(fid, '%s','(',num2str(R_matrix(i,j,1)),',',num2str(R_matrix(i,j,2)),',',num2str(R_matrix(i,j,3)),')');
%     end
% end
% fclose(fid);


end

