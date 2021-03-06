function [ Rr_3D, Er_vec, Er_min_vec] = ConcavityRrOnX()
%MAINFUNCTION Summary of this function goes here
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

theta_bar_vec = linspace(0.3,0.3,1);
% theta_bar_vec = linspace(0.1,0.5,1);
Initial = zeros(length(theta_bar_vec),1);

Rs_vec = Initial;

Ts_vec = Initial;
Rr_vec = Initial;
Tr_vec = Initial;
Reg_up_vec = Initial;
Reg_down_vec = Initial;
Alpha_s_vec = Initial;
Alpha_r_vec = Initial;
Uber_vec = Initial;
Er_vec = Initial;
Er_min_vec = Initial;
thr_Er_nash_vec = Initial;

x_opt_vec = Initial;

rd_vec = linspace(0,1,101);
% rd_vec = linspace(0.46,0.46,1);
Rr_3D = zeros(length(rd_vec),101);

% \regcharging station revenue function
    function [r] = R_reg(tr, ts)
        if tr/p_a<ts/p_d
            r = (exp(-tr*const_r)-exp(-(ts-tr)*const_s))*(tr+Er)*c_b;
        else
            r = 0;
        end
    end

% \simplecharging station revenue funtion
    function [r] = R_simple(tr,ts)
        if tr/p_a < ts/p_d
            r = exp(-(ts-tr)*const_s)*(ts-t)*c_b;
        else
            r = exp(-ts*const)*(ts-t)*c_b;
        end
    end

% \simplecharging station revenue when regulation station is bankrupt 
    function [r] = R_simple_noreg(ts)
        r = exp(-ts*const)*(ts-t)*c_b;
    end

for idx_rd = 1:1:length(rd_vec)
    rd = rd_vec(idx_rd);
    ru = rd*1.1+0.94;
% for tht_id = 1:1:length(theta_bar_vec)
%     theta_bar = theta_bar_vec(tht_id);
    tht_id = 1;
    theta_bar = 0.3;
    Rr_temp = 0;
    
    x_vec = linspace(0,1,101);
    for i = 1:1:101
        x = x_vec(i);
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
        
        % revenue from recharging one kWh of energy into a regulating EV:
        
        Er = t*p_d*(rho_u*ru*x-rho_d*(1-rd)*(1-x)-x)/p_bar; % unit is $/kWh
        Er_min = -(t+(p_d-p_a)*theta_bar/c_b)*p_a/p_d; % minimum acceptable Er
        
        ts_temp = t+(p_d-p_a)*theta_bar/c_b;
        % Er threshold, above which the regulation station offers free recharging
        thr_Er = (1-exp(-const_s*ts_temp))/(const_r+const_s*exp(-const_s*ts_temp));
        
        Rr_matrix_tr0(tht_id,i) = R_reg(0, ts_temp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----different cases depending on Er -----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        if Er <= Er_min % reg charging canNOT achieve positive profit
            Ts_matrix(tht_id,i) = t+p_d*theta_bar/c_b;
            Tr_matrix(tht_id,i) = NaN;
            x_matrix(tht_id,i) = NaN;
            
            Rs_matrix(tht_id,i) = R_simple_noreg(Ts_matrix(tht_id,i));
            Rr_matrix(tht_id,i) = 0;
            Alpha_s_matrix(tht_id,i) = exp(-Ts_matrix(tht_id,i)*c_b/(p_d*theta_bar));
            Alpha_r_matrix(tht_id,i) = 0;
            Reg_up_matrix(tht_id,i) = 0;
            Reg_down_matrix(tht_id,i) = 0;
            
        else % reg charging CAN achieve positive profit
            if Er <= thr_Er 
                fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts_temp-tr)*const_s)...
                            -(exp(-tr*const_r)*const_r+exp(-(ts_temp-tr)*const_s)*const_s)*(tr+Er);
                tr_temp = fzero(fun_tr_nash, [0,0.5]);
            
            else % free reg-charging is offerred 
                tr_temp =0;
            end
            
            Ts_matrix(tht_id,i) = ts_temp;
            Tr_matrix(tht_id,i) = tr_temp;
            x_matrix(tht_id,i) = x;
            
            Rs_matrix(tht_id,i) = R_simple(tr_temp,ts_temp);
            Rr_matrix(tht_id,i) = R_reg(tr_temp,ts_temp);
            Alpha_s_matrix(tht_id,i) = exp(-(ts_temp-tr_temp)*c_b/((p_d-p_a)*theta_bar));
            Alpha_r_matrix(tht_id,i) = exp(-tr_temp*c_b/(p_a*theta_bar)) - Alpha_s_matrix(tht_id,i); 
            Reg_up_matrix(tht_id,i) = Alpha_r_matrix(tht_id,i)*p_d*x;
            Reg_down_matrix(tht_id,i) = Alpha_r_matrix(tht_id,i)*p_d*(1-x);
            Uber_matrix(tht_id,i) = theta_bar*(Alpha_s_matrix(tht_id,i)*p_d+Alpha_r_matrix(tht_id,i)*p_a);
        end    
       

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %
%         %-----different cases depending on Er -----%
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         Er_vec(idx_ru,i) = Er;
%         Er_min_vec(idx_ru,i) = Er_min;
%         
%         if Er <= Er_min % reg charging canNOT achieve positive profit
%             Ts_vec(tht_id) = t+p_d*theta_bar/c_b;
%             Tr_vec(tht_id) = NaN;
%             x_opt_vec(tht_id) = NaN;
%             
%             Rs_vec(tht_id) = R_simple_noreg(Ts_vec(tht_id));
%             Rr_vec(tht_id) = 0;
%             Alpha_s_vec(tht_id) = exp(-Ts_vec(tht_id)*c_b/(p_d*theta_bar));
%             Alpha_r_vec(tht_id) = 0;
%             Reg_up_vec(tht_id) = 0;
%             Reg_down_vec(tht_id) = 0;
%             Uber_vec(tht_id) = theta_bar*Alpha_s_vec(tht_id)*p_d;
%             
%         else % reg charging CAN achieve positive profit
%             if Er <= thr_Er 
%                 fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts_temp-tr)*const_s)...
%                             -(exp(-tr*const_r)*const_r+exp(-(ts_temp-tr)*const_s)*const_s)*(tr+Er);
%                 tr_temp = fzero(fun_tr_nash, [0,0.5]);
%             
%             else % free reg-charging is offerred 
%                 tr_temp =0;
%             end
%             
%             Rr_vec(tht_id) = R_reg(tr_temp,ts_temp);
%             if Rr_vec(tht_id)>Rr_temp
% 
%                 % update the temporary maxima
%                 Rr_temp = Rr_vec(tht_id);
% 
%                 % update some interesting paramaters
%                 x_opt_vec(tht_id) = x;
%                 Er_vec(tht_id) = Er;
%                 Er_min_vec(tht_id) = Er_min;
%                 thr_Er_nash_vec(tht_id) = thr_Er;
% 
%                 % update the outputs
%                 Ts_vec(tht_id) = ts_temp;
%                 Tr_vec(tht_id) = tr_temp;
% 
%                 Rs_vec(tht_id) = R_simple(tr_temp,ts_temp);
%                 Rr_vec(tht_id) = R_reg(tr_temp,ts_temp); 
%                 Alpha_s_vec(tht_id) = exp(-(ts_temp-tr_temp)*c_b/((p_d-p_a)*theta_bar));
%                 Alpha_r_vec(tht_id) = exp(-tr_temp*c_b/(p_a*theta_bar)) - Alpha_s_vec(tht_id);
%                 Reg_up_vec(tht_id) = Alpha_r_vec(tht_id)*p_d*x;
%                 Reg_down_vec(tht_id) = Alpha_r_vec(tht_id)*p_d*(1-x);
%                 Uber_vec(tht_id) = theta_bar*(Alpha_s_vec(tht_id)*p_d+Alpha_r_vec(tht_id)*p_a);
%             end
%         end
        
       % Rr_3D (idx_ru,i)=Rr_vec;
        
    end
    Rr_3D (idx_rd,:)=Rr_matrix;
    
 end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % write the data into HUB_nash.tex
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fun_name = 'MainFunction_Nash';
% file_name = 'HUB.tex';
% file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition model/Graphics/' ;
% fid = fopen([file_path file_name],'w');
% fprintf(fid,'%s\n',['% from' fun_name '; parameters: t = ' num2str(t) '; c_b = '...
%        num2str(c_b) '; rho_d = rho_u = ' num2str(rho_d) '; gamma = ' num2str(gamma)...
%        '; rd = ' num2str(rd) '; ru = ' num2str(ru)]);
% 
%    
% theta_base_vec = theta_bar_vec;
% % comment the following line if want to show the relative value over theta
% theta_base_vec = ones(length(theta_bar_vec));  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- the revenue of the two stations plus user welfare
% 
% % simple charging station revenue
% fprintf(fid,'%s\n','% simple charging station revenue');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str(Rs_vec(tht_id)/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % reg charging station revenue over that of the simple charging station 
% fprintf(fid,'%s\n','% reg charging station revenue over simple station');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str((Rr_vec(tht_id)+Rs_vec(tht_id))/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% 
% % user welfare on top of sum of the revenues 
% fprintf(fid,'%s\n','% user welfare on top of sum of the revenues');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str((Rr_vec(tht_id)+Rs_vec(tht_id)+Uber_vec(tht_id))/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% 
% % simple charging station revenue again with index descending,
% % in order to fill the reg station revenue margin with color
% fprintf(fid,'%s\n','% reverse order of simple charging station revenue');
% for fword = 1:1:length(theta_bar_vec)
%     tht_id=length(theta_bar_vec)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str(Rs_vec(tht_id)/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % sum of the two stations revenue again with index descending, 
% % in order to fill the userwelfare margin in color
% fprintf(fid,'%s\n','% reverse order of revenue sum ');
% for fword = 1:1:length(theta_bar_vec)
%     tht_id=length(theta_bar_vec)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str((Rr_vec(tht_id)+Rs_vec(tht_id))/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%%--------- the participation probability of the two stations
% 
% % probability of chosing the simple charging station
% fprintf(fid,'%s\n','% probability of simple charging station chosen');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',num2str(Alpha_s_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % probability of chosing the reg charging station over that of the simple
% fprintf(fid,'%s\n','% probability of reg charging station chosen over simple');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str(Alpha_r_vec(tht_id)+Alpha_s_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % probability simple charging station again with index descending,
% % in order to fill the probability of reg station margin in color
% fprintf(fid,'%s\n','% probability of simple charging station chosen in reverse order');
% for fword = 1:1:length(theta_bar_vec)
%     tht_id=length(theta_bar_vec)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',num2str(Alpha_s_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%%--------- the electricity price of the two stations
% 
% % electricity price of the reg charging station
% fprintf(fid,'%s\n','% electricity price of the reg charging station');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str(Tr_vec(tht_id)/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % electricity price of the simple charging station
% fprintf(fid,'%s\n','% electricity price of the simple charging station');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_vec)
%     fprintf(fid,'%s','(',num2str(theta_bar_vec(tht_id)),',',...
%             num2str(Ts_vec(tht_id)/theta_base_vec(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% fclose(fid);
 end
