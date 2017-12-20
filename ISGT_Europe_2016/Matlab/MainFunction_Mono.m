function [ x_opt_mono,Rs_mono,Rr_mono,Ts_mono,Tr_mono,Reg_up_mono,Reg_down_mono,...
         Alpha_s_mono,Alpha_r_mono,Er_mono,Uber_mono,R_sim] = MainFunction_Mono( ru, rd )
%MAINFUNCTION_MONO Summary of this function goes here
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

 theta_bar_mono = linspace(0.1,0.1,1);
% theta_bar_mono = linspace(0.1,0.5,41);
x_vec = linspace(0,1,101);
Initial = zeros(length(theta_bar_mono),1);
R_sim = Initial;

Rs_mono = Initial;
Ts_mono = Initial;
Rr_mono = Initial;
Tr_mono = Initial;
Reg_up_mono = Initial;
Reg_down_mono = Initial;
Alpha_s_mono = Initial;
Alpha_r_mono = Initial;
x_opt_mono = Initial;
Uber_mono = Initial;
Er_mono = Initial;


% revenue from the \simplecharging section
    function [r] = R_simple(tr,ts)
        if tr/p_a < ts/p_d
            r = exp(-(ts-tr)*const_s)*(ts-t)*c_b;
        else
            r = exp(-ts*const)*(ts-t)*c_b;
        end
    end

% revenue from the \regcharging section
    function [r] = R_reg(tr, ts)
        if tr/p_a<ts/p_d
            r = (exp(-tr*const_r)-exp(-(ts-tr)*const_s))*(tr+Er)*c_b;
        else
            r = 0;
        end
    end
% revenue when regulation section is bankrupt 
    function [r] = R_simple_noreg(ts)
        r = exp(-ts*const)*(ts-t)*c_b;
    end

% compute for each theta_bar from 0.1 to 0.5
for tht_id = 1:1:length(theta_bar_mono)
    theta_bar = theta_bar_mono(tht_id);
    R_temp = 0;
    
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
        
        % electricity price for the two charging options
        ts_temp = t+p_d*theta_bar/c_b;
        tr_temp = p_a*theta_bar/c_b-Er;
        
        Rs_temp = R_simple(tr_temp,ts_temp);
        Rr_temp = R_reg(tr_temp,ts_temp);
        
        if Rs_temp + Rr_temp > R_temp
            
            % update the temporary maxima
            R_temp = Rs_temp + Rr_temp;
            
            % update some interesting paramaters
            x_opt_mono(tht_id) = x;
            Er_mono(tht_id) = Er;
            
            % update the outputs
            Tr_mono(tht_id) = tr_temp;
            Ts_mono(tht_id) = ts_temp;
            
            Rs_mono(tht_id) = Rs_temp;
            Rr_mono(tht_id) = Rr_temp; 
            Alpha_s_mono(tht_id) = exp(-(ts_temp-tr_temp)*c_b/((p_d-p_a)*theta_bar));
            Alpha_r_mono(tht_id) = exp(-tr_temp*c_b/(p_a*theta_bar)) - Alpha_s_mono(tht_id);
            Reg_up_mono(tht_id) = Alpha_r_mono(tht_id)*p_d*x;
            Reg_down_mono(tht_id) = Alpha_r_mono(tht_id)*p_d*(1-x);
            Uber_mono(tht_id) = theta_bar*(Alpha_s_mono(tht_id)*p_d+Alpha_r_mono(tht_id)*p_a);
        end
    end
    
% with 'x' optimized, compare the both-choice option with uni-choice option
    R_sim(tht_id) = R_simple_noreg(Ts_mono(tht_id));
    if Rs_mono(tht_id)+Rr_mono(tht_id) <= R_sim(tht_id)
        
        Tr_mono(tht_id) = NaN;
        x_opt_mono(tht_id) = NaN;
        
%         Rs_mono(tht_id) = R_simple_noreg(Ts_mono(tht_id));
        Rr_mono(tht_id) = NaN;
%         Alpha_s_mono(tht_id) = exp(-(Ts_mono(tht_id))*c_b/(p_d*theta_bar));
        Alpha_r_mono(tht_id) = NaN;
        Reg_up_mono(tht_id) = NaN;
        Reg_down_mono(tht_id) = NaN;
        Uber_mono(tht_id) = theta_bar*Alpha_s_mono(tht_id)*p_d;
    end
        

end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % write the data into HUB_mono.tex
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fun_name = 'MainFunction_Mono';
% file_name = 'HUB_mono.tex';
% file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition model/Graphics/' ;
% fid = fopen([file_path file_name],'w');
% fprintf(fid,'%s\n',['% from' fun_name '; parameters: t = ' num2str(t) '; c_b = '...
%        num2str(c_b) '; rho_d = rho_u = ' num2str(rho_d) '; gamma = ' num2str(gamma)...
%        '; rd = ' num2str(rd) '; ru = ' num2str(ru)]);
% 
% theta_base_mono = theta_bar_mono;
% % comment the following line if want to show the relative value over theta
% theta_base_mono = ones(length(theta_bar_mono));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- the revenue of the two sections plus user welfare
% 
% % simple charging section revenue
% fprintf(fid,'%s\n','% simple charging section revenue');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str(Rs_mono(tht_id)/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % reg charging section revenue over that of the simple charging section 
% fprintf(fid,'%s\n','% reg charging section revenue over simple section');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str((Rr_mono(tht_id)+Rs_mono(tht_id))/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % user welfare on top of sum of the revenues 
% fprintf(fid,'%s\n','% user welfare on top of sum of the revenues');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str((Rr_mono(tht_id)+Rs_mono(tht_id)+Uber_mono(tht_id))/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % simple charging section revenue again with index descending,
% % in order to fill the reg section revenue margin with color
% fprintf(fid,'%s\n','% reverse order of simple charging section revenue');
% for fword = 1:1:length(theta_bar_mono)
%     tht_id=length(theta_bar_mono)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str(Rs_mono(tht_id)/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % sum of the two sections revenue again with index descending, 
% % in order to fill the userwelfare margin in color
% fprintf(fid,'%s\n','% reverse order of revenue sum ');
% for fword = 1:1:length(theta_bar_mono)
%     tht_id=length(theta_bar_mono)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str((Rr_mono(tht_id)+Rs_mono(tht_id))/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%%--------- the participation probability of the two sections
% 
% % probability of chosing the simple charging section
% fprintf(fid,'%s\n','% probability of simple charging section chosen');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',num2str(Alpha_s_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % probability of chosing the reg charging section over that of the simple
% fprintf(fid,'%s\n','% probability of reg charging section chosen over simple');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str(Alpha_r_mono(tht_id)+Alpha_s_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % probability simple charging section again with index descending,
% % in order to fill the probability of reg section margin in color
% fprintf(fid,'%s\n','% probability of simple charging section chosen in reverse order');
% for fword = 1:1:length(theta_bar_mono)
%     tht_id=length(theta_bar_mono)+1-fword;
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',num2str(Alpha_s_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%%--------- the electricity price of the two sections
% 
% % electricity price of the reg charging section
% fprintf(fid,'%s\n','% electricity price of the reg charging section');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str(Tr_mono(tht_id)/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% % electricity price of the simple charging section
% fprintf(fid,'%s\n','% electricity price of the simple charging section');
% fprintf(fid,'%s\n','\addplot coordinates{');
% for tht_id = 1:1:length(theta_bar_mono)
%     fprintf(fid,'%s','(',num2str(theta_bar_mono(tht_id)),',',...
%             num2str(Ts_mono(tht_id)/theta_base_mono(tht_id)),')');
% end
% fprintf(fid,'%s\n', '};');
% 
% fclose(fid);
end

