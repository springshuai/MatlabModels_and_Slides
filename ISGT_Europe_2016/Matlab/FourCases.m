function [Tr_vec,Ts_vec,Er,Er_th1,Er_th2] = FourCases( rd, ru, theta_bar)
%UNTITLED Summary of this function goes here
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

ts_vec = linspace(0.04,0.12,81);
tr_vec = linspace(-0.05,0.05,101);

    function [tr]=Tr_Br (ts,Er)
            Er_th1 = (1-exp(-const_s*ts))/(const_r+const_s*exp(-const_s*ts));
            Er_th2 = Er_th1*(1+(p_d/p_a-1)*exp(ts*const_s));
%             Er_th2 = 1/const_s*(exp(const_s*ts)-1);
            
        if ts < -Er*p_d/p_a
            tr = ts*p_a/p_d;
        elseif Er <= Er_th1
            fun_tr_nash = @(tr) exp(-tr*const_r)*(1-const_r*(tr+Er))...
                -exp(-(ts-tr)*const_s)*(1+(tr+Er)*const_s);
%             fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts-tr)*const_s)...
%             -(exp(-tr*const_r)*const_r+exp(-(ts-tr)*const_s)*const_s)*(tr+Er);
            tr = fzero(fun_tr_nash, [0,0.5]);
        elseif Er <= Er_th2
            tr = 0;
        else
            fun_tr_prime = @(tr)1-exp(-(ts-tr)*const_s)*(1+(tr+Er)*const_s);
%             fun_tr_prime = @(tr)1-(exp(-tr*const_r)*const_r+exp(-(ts-tr)*const_s)*const_s)*(tr+Er);
            tr=fzero(fun_tr_prime,[-1,0]);
        end
    end
        
    function [ts]=Ts_Br (tr)
        if tr < (t+(p_d-p_a)*theta_bar/c_b)*p_a/p_d
            ts = t+(p_d-p_a)*theta_bar/c_b;
        elseif tr <= (t+p_d*theta_bar/c_b)*p_a/p_d
            ts = tr*p_d/p_a;
        else
            ts = t+p_d*theta_bar/c_b;
        end
    end
           

    x = 0.5;
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
        
        % revenue from recharging one kWh of energy into a regulating EV
        
        Er = t*p_d*(rho_u*ru*x-rho_d*(1-rd)*(1-x)-x)/p_bar; % unit is $/kWh
        Er_12 = -(t+p_d*theta_bar/c_b)*p_a/p_d; %
        Er_23 = -(t+(p_d-p_a)*theta_bar/c_b)*p_a/p_d; %
%         Er_34 = (1-exp(-const_s*ts_temp))/(const_r+const_s*exp(-const_s*ts_temp));
        
%         ts_1 = t+p_d*theta_bar/c_b;
%         ts_2 = -Er*p_d/p_a;
%         ts_3 = t+(p_d-p_a)*theta_bar/c_b;
%         ts_4 = t+(p_d-p_a)*theta_bar/c_b;
%         
%         tr_1 = (p_a/p_d)*ts_1+0.01;
%         tr_2 = -Er;
%         fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts_3-tr)*const_s)...
%         -(exp(-tr*const_r)*const_r+exp(-(ts_3-tr)*const_s)*const_s)*(tr+Er);
%         tr_3 = fzero(fun_tr_nash, [0,0.5]);
%         tr_4 = 0;

for i = 1:1:length(ts_vec)
    ts = ts_vec(i);
    Tr_vec(i) = Tr_Br (ts,Er);
end

for j = 1:1:length(tr_vec)
    tr = tr_vec(j);
    Ts_vec(j) = Ts_Br (tr);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % write the data into HUB_nash.tex
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun_name = 'FourCases';
file_name = 'FourCases.tex';
file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition_model/Graphics/' ;
fid = fopen([file_path file_name],'w');
fprintf(fid,'%s\n',['% from' fun_name '; parameters: theta = ' num2str(theta_bar) '; rd = '...
       num2str(rd) '; ru = ' num2str(ru)]);
   
fprintf(fid,'%s\n','% best-response Tr');
fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:length(ts_vec)   
    fprintf(fid,'%s','(',num2str(ts_vec(i)),',',num2str(Tr_vec(i)),')');
end
fprintf(fid,'%s\n', '};');  

fprintf(fid,'%s\n','% best-response Ts');
fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:length(tr_vec)   
    fprintf(fid,'%s','(',num2str(Ts_vec(i)),',',num2str(tr_vec(i)),')');
end
fprintf(fid,'%s\n', '};');

end


