function [Rr_matrix, Rr_max, Tr_opt, X_opt, Tr_theo, Rr_theo, X_theo] = RrResponseToTs( ru, rd )
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

theta_bar = 0.3;

x_vec = linspace(0,1,101);
Ts_vec = linspace(0.0954,0.0954,1);
Tr_vec = linspace(0,0.1,101);
Rr_matrix = zeros(length(Ts_vec),length(Tr_vec));
Rr_theo = zeros(length(Ts_vec),1);
Tr_opt = zeros(length(Ts_vec),1);
Tr_theo = zeros(length(Ts_vec),1);
X_opt = zeros(length(Ts_vec),1);
X_theo = NaN(length(Ts_vec),1);
x_ts_tr = NaN(length(Ts_vec),length(Tr_vec));
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


for ts_id = 1:1:length(Ts_vec)
    ts = Ts_vec(ts_id);
    
    for tr_id = 1:1:length(Tr_vec)
        tr = Tr_vec(tr_id);
        
        Rr_temp = -10; % a sufficiently small initial revenue, to compare 
        Rr_theo_temp = -10; % with current revenue, in order to find the 'x'
                       % that maximize the revenue for a particular tr,ts
        
        for x_id = 1:1:length(x_vec)
            x = x_vec(x_id);
            
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----compute the theoretical Tr^br and corresponding Rr-----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        Er_min = -(t+p_d*theta_bar/c_b)*p_a/p_d; % minimum acceptable Er
        thr_Er = (1-exp(-const_s*ts))/(const_r+const_s*exp(-const_s*ts));
        
        if Er <= Er_min % reg charging canNOT achieve positive profit
            Tr_theo(ts_id)=NaN;
            Rr_theo(ts_id)=NaN;
            
        else % reg charging CAN achieve positive profit
            if Er <= thr_Er 
                fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts-tr)*const_s)...
                            -(exp(-tr*const_r)*const_r+exp(-(ts-tr)*const_s)*const_s)*(tr+Er);
                Tr_theo(ts_id) = fzero(fun_tr_nash, [0,0.1]);
            
            else % free reg-charging is offerred 
                Tr_theo(ts_id) = 0;
            end
           Rr_theo(ts_id) = R_reg(Tr_theo(ts_id), ts);
        end
        
        
        Rr_matrix(ts_id,tr_id) = R_reg(tr, ts);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----find the x that maximize the Rr_matrix -----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            if Rr_matrix(ts_id,tr_id) > Rr_temp;
                Rr_temp = Rr_matrix(ts_id,tr_id);
                x_ts_tr(ts_id, tr_id) = x;
            else
                Rr_matrix(ts_id,tr_id) = Rr_temp;
            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %-----find the x that maximize the Rr_theo -----%
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            if Rr_theo(ts_id) > Rr_theo_temp;
                Rr_theo_temp = Rr_theo(ts_id);
                X_theo(ts_id) = x;
            else
                Rr_matrix(ts_id,tr_id) = Rr_temp;
            end
        end
    end
end

[Rr_max, tr_opt_id]=max(Rr_matrix,[],2);
for j = 1:1:length(Tr_opt)
    Tr_opt(j)= Tr_vec(tr_opt_id(j));
    X_opt(j)=x_ts_tr(j,tr_opt_id(j));
end


end
