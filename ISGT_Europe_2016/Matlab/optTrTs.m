function [t_nash,R_nash,A_nash,Reg_nash,uber_nash,t_mono,R_mono,A_mono,Reg_mono,uber_mono] = optTrTs(theta_bar, rd, ru)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global t c_b rho_d rho_u rho_n gamma 

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----powers and depending parameters -----%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  p_d p_n p_bar p_a const const_r const_s

x=1;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%---exhausging the best responding tr---%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the best reponding tr, as ts goes from 0 to 0.15;

% ts_idx = linspace(t,0.15,13); % simplecharging price higher than wholesale price 't'
% RrMax_vec = zeros(length(ts_idx),2);
% for i = 1:1:length(ts_idx)
%     fun_tr = @(tr) exp(-tr*const_r)-exp(-(ts_idx(i)-tr)*const_s)...
%         -(exp(-tr*const_r)*const_r+exp(-(ts_idx(i)-tr)*const_s)*const_s)*(tr+Er);
%     
% % 'thr_Er'is the threshold above which FREE regcharging is the optimal choice;
% % this threshold increases as 'ts' increase.
%       thr_Er(i) = (1-exp(-const_s*ts_idx(i)))/(const_r+const_s*exp(-const_s*ts_idx(i)));
%       
% % the solutions in 'tr_vec' should give the theoretical best-responding tr
%       tr_vec(i) = fzero(fun_tr,[0,1]);
%       if  tr_vec(i)/p_a > ts_idx(i)/p_d
%           tr_vec(i) = NaN;
%       end
% % 'Rr_vec' gives the revenue at the theoretical best-responding 'tr_vec'
%       Rr_vec(i) = R_reg(tr_vec(i), ts_idx(i));
% 
% % exhausive search for the best resonding tr and corresponding revenue
% % 'RrMax_vec(:,2)' gives the exhausive searched best-responding 'tr'
% % 'RrMax_vec(:,1)' gives the \regcharging revenue at 'RrMax_vec(:,2)'
%     tr = linspace(0,0.05,51);
%     for j = 1:1:length(tr)
%         Rr_surf(i,j)= R_reg(tr(j),ts_idx(i));
%         if Rr_surf(i,j) > RrMax_vec(i,1) 
%             RrMax_vec(i,1) = Rr_surf(i,j);
%             RrMax_vec(i,2) = tr(j);
%         end
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----exhausging the best responding ts-----%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the best reponding ts as tr goes from 0 to 0.5;
% tr_idx = linspace(0.04,0.08,41);
% RsMax_vec = zeros(length(tr_idx),2);
% for i=1:1:length(tr_idx)
%     
% % 'ts_vec'gives the theoretical best-responding 'ts'
%     ts_cd1 = t+(p_d-p_a)*theta_bar/c_b; % ts candidate 1
%     ts_cd2 = t+p_d*theta_bar/c_b; % ts candidate 2
%     if ts_cd1/p_d - tr_idx(i)/p_a > 0
%         ts_vec(i) = ts_cd1;
%     elseif ts_cd2/p_d - tr_idx(i)/p_a < 0
%         ts_vec(i) = ts_cd2;
%     else 
%         ts_vec(i) = p_d*tr_idx(i)/p_a;
%     end
% % 'Rs_vec' gives the revenue at the theoretical best-responding 'ts_vec'
%     Rs_vec(i) = R_simple(tr_idx(i),ts_vec(i));
%     
% % exhausive search for the best responding 'ts'
% % 'RsMax_vec(:,2)' gives the exhausive searched best-responding 'ts'
% % 'RsMax_vec(:,1)' gives the \simplecharging revenue at 'RsMax_vec(:,2)'
%     ts = linspace(0.09,0.18,101);
%     for j=1:1:length(ts)
%         
%         RegValid(i,j) = ts(j)/p_d - tr_idx(i)/p_a;
%         Rs_surf(i,j)=R_simple(tr_idx(i),ts(j));
%         if Rs_surf(i,j) > RsMax_vec(i,1)
%             RsMax_vec(i,1) = Rs_surf(i,j);
%             RsMax_vec(i,2) = ts(j);
%         end
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----------the Nash equilibrium---------%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two prices at the Nash equilibrium:
fun_tr_nash = @(tr) exp(-tr*const_r)-exp(-(ts_nash-tr)*const_s)...
        -(exp(-tr*const_r)*const_r+exp(-(ts_nash-tr)*const_s)*const_s)*(tr+Er);
tr_nash = fzero(fun_tr_nash, [0,0.5]);
% Two revenues at the Nash equilibrium:
Rs_nash = R_simple(tr_nash,ts_nash);
Rr_nash = R_reg(tr_nash,ts_nash);
% User welfare at the Nash equilibrium:
% uber_nash = User_welfare( theta_bar, c_b, p_d, p_a, ts_nash, tr_nash );
uber_nash = theta_bar*p_a*(exp(-tr_nash*const_r)-exp(-(ts_nash-tr_nash)*const_s))...
               + theta_bar*p_d*exp(-(ts_nash-tr_nash)*const_s);


% Market share of the two stations at the Nash equilibrium:
As_nash = exp(-(ts_nash-tr_nash)*c_b/((p_d-p_a)*theta_bar));
Ar_nash = exp(-tr_nash*c_b/(p_a*theta_bar)) - As_nash;
% Regulation up and down output
Reg_u_nash = Ar_nash*p_d*x;
Reg_d_nash = Ar_nash*p_d*(1-x);

t_nash(1)=ts_nash;
t_nash(2)=tr_nash;
R_nash(1)=Rs_nash;
R_nash(2)=Rr_nash;
A_nash(1)=As_nash;
A_nash(2)=Ar_nash;
Reg_nash(1)=Reg_u_nash;
Reg_nash(2)=Reg_d_nash;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%--the Monopoly case as a benchmark--%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two prices at the Monopoly case, referring to ISGT-Asia-2015 paper
ts_mono = t + p_d*theta_bar/c_b;
tr_mono = p_a*theta_bar/c_b-Er;

% Reg output at the Monopoly:
As_mono = exp(-(ts_mono-tr_mono)*c_b/((p_d-p_a)*theta_bar));
Ar_mono = exp(-tr_mono*c_b/(p_a*theta_bar)) - As_mono;
Reg_u_mono = Ar_mono*p_d*x;
Reg_d_mono = Ar_mono*p_d*(1-x);
% Revenue at the Monopoly case
R_mono = Ar_mono*c_b*(tr_mono+Er) + As_mono*c_b*(ts_mono - t);
% User welfare at the Monopoly case
% uber_mono = User_welfare( theta_bar, c_b, p_d, p_a, ts_mono, tr_mono );
uber_mono = theta_bar*p_a*(exp(-tr_mono*const_r)-exp(-(ts_mono-tr_mono)*const_s))...
               + theta_bar*p_d*exp(-(ts_mono-tr_mono)*const_s);

t_mono(1)=ts_mono;
t_mono(2)=tr_mono;  
A_mono(1)=As_mono;
A_mono(2)=Ar_mono;
Reg_mono(1)=Reg_u_mono;
Reg_mono(2)=Reg_d_mono;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%---data output to a relay file: HUB.tex---%
%-- then move the data out and add preamble--%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file = 'HUB.tex';
% fid = fopen(['/Users/shuaiwenjing/Dropbox/Wenjing/Competition model/Graphics/' file],'w');
% parameters = ['%$t$=',num2str(t),', $r_d$=',num2str(rd),', $r_u$=',num2str(ru),', $C_B$=',num2str(c_b),', $\rho_d=$',num2str(rho_d),', $\rho_u=$',num2str(rho_u),', $\gamma=$',num2str(gamma),', $\bar\theta=$',num2str(theta_bar),', $x=$',num2str(x)];
% fprintf(fid,'%s\n',parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%---to draw a 3-D vision of \regcharging{} station revenue and T_r^{br}--%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% print the \regcharging{} station revenue on a ts-tr plane
% fprintf(fid,'%s\n','\addplot3 [surf] coordinates{');
% for i=1:1:length(ts_idx)
%     for j=1:1:length(tr)
%      fprintf(fid,'%s','(',num2str(ts_idx(i)),',',num2str(tr(j)),',',num2str(Rr_surf(i,j)),')');
%     end
%     fprintf (fid,'%s\n','');
%     fprintf (fid,'%s\n','');
% end
% fprintf(fid,'%s\n','};');
% 
% % print the regcharging station revenue at best responding tr
% fprintf(fid,'%s\n','\addplot3 coordinates{');
% for i=1:1:length(ts_idx)
%     fprintf(fid,'%s','(',num2str(ts_idx(i)),',',num2str(tr_vec(i)),',',num2str(Rr_vec(i)),')');
% end
% fprintf(fid,'%s\n','};');
% 
% % print the feasible region of tr
% fprintf(fid,'%s\n','\addplot3 coordinates{');
% for i=1:1:length(ts_idx)
%     fprintf(fid,'%s','(',num2str(ts_idx(i)),',',num2str(ts_idx(i)*p_a/p_d),',',num2str(-0.1),')');
% end
% fprintf(fid,'%s\n','};');
% 
% fclose(fid);
end

