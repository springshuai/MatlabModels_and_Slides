function [ Save_vec_pat_1, Save_vec_pat_2, OptK_vec_1, OptK_vec_2 ] = Exhaustive_k_lambda()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eta = 0.90;
epsilon = 0.015;
chi = 0.85;
g = 0.5;
Unplug = 0.2;
A_d = 0.2;
v_0 = 0.1;
alpha = 0.02;
mu = 1;
K_max = 50;
P_m = 200;
% P_m_g = 200;

L_max = 40;
L_min = 10;
L_fix = 30;
nb_Lambda = L_max - L_min +1;
L_vec = linspace (L_min, L_max, nb_Lambda);

Tht_max = 0.25;
Tht_min =0.05;
nb_theta = 21;
Tht_fix = 0.1;
Tht_vec = linspace(Tht_max, Tht_min, nb_theta);


%%%%%%%%%%%% Using Wenjinig's alpha pricing policy %%%%%%%%%%
p_star_wenj = (g-v_0)/(2*alpha);
c_wenj = P_m/p_star_wenj;
function [ p_EV ] = p_EV_wenj(p)
    p_EV = v_0+alpha*p;
end
function [c_e] = cost_elec_wenj (n)
    if n<=c_wenj
        c_e = g*(P_m-n*p_star_wenj) + n*p_star_wenj*p_EV_wenj(p_star_wenj);
    else
        c_e = P_m*p_EV_wenj(P_m/n);
    end
end


%%%%%%%%%%% Using Patrick's epsilon pricing policy %%%%%%%%%%%
p_star_pat = (1-sqrt(v_0/g))/epsilon;
c_pat = P_m/p_star_pat;
function [ p_EV ] = p_EV_pat(p)
    p_EV = v_0/(1-epsilon*p);
end
function [c_e] = cost_elec_pat (n)
    if n<=c_pat
        c_e = g*(P_m-n*p_star_pat) + n*p_star_pat*p_EV_pat(p_star_pat);
    else
        c_e = P_m*p_EV_pat(P_m/n);
    end
end


% %%%%%%%%%%%%%%% the cost analysis using Wenjing's pricing policiy 
% %%%%%%%%%%%%%%%
% 
% Cost_matrix_1 = zeros(nb_Lambda, K_max);
% Cost_matrix_2 = zeros(nb_Lambda, K_max);
% Cost_vec_1 = zeros(2,nb_Lambda); % minimized from Cost_matrix_1, for each lambda
% Cost_vec_2 = zeros(2,nb_Lambda); % minimized from Cost_matrix_2, for each lambda
% Bound_k_1 = zeros(2,nb_Lambda); % upper and lower bound in scheme 1
% Bound_k_2 = zeros(2,nb_Lambda); % upper and lower bound in scheme 2
% 
% for i=1:1:nb_Lambda
%     l = L_vec(i);
%         % initialize the blocking probability
%         Block_1_pre = 1;
%         Block_2_pre = 1;
%     for k=1:1:K_max
%         Park_vec = MMCK_classic( k, l, Tht_fix, p_star_wenj, P_m );
%         d_vec = Steady_state_distribution( k, l, mu, Tht_fix, p_star_wenj, P_m );
%         d_prime_vec = MMCK_pseudo( k, l, mu, Tht_fix, p_star_wenj, P_m );
%         for j = 1:1:k+1
%             Cost_matrix_1(i,k)=Cost_matrix_1(i,k)+d_vec(j)*cost_elec_wenj (j-1);
%             Cost_matrix_2(i,k)=Cost_matrix_2(i,k)+d_prime_vec(j)*cost_elec_wenj (j-1);
%         end
%         Cost_matrix_1(i,k) = Cost_matrix_1(i,k) + k*A_d;
%         workload = l*(Park_vec(k+1)-d_prime_vec(k+1))*Unplug;
%         Cost_matrix_2(i,k) = Cost_matrix_2(i,k) + k*A_d + workload;
%         
% %         % look for the upper and lower bounds 
% %         Block_1 = Park_vec(k+1);
% %         Block_2 = d_prime_vec (k+1);        
% %         % upper bound for the first scheme
% %         if l*(Block_1_pre-Block_1)*P_star/(mu+Tht_fix*P_star)*(g-(v_0+alpha*(P_m/k)))>=A_d
% %             Bound_k_1(1,i)=k;
% %         else
% %         end
% %         % lower bound for the first scheme
% %         if l*(Block_1_pre-Block_1)*(P_m/k)/(mu+Tht_fix*(P_m/k))*(g-(v_0+alpha*P_star))>=A_d
% %             Bound_k_1(2,i)=k;
% %         else
% %         end
% %         % upper bound for the second scheme
% %         if l*(Block_2_pre-Block_2)*P_star/(mu+Tht_fix*P_star)*(g-(v_0+alpha*(P_m/k)))>=A_d
% %             Bound_k_2(1,i)=k;
% %         else
% %         end
% %         % lower bound for the second scheme
% %         if l*(Block_2_pre-Block_2)*(P_m/k)/(mu+Tht_fix*(P_m/k))*(g-(v_0+alpha*P_star))>=A_d
% %             Bound_k_2(2,i)=k;
% %         else
% %         end       
% %         
% %         Block_1_pre = Block_1;
% %         Block_2_pre = Block_2;
%     end
% end
% [Cost_vec_1(1,:),Cost_vec_1(2,:)]=min(Cost_matrix_1,[],2);
% [Cost_vec_2(1,:),Cost_vec_2(2,:)]=min(Cost_matrix_2,[],2);
% Save_vec_wenj_1(1,:) = 100*(1-Cost_vec_1(1,:)/(P_m*g));
% Save_vec_wenj_2(1,:) = 100*(1-Cost_vec_2(1,:)/(P_m*g));


%%%%%%%%%%%%%%%  cost analysis using Patrick's pricing policiy 
%%%%%%%%%%%%%%%

Cost_matrix_1 = zeros(nb_Lambda, K_max);
Cost_matrix_2 = zeros(nb_Lambda, K_max);
Cost_vec_1 = zeros(2,nb_Lambda); % minimized from Cost_matrix_1, for each lambda
Cost_vec_2 = zeros(2,nb_Lambda); % minimized from Cost_matrix_2, for each lambda
Bound_k_1 = zeros(2,nb_Lambda); % upper and lower bound in scheme 1
Bound_k_2 = zeros(2,nb_Lambda); % upper and lower bound in scheme 2

for i=1:1:nb_Lambda
    l = L_vec(i);
        % initialize the blocking probability
        Block_1_pre = 1;
        Block_2_pre = 1;
    for k=1:1:K_max
%        Park_vec = MMCK_classic( k, L_fix, t, p_star_pat, P_m );
        Park_vec = MMKK(k,l,mu);
        d_vec = Steady_state_distribution( k, l, mu, Tht_fix, p_star_pat, P_m );
        d_prime_vec = MMCK_pseudo( k, l, mu, Tht_fix, p_star_pat, P_m );
        for j = 1:1:k+1
            Cost_matrix_1(i,k)=Cost_matrix_1(i,k)+d_vec(j)*cost_elec_pat (j-1);
            Cost_matrix_2(i,k)=Cost_matrix_2(i,k)+d_prime_vec(j)*cost_elec_pat (j-1);
        end
        Cost_matrix_1(i,k) = Cost_matrix_1(i,k) + k*A_d;
        workload = l*(Park_vec(k+1)-d_prime_vec(k+1))*Unplug;
        Cost_matrix_2(i,k) = Cost_matrix_2(i,k) + k*A_d + workload;
        
%         % look for the upper and lower bounds 
%         Block_1 = Park_vec(k+1);
%         Block_2 = d_prime_vec (k+1);        
%         % upper bound for the first scheme
%         if l*(Block_1_pre-Block_1)*P_star/(mu+Tht_fix*P_star)*(g-(v_0+alpha*(P_m/k)))>=A_d
%             Bound_k_1(1,i)=k;
%         else
%         end
%         % lower bound for the first scheme
%         if l*(Block_1_pre-Block_1)*(P_m/k)/(mu+Tht_fix*(P_m/k))*(g-(v_0+alpha*P_star))>=A_d
%             Bound_k_1(2,i)=k;
%         else
%         end
%         % upper bound for the second scheme
%         if l*(Block_2_pre-Block_2)*P_star/(mu+Tht_fix*P_star)*(g-(v_0+alpha*(P_m/k)))>=A_d
%             Bound_k_2(1,i)=k;
%         else
%         end
%         % lower bound for the second scheme
%         if l*(Block_2_pre-Block_2)*(P_m/k)/(mu+Tht_fix*(P_m/k))*(g-(v_0+alpha*P_star))>=A_d
%             Bound_k_2(2,i)=k;
%         else
%         end       
%         
%         Block_1_pre = Block_1;
%         Block_2_pre = Block_2;
    end
end
[Cost_vec_1(1,:),Cost_vec_1(2,:)]=min(Cost_matrix_1,[],2);
[Cost_vec_2(1,:),Cost_vec_2(2,:)]=min(Cost_matrix_2,[],2);
Save_vec_pat_1(1,:) = 100*(1-Cost_vec_1(1,:)/(P_m*g));
Save_vec_pat_2(1,:) = 100*(1-Cost_vec_2(1,:)/(P_m*g));
OptK_vec_1(1,:)=Cost_vec_1(2,:);
OptK_vec_2(1,:)=Cost_vec_2(2,:);





%%%%%%%%%%%%%%%%%%%%% print the TEX files %%%%%%%%%%%%%%%%%%%%%%%%%
%
% first print the saving vs. lambda, with both schemes
% cap_lambda = ['\caption{Cost saved as coming rate increases}'];
% string_parameters=['$\theta$ goes from ' num2str(theta_min) ' to '
% num2str(theta_max) ',$\mu=$' num2str(mu) ',$\lambda=$' num2str(lambda)];
% parameters = [' $\chi$ = 0.85; $g$ = 0.5; $U$ = 0.2; $A_d$ = 0.02; $v_0$ = 0.1; $\alpha$ = 0.02; $\mu$ = 1; $\theta$ = 0.1; P = 10; $P_m$ = 200; '];
 axis_options_lambda = 'legend entries={Saving in Scheme 1,Saving in Scheme 2}';
 x_label_lambda = 'Number of EVs coming per hour';
 y_label_lambda = 'Saving in \%';
 filename_lambda = 'Eco_lambda.tex';
 fid_lambda=fopen(['Data_in_TEX/' filename_lambda],'w');
%  c=clock;
%  infos=['%Obtained from funtion Eco_lambda with parameters of'
% string_parameters ', run on ' int2str(c(3)) '/' int2str(c(2)) '/'
% int2str(c(1)) ' at ' int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6))];
% fprintf(fid_lambda,'%s\n',parameters);
 fprintf(fid_lambda,'%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={(1,1.03)},anchor=south east},width=\figwidth,height=\figheight,cycle list name=\mylist,every axis legend/.append style={nodes={right}},xlabel=' x_label_lambda ',ylabel=' y_label_lambda ',' axis_options_lambda ']']);
 
 fprintf(fid_lambda,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_Lambda
     fprintf(fid_lambda,'%s','(',num2str(L_vec(i)),',',num2str(Save_vec_pat_1(1,i)), ')');
 end
 fprintf(fid_lambda, '%s\n', '};');
 
 fprintf(fid_lambda,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_Lambda
     fprintf(fid_lambda,'%s','(',num2str(L_vec(i)),',',num2str(Save_vec_pat_2(1,i)), ')');
 end
 fprintf(fid_lambda, '%s\n', '};');
 
 fprintf(fid_lambda,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_Lambda
     fprintf(fid_lambda,'%s','(',num2str(L_vec(i)),',',num2str(OptK_vec_1(1,i)), ')');
 end
 fprintf(fid_lambda, '%s\n', '};');
 
 fprintf(fid_lambda,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_Lambda
     fprintf(fid_lambda,'%s','(',num2str(L_vec(i)),',',num2str(OptK_vec_2(1,i)), ')');
 end
 fprintf(fid_lambda, '%s\n', '};');
 
% fprintf(fid_lambda, '%s\n', '};');
% 
 fprintf(fid_lambda,'%s\n','\end{axis}\end{tikzpicture}}');
% fprintf(fid_lambda,'\n%s\n',cap_lambda);
% fclose(fid_lambda);
% end
% 
% print the upper-lower and actual optimum of k for scheme 1
% cap_bound_1 = [' \caption{Upper-and Lower-bounds of $k_1^*$}'];
% axis_options_k_1 = 'legend entries={Upper bound, Lower bound, Real value}';
% x_label_k_1 = 

end

