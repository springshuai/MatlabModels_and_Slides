function [ Save_vec_M_1, Save_vec_M_2, Save_vec_M_3, OptK_vec_M_1, OptK_vec_M_2, OptK_vec_M_3, Cost_matrix_1, Cost_matrix_2, Cost_matrix_3 ] = Exhaustive_k_theta_Myopic()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

run Parameters.m

L_fix = 20;

Tht_max = 0.25;
Tht_min =0.05;
nb_theta = 21;
Tht_vec = linspace(Tht_max, Tht_min, nb_theta);

% %%%%%%%%%%%% Using Wenjinig's alpha pricing policy %%%%%%%%%%
% p_star_wenj = (g-v_0)/(2*alpha);
% c_wenj = P_m/p_star_wenj;
% function [ p_EV ] = p_EV_wenj(p)
%     p_EV = v_0+alpha*p;
% end
% function [c_e] = cost_elec_wenj (n)
%     if n<=c_wenj
%         c_e = g*(P_m-n*p_star_wenj) + n*p_star_wenj*p_EV_wenj(p_star_wenj);
%     else
%         c_e = P_m*p_EV_wenj(P_m/n);
%     end
% end
% 
% 
% %%%%%%%%%%% Using Patrick's epsilon pricing policy %%%%%%%%%%%
% p_star_pat = (1-sqrt(v_0/g))/epsilon;
% c_pat = P_m/p_star_pat;
% function [ p_EV ] = p_EV_pat(p)
%     p_EV = v_0/(1-epsilon*p);
% end
% function [c_e] = cost_elec_pat (n)
%     if n<=c_pat
%         c_e = g*(P_m-n*p_star_pat) + n*p_star_pat*p_EV_pat(p_star_pat);
%     else
%         c_e = P_m*p_EV_pat(P_m/n);
%     end
% end


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
c_vec_M = zeros(1,nb_theta);
Cost_matrix_1 = zeros(nb_theta, K_max);
Cost_matrix_2 = zeros(nb_theta, K_max);
Cost_matrix_3 = zeros(nb_theta, K_max);
Cost_vec_1 = zeros(2,nb_theta); % minimized from Cost_matrix_1, for each lambda
Cost_vec_2 = zeros(2,nb_theta); % minimized from Cost_matrix_2, for each lambda
Cost_vec_3 = zeros(2,nb_theta); % minimized from Cost_matrix_2, for each lambda
Bound_k_1 = zeros(2,nb_theta); % upper and lower bound in scheme 1
Bound_k_2 = zeros(2,nb_theta); % upper and lower bound in scheme 2

for i=1:1:nb_theta
    t = Tht_vec(i);
        % initialize the blocking probability
        Block_1_pre = 1;
        Block_2_pre = 1;
    for k=1:1:K_max
%        Park_vec = MMCK_classic( k, L_fix, t, p_star_pat, P_m );
        [ P0_vec, Cost_ins, c_vec_M(i) ] = P0_vec_Myopic( k, g, v_0, epsilon, P_m );
        
        Park_vec = MMKK(k,L_fix,mu);
        d_vec = Steady_state_distribution( k, L_fix, mu, t, P0_vec(1,1), P_m );
        d_prime_vec = MMCK_pseudo( k, L_fix, mu, t, P0_vec(1,1), P_m );
        
%         for j = 1:1:k+1
%             Cost_matrix_1(i,k)=Cost_matrix_1(i,k)+d_vec(j)*cost_elec_pat (j-1);
%             Cost_matrix_2(i,k)=Cost_matrix_2(i,k)+d_prime_vec(j)*cost_elec_pat (j-1);
%         end
        Cost_matrix_1(i,k)=dot(Cost_ins,d_vec);
        Cost_matrix_2(i,k)=dot(Cost_ins,d_prime_vec);
        Cost_matrix_3(i,k) = Cost_matrix_2(i,k) + k*A_d;
        Cost_matrix_1(i,k) = Cost_matrix_1(i,k) + k*A_d;
        workload = L_fix*(Park_vec(k+1)-d_prime_vec(k+1))*Unplug;
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
[Cost_vec_3(1,:),Cost_vec_3(2,:)]=min(Cost_matrix_3,[],2);
Save_vec_M_1(1,:) = 100*(1-Cost_vec_1(1,:)/(P_m*g));
Save_vec_M_2(1,:) = 100*(1-Cost_vec_2(1,:)/(P_m*g));
Save_vec_M_3(1,:) = 100*(1-Cost_vec_3(1,:)/(P_m*g));
OptK_vec_M_1(1,:)=Cost_vec_1(2,:);
OptK_vec_M_2(1,:)=Cost_vec_2(2,:);
OptK_vec_M_3(1,:)=Cost_vec_3(2,:);





%%%%%%%%%%%%%%%%%%%%% print the TEX files %%%%%%%%%%%%%%%%%%%%%%%%%
%
 axis_options_theta = 'legend entries={Saving in Scheme 1,Saving in Scheme 2,Saving in Scheme 2 with no unplugging cost}';
 x_label_theta = 'Surplus energy in one EV';
 y_label_theta = 'Saving in \%';
 filename = 'Eco_theta_M.tex';
 fid=fopen(['Data_in_TEX/' filename],'w');
 time =clock;
 infos = ['% Obtained from Exhaustive_k_theta_Myopic(), run on ', int2str(time(3)), '/', int2str(time(2)), '/',int2str(time(1)), ' at ', int2str(time(4)), ':', int2str(time(5)), ':', int2str(time(6)), 'when $\lambda$ = ', num2str(L_fix)]; 
 fprintf(fid,'%s\n',infos);
 fprintf(fid,'%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={(1,1.03)},anchor=south east},width=\figwidth,height=\figheight,cycle list name=\mylist,every axis legend/.append style={nodes={right}},xlabel=' x_label_theta ',ylabel=' y_label_theta ',' axis_options_theta ']']);
 
 fprintf(fid, '%s\n', '% Saving in Scheme 1');
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(Save_vec_M_1(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');
 
 fprintf(fid, '%s\n', '% Saving in Scheme 2'); 
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(Save_vec_M_2(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');

 fprintf(fid, '%s\n', '% Saving in Scheme 2 without unplugging cost');
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(Save_vec_M_3(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');
 
 fprintf(fid,'%s\n', '% Optimal k in scheme 1');
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(OptK_vec_M_1(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');
 
 fprintf(fid,'%s\n', '% Optimal k in scheme 2');
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(OptK_vec_M_2(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');
 
 fprintf(fid,'%s\n', '% Optimal k in scheme 2 without unplugging cost');
 fprintf(fid,'%s\n','\addplot coordinates{');
 for i = 1:1:nb_theta
     fprintf(fid,'%s','(',num2str(1/Tht_vec(i)),',',num2str(OptK_vec_M_3(1,i)), ')');
 end
 fprintf(fid, '%s\n', '};');
 
 fprintf(fid,'%s\n','\end{axis}\end{tikzpicture}}');

 fclose(fid);


end

