function [ ratio_R,ratio_U,log_R_mono,log_R_nash,Uber_vec,Uber_mono,...
    Alpha_s_vec,Alpha_r_vec,Alpha_s_mono,Alpha_r_mono] = Nash_over_Mono(  )
%NASH_OVER_MONO Summary of this function goes here
%   Detailed explanation goes here

theta_bar = linspace(0.1,0.5,41);
ru=1.7;
rd=0.6;

[ x_opt_vec,Rs_vec,Rr_vec,Ts_vec,Tr_vec,Reg_up_vec,Reg_down_vec,...
         Alpha_s_vec,Alpha_r_vec,Uber_vec,Er_vec,Er_min_vec,thr_Er_nash_vec] = MainFunction_Nash( ru, rd );
     
[ x_opt_mono,Rs_mono,Rr_mono,Ts_mono,Tr_mono,Reg_up_mono,Reg_down_mono,...
         Alpha_s_mono,Alpha_r_mono,Er_mono,Uber_mono,R_sim] = MainFunction_Mono( ru, rd );
 
ratio_R = (Rr_vec+Rs_vec)./(Rr_mono+Rs_mono);
ratio_U = Uber_vec./Uber_mono;
ratio_Tr = Tr_vec./Tr_mono;
ratio_Ts = Ts_vec./Ts_mono;
ratio_social = (Rr_vec+Rs_vec+Uber_vec)./(Rr_mono+Rs_mono+Uber_mono);

log_R_mono = log(Rr_mono+Rs_mono)/transpose(log(theta_bar));
log_R_nash = log(Rr_vec+Rs_vec)/transpose(log(theta_bar));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % write the data into HUB_nash_over_Mono.tex
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun_name = 'HUB_nash_over_Mono';
file_name = 'HUB_nash_over_Mono.tex';
file_path = '/Users/shuaiwenjing/Dropbox/Wenjing/Competition model/Graphics/' ;
fid = fopen([file_path file_name],'w');
% fprintf(fid,'%s\n',['% from' fun_name '; parameters: t = ' num2str(t) '; c_b = '...
%        num2str(c_b) '; rho_d = rho_u = ' num2str(rho_d) '; gamma = ' num2str(gamma)...
%        '; rd = ' num2str(rd) '; ru = ' num2str(ru)]);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- revenue ratio nash_over_Mono
fprintf(fid,'%s\n','% revenue ratio nash_over_Mono');
fprintf(fid,'%s\n','\addplot coordinates{');
for tht_id = 1:1:length(theta_bar)
    fprintf(fid,'%s','(',num2str(theta_bar(tht_id)),',',...
            num2str(ratio_R(tht_id)),')');
end
fprintf(fid,'%s\n', '};');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- user welfare ratio nash_over_Mono
fprintf(fid,'%s\n','% user welfare ratio nash_over_Mono');
fprintf(fid,'%s\n','\addplot coordinates{');
for tht_id = 1:1:length(theta_bar)
    fprintf(fid,'%s','(',num2str(theta_bar(tht_id)),',',...
            num2str(ratio_U(tht_id)),')');
end
fprintf(fid,'%s\n', '};');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- social welfare ratio nash_over_Mono
fprintf(fid,'%s\n','% social welfare ratio nash_over_Mono');
fprintf(fid,'%s\n','\addplot coordinates{');
for tht_id = 1:1:length(theta_bar)
    fprintf(fid,'%s','(',num2str(theta_bar(tht_id)),',',...
            num2str(ratio_social(tht_id)),')');
end
fprintf(fid,'%s\n', '};');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- R_charging price ratio nash_over_Mono
fprintf(fid,'%s\n','% R_charging price ratio nash_over_Mono');
fprintf(fid,'%s\n','\addplot coordinates{');
for tht_id = 1:1:length(theta_bar)
    fprintf(fid,'%s','(',num2str(theta_bar(tht_id)),',',...
            num2str(ratio_Tr(tht_id)),')');
end
fprintf(fid,'%s\n', '};');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %%%%%%%%------------- S_charging price ratio nash_over_Mono
fprintf(fid,'%s\n','% S_charging price ratio nash_over_Mono');
fprintf(fid,'%s\n','\addplot coordinates{');
for tht_id = 1:1:length(theta_bar)
    fprintf(fid,'%s','(',num2str(theta_bar(tht_id)),',',...
            num2str(ratio_Ts(tht_id)),')');
end
fprintf(fid,'%s\n', '};');




fclose(fid);
end

