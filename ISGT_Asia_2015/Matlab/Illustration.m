function [Tr_theo_rd_ru,R_theo_rd_ru,U_bar_theo] = Illustration(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Probabilities;
rho_u = 0.5;% the probability of receiving an regulation "up" demand 
rho_d = 0.5;%..."down"...
rho_n = 1-rho_u-rho_d; %..."null"..., 

gamma = 0;

t=0.03;
r_u_min=2-rho_u+gamma*rho_u^(-0.5)*(1-rho_u)^1.5;
r_d_min=1-(rho_d-gamma*sqrt(rho_d-(rho_d^2)));
% rd_Step = 0.1*r_d_min;
% ru_Step = 0.1*r_u_min;

rd_Step = 0.1;
ru_Step = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%-----print the experimantal revenue based in .tex----%%%%%%
% file_convex = 'ShowConvex.tex';
% fid_convex = fopen (['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' file_convex], 'w');
% x_convex = 'x';
% y_convex = 'Expected revenue brought per EV';
% fprintf(fid_convex, '%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
%     x_convex ',ylabel=', y_convex ']']);
% 
% file_welfare = 'UserWelfare.tex';
% fid_welfare = fopen (['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' file_welfare], 'w');
% x_welfare = 'x';
% y_welfare = 'User welfare';
% fprintf(fid_welfare, '%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
%     x_welfare ',ylabel=', y_welfare ']']);

% file_P0_Pn = 'P0_Pn.tex';
% fid_P0_Pn = fopen (['/Users/shuaiwenjing/Dropbox/Wenjing/August_ISGT_asia_2015/Graphics/' file_P0_Pn], 'w');
% x_P0_Pn = '$\r_u$';
% y_P0_Pn = '$\r_d$';
% fprintf(fid_welfare, '%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={1,1},anchor=south west},width=\figwidth,height=\figheight,cycle list name=\mylist,xlabel=' ...
%     x_P0_Pn ',ylabel=',y_P0_Pn ']']);


% for cbn_rd_ru = 0:1:3
%     switch cbn_rd_ru
%         case 0
%             rd_Start = 0;
%             rd_End = 0;
%             rd_Num = length(rd_Start:rd_Step:rd_End);
% 
%             ru_Start = 0;
%             ru_End = 0;
%             ru_Num = length(ru_Start:ru_Step:ru_End);
%         case 1
%             rd_Start = 1.5*r_d_min;
%             rd_End = 1.5*r_d_min;
%             rd_Num = length(rd_Start:rd_Step:rd_End);
% 
%             ru_Start = 0;
%             ru_End = 0;
%             ru_Num = length(ru_Start:ru_Step:ru_End);
%         case 2
%             rd_Start = 0;
%             rd_End = 0;
%             rd_Num = length(rd_Start:rd_Step:rd_End);
% 
%             ru_Start = 1.2*r_u_min;
%             ru_End = 1.2*r_u_min;
%             ru_Num = length(ru_Start:ru_Step:ru_End);
%         case 3
%             rd_Start = 1.4*r_d_min;
%             rd_End = 1.6*r_d_min;
%             rd_Num = length(rd_Start:rd_Step:rd_End);
% 
%             ru_Start = 1.1*r_u_min;
%             ru_End = 1.3*r_u_min;
%             ru_Num = length(ru_Start:ru_Step:ru_End);
%     end
            rd_Start = 0.5;
            rd_End =0.8;
            rd_Num = length(rd_Start:rd_Step:rd_End);

            ru_Start = 1.5;
            ru_End = 2.1;
            ru_Num = length(ru_Start:ru_Step:ru_End);
%             
% P_x = zeros(rd_Num,ru_Num);
    for d=1:1:rd_Num
    r_d(d)=rd_Start+(d-1)*rd_Step;
    
        for u=1:1:ru_Num
        r_u(u)=ru_Start+(u-1)*ru_Step;
        [ R_compare, slope, fun1, fun2, U_bar, w_P_vec] = R_Tc_Tr_Correction( t, gamma, r_u(u), r_d(d), rho_u, rho_d);
        
%         Tc_theo_rd_ru(:,d,u)=R_compare(1,1,:);
        Tr_theo_rd_ru(:,d,u)=R_compare(1,2,:);
        R_theo_rd_ru(:,d,u)=R_compare(1,3,:);
        U_bar_theo(:,d,u)=U_bar(1,:);
%         subplot(3,2,1);plot(0:0.05:1,transpose(R_theo_rd_ru(:,d,u))); hold on
%         subplot(3,2,3);plot(0:0.1:1,transpose(Tr_theo_rd_ru(:,d,u))); hold on
%         subplot(3,2,5);plot(0:0.05:1,transpose(U_bar_theo(:,d,u)));hold on
%         if (r_d(d)>=r_d_min)||(r_u(u)>=r_u_min)
%          if R_theo_rd_ru(1,d,u)<R_theo_rd_ru(length(w_P_vec),d,u);
%              P_x(d,u)=1;
%          else
%              P_x(d,u)=0;
%          end
%         else
%              P_x(d,u)=NaN;
%         end
%              
            
%         %%%%%%% ------Start printing the revenues
%         fprintf(fid_convex, '%s\n', ['% revenue; r_d = ' num2str(r_d(d)) 'r_u = ' num2str(r_u(u))]);
%         fprintf(fid_convex, '%s\n', '\addplot coordinates{');
%         for id_x = 1:1:length(w_P_vec);
%             fprintf(fid_convex, '%s', '(',num2str(w_P_vec(id_x)),',',num2str(R_theo_rd_ru(id_x,d,u)),')');
%         end
%         fprintf(fid_convex, '%s\n', '};');
%         
%         %%%%%%% ------Start printing the User welfare
%         fprintf(fid_welfare, '%s\n', ['% welfare;% r_d = ' num2str(r_d(d)) 'r_u = ' num2str(r_u(u))]);
%         fprintf(fid_welfare, '%s\n', '\addplot coordinates{');
%         for id_x = 1:1:length(w_P_vec);
%             fprintf(fid_welfare, '%s', '(',num2str(w_P_vec(id_x)),',',num2str(U_bar_theo(id_x,d,u)),')');
%         end
%         fprintf(fid_welfare, '%s\n', '};');
%         %%%%%%%%%%%----write in P0_Pn.tex-----%%%%%%%%%%%
%         fprintf(fid_P0_Pn,'%s\n',['(',num2str(r_u(u)),',',num2str(r_d(d)),')','[',num2str(P_x(d,u)),']']);

%         Tc_bench_rd_ru(:,d,u)=R_compare(2,1,:);
%         Tr_bench_rd_ru(:,d,u)=R_compare(2,2,:);        
%         R_bench_rd_ru(:,d,u)=R_compare(2,3,:);
%         U_bar_bench(:,d,u)=U_bar(2,:);
%         subplot(3,2,5);plot(0:0.1:1,transpose(U_bar_bench(:,d,u)));hold on
        
%         Tc_real_rd_ru(:,d,u)=R_compare(3,1,:);
%         Tr_real_rd_ru(:,d,u)=R_compare(3,2,:);
%         R_real_rd_ru(:,d,u)=R_compare(3,3,:);
%         U_bar_real(:,d,u)=U_bar(3,:);
        
%         subplot(3,2,2);plot(0:0.1:1,transpose(R_real_rd_ru(:,d,u))); hold on
%         subplot(3,2,4);plot(0:0.1:1,transpose(Tr_real_rd_ru(:,d,u))); hold on
%         subplot(3,2,6);plot(0:0.1:1,transpose(U_bar_real(:,d,u)));hold on
%         subplot(3,2,6);plot(0:0.1:1,transpose(U_bar_bench(:,d,u)));hold on
        
        end

    end
% end
% hold off

% fprintf (fid_convex, '%s\n', '\end{axis}\end{tikzpicture}}');
% fprintf (fid_convex, '%s\n','\caption{Aggregator Revenue with multiple combinations of $r_d$ and $r_u$}');
% fclose (fid_convex);
% 
% fprintf (fid_welfare, '%s\n', '\end{axis}\end{tikzpicture}}');
% fprintf (fid_welfare, '%s\n','\caption{User welfare with multiple combinations of $r_d$ and $r_u$}');
% fclose (fid_welfare);
% fprintf(fid_P0_Pn,'%s\n','};\end{axis}\end{tikzpicture}');
% fclose(fid_P0_Pn);

end

