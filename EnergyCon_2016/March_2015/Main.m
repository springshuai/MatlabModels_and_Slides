function [ Cost_compare ] = Main(unPlug, BC, BL)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
load('Variation_K.mat','BatteryCapital','BatteryLifespan','EUI','S','alpha','g','v_0');

% BatteryCapital is the oneshot cost containing parcharsing and instalment,
% the unit is euro per kWh, BatteryLifespan is the life of the battery in
% terms of cycles, so CKC is the Cost per Kilowatt per cycle. 
 CKC = BC/BL;


% P is the optimal discharging point considering the electricity price of
% the grid and the pricing function of EVs
P = (g-v_0)/(2*alpha);

% Maximum enery comsumpiton per hour is the Energy Use Index multiply by
% the area of the mall
P_m = EUI*S;
P_m = max(P_m,P);
% Maximum simultanuous dishcharging number of cars
c = floor(P_m/P);
K_m = 50;
lambda = 30;
mu = 1;
theta = 0.05;
A_d = 0.1;
k_vector=zeros(1,K_m);
global string_parameters;
string_parameters = ['$\lambda = ' num2str(lambda) ',\mu =' num2str(mu) ' ,\theta = ' num2str(theta) ',A_d = ' num2str(A_d) ',U = ' num2str(unPlug) ',C_B = ' num2str(CKC) '$'];
% the vectors Pr, Output, Complements and Cost are the probabilitiy, EVs
% output, Grid compelent, and average cost respectivly. the each contains
% five lines corresponding to case 1 to 3, and 5 which is grid only, as a
% benchmark for comparison, the fifth line is the cost of using battery.

% Case 1 stands for without removing
% Case 2 means with removoing
% Cost_case_1 = zeros (1,K_m);
% Cost_case_2 = zeros (1,K_m);
% Cost_case_3 = zeros (1,K_m);
% Cost_case_Battery = zeros (1,K_m);
% Cost_case_Grid = zeros (1,K_m);

for k = 1:1:K_m;
    Pr = zeros(2,k+1);
    Cost = zeros(4,1);
    Output = zeros(2,k+1);
    Complement = zeros(2,k+1);

    % the steady state probability distribution of a MMKK queue
    % Pr(1,:) = MMKK( k, lambda, mu);

    % the steady state probability distribution of a MMCK queue
    Pr(2,:) = MMCK_pseudo( k, lambda, mu, theta, P, P_m );

    % the steady state probability distribution of a mingled MMKK and MMCK queue
    Pr(1,:) = Steady_state_distribution(k, lambda, mu, theta, P, P_m );

    % Battery supporing the mall
%     Cost_case_Battery(1,:) = P_m*(CKC+v_0);
    % without any EV nor battery
    Cost_case_Grid(1,:) = P_m*g;

    % 'j' is the index of cases and 'i' is the number of discharging cars
    for j = 1:1:2
        Cost(j,1) = 0;
        for i=0:1:k
        if i<c+1
            Output (j,i+1) = i*P;
            Complement (j,i+1) = P_m - Output(j,i+1);
            Cost(j,1) = Cost(j,1) + (Output (j,i+1)*(alpha*P+v_0) + Complement (j,i+1)*g)*Pr(j,i+1);
        else
            Output (j,i+1) = P_m;
            Complement (j,i+1) = P_m - Output(j,i+1);
            Cost(j,1) = Cost(j,1) + (Output (j,i+1)*(alpha*(P_m/(i+1))+v_0) + Complement (j,i+1)*g)*Pr(j,i+1);
        end
        end
    end
%     prob = MMCK_classic (k, lambda, theta, P, P_m);
%     block = prob (1,k+1);
    Pr_ref = MMKK(k, lambda, mu);
    NeedPush = Pr_ref(1,k+1)-Pr(2,k+1); 
    workload = lambda*NeedPush*unPlug;

    Cost_case_1(1,k)=Cost(1,1)+k*A_d;
    Cost_case_2(1,k)=Cost(2,1)+k*A_d;
    Cost_case_3(1,k)=Cost(2,1)+k*A_d+workload;
    
    k_vector(k)=k;
end
% Cost_compare = zeros(5,K_m);
% Cost_compare(1,:) = Cost_case_1;
% Cost_compare(2,:) = Cost_case_2;
% Cost_compare(3,:) = Cost_case_3;
% Cost_compare(4,:) = Cost_case_Battery;
% Cost_compare(5,:) = Cost_case_Grid;

Cost_compare = zeros(3,K_m);
Cost_compare(1,:) = Cost_case_1;
Cost_compare(2,:) = Cost_case_3;
Cost_compare(3,:) = Cost_case_Grid;

axis_options = 'legend entries={No unplugging,With unplugging,Grid only}';
x_label = 'Number of spots $k$';
y_label = 'Cost per hour';
filename = 'Variation_K.tex';
fid=fopen(['Data_in_TEX/' filename],'w');
c=clock;
infos=['%Obtained from function Main with parameters of' string_parameters ', run on ' int2str(c(3)) '/' int2str(c(2)) '/' int2str(c(1)) ' at ' int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6))];
fprintf(fid,'%s\n',infos);
fprintf(fid,'%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={(1,1.03)},anchor=south east},width=\figwidth,height=\figheight,cycle list name=\mylist,every axis legend/.append style={nodes={right}},xlabel=' x_label ',ylabel=' y_label ',' axis_options ']']);

for i=1:1:size(Cost_compare,1)
    fprintf(fid,'%s\n','\addplot coordinates{');
    for j = 1:1:size(Cost_compare,2)
        fprintf(fid,'%s','(',num2str(j),',',num2str(Cost_compare(i,j)), ')');
    end
    fprintf(fid, '%s\n', '};');
end
fprintf(fid,'%s\n','\end{axis}\end{tikzpicture}}');
fprintf(fid,'%s\n',string_parameters);
fclose(fid);
end
