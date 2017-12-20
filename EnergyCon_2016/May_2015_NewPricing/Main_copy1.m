function [ Cost_case_1 ] = Main_copy1(lambda, ok)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
load('Variation_K.mat','BatteryCapital','BatteryLifespan','EUI','S','alpha','g','v_0');

% BatteryCapital is the oneshot cost containing parcharsing and instalment,
% the unit is euro per kWh, BatteryLifespan is the life of the battery in
% terms of cycles, so CKC is the Cost per Kilowatt per cycle. 
CKC = BatteryCapital/BatteryLifespan;
unPlug = 0.2;

% P is the optimal discharging point considering the electricity price of
% the grid and the pricing function of EVs
P = (g-v_0)/(2*alpha);

% Maximum enery comsumpiton per hour is the Energy Use Index multiply by
% the area of the mall
P_m = EUI*S;
P_m = max(P_m,P);
% Maximum simultanuous dishcharging number of cars
c = floor(P_m/P);
K_m = 5;
% lambda = 20;
mu = 1;
theta = 0.2;
Ad = 0.2;

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

for count = 1:1:K_m;
    k = count + ok - 3;
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
    Cost_case_Battery(1,:) = P_m*(CKC+v_0);
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
%     NeedPush = Pr(1,k+1); 
%     workload = lambda*NeedPush*unPlug;

    Cost_case_1(1,count)=Cost(1,1)+k*Ad;
%     Cost_case_2(1,k)=Cost(2,1)+k*Ad;
%     Cost_case_3(1,k)=Cost(2,1)+k*Ad+workload;
end
% Cost_compare = zeros(5,K_m);
% Cost_compare(1,:) = Cost_case_1;
% Cost_compare(2,:) = Cost_case_2;
% Cost_compare(3,:) = Cost_case_3;
% Cost_compare(4,:) = Cost_case_Battery;
% Cost_compare(5,:) = Cost_case_Grid;
end

