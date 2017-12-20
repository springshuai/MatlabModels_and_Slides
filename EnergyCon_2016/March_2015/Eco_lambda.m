function [ Save, opt_K_U, opt_K_L, opt_K_R ] = Eco_lambda(mu, alpha, P_m, A_d, g, v_0 )
% this function finds the optimal number of spots for each lambda, and the
% saving made by those spots

%   load('Variation_K.mat','BatteryCapital','BatteryLifespan','EUI','S','alpha','g','v_0');

K = 40; % a sufficiently large number, larger than possible optimal k
P_star = (g-v_0)/(2*alpha); % P_star is the optimal discharging rate
c = P_m/P_star; % c is the minimum number of discharging EVs that can enebles the mall to cut-off from the grid
v = alpha*P_star + v_0; % v is the unit electricity price offered to 
Pr_Park = zeros(3,K+1); % Pr(k+1) is the number of 'k' cars parked in the parking lot
unPlug = 0.02; % unplug is the cost of unplugging one car
Block = ones (3,K+1); % Block(k+1) is the blocking probability when having 'k' spots installe
% Gain = zeros(1,K); % Gain(k) is the marginal saving made by a working 'k-th' spot
theta = 0.1;

RealCost_min = P_m*g;
lambda_max = 40;
lambda_min = 10;
lambda_vec = linspace(lambda_min, lambda_max, 31);
nb_lambda = size(lambda_vec,2);
Save = zeros(3,nb_lambda); % saving for nb_lambda values of lambda, in three cases
opt_K_U = zeros(3,nb_lambda);% the upper-bound of optimal number of spots for nb_lambda values of lambda, in three cases
opt_K_L = zeros(3,nb_lambda);% the upper-bound of optimal number of spots for nb_lambda values of lambda, in three cases
opt_K_R = zeros(3,nb_lambda);% the real optimal number of spots for nb_lambda values of lambda, in three cases


for idx = 1:1:nb_lambda
    lambda = lambda_vec(idx);
    Pr_Park(1,:) = MMKK(K, lambda, mu);% for the blocking probability 
    Pr_Park(2,:) = MMCK_pseudo(K, lambda, mu, theta, P_star, P_m );
    
%%%%%%%%%%%%%%%%%%%%%%%%-----------First Scheme, without removing -----------%%%%%%%%%%%%%   
 
    ca = 1;
        Marginal_gain_U = zeros (1,K);
        Marginal_gain_L = zeros (1,K);
        Gain_max = zeros(3,1);
        for k = 0:1:K
            Block(ca,k+1) = Pr_Park(ca,k+1)/sum(Pr_Park(ca,1:1:k+1));% Block(k+1) is the blocking prob when having k spots
        end
% find the upper- and lower- bounds of optmimal k
        for k = 1:1:K
            if k<=c
                Marginal_gain_U(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-v);
                Marginal_gain_L(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-v);
            else
                Marginal_gain_U(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-(v_0+alpha*P_m/k));
                Marginal_gain_L(k) = lambda/(mu+theta*(P_m/k))*(Block(ca,k) - Block(ca,k+1))*(P_m/k)*(g-(g+v_0)/2);
            end

            if A_d <= Marginal_gain_U(k)
               opt_K_U(ca,idx) = k;
            end
            if A_d <= Marginal_gain_L(k)
               opt_K_L(ca,idx) = k;
               Gain_max(ca,1) = Gain_max(ca,1) + (Marginal_gain_L(k)-A_d);
            else
            end
        end
% find the optimal k between the upper- and lower- bounds         
        if opt_K_U(ca,idx)~=opt_K_L(ca,idx)
            opt_K_U_L_gap = opt_K_U(ca,idx)-opt_K_L(ca,idx)+1;
            RealCost = zeros(opt_K_U_L_gap,1);
            RealCost_min = P_m*g;
            for i = opt_K_L(ca,idx):1:opt_K_U(ca,idx);% 'i' is the candidate of the optimal 'k'
                Pr_Disch = Steady_state_distribution( i, lambda, mu, theta, P_star, P_m );
                for j = 1:1:i+1
                    RealCost(i-opt_K_L(ca,idx)+1,1) = RealCost(i-opt_K_L(ca,idx)+1,1)+Pr_Disch(j)*Cost( j-1, c, g, v_0, P_m, alpha);
                end
                RealCost(i-opt_K_L(ca,idx)+1,1)=RealCost(i-opt_K_L(ca,idx)+1,1)+i*A_d;
                if RealCost(i-opt_K_L(ca,idx)+1,1)<RealCost_min
                    RealCost_min = RealCost(i-opt_K_L(ca,idx)+1,1);
                    opt_K_R(ca,idx) = i;
                else
                end
            end
        else
            opt_K_R(ca,idx)=opt_K_L(ca,idx);
        end

        if opt_K_R(ca,idx)<=c
            Save (ca,idx) = Gain_max(ca,1)/(P_m*g*0.85);
        else
            Save (ca,idx)= RealCost_min/(P_m*g*0.85);
            Save (ca,idx)=1-Save (ca,idx);
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%-----------Second Scheme, free removing -----------%%%%%%%%%%%%%   

%         ca = 2;
%         Marginal_gain = zeros (1,K);
%         for k = 0:1:K
%             Block(ca,k+1) = Pr(ca,k+1)/sum(Pr(ca,1:1:k+1));% Block(k+1) is the blocking prob when having k spots
%         end
%         for k = 1:1:K
%             if k<=c
%                 Marginal_gain(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*Gain(1,k);
%                 %Marginal_gain(k) = (Block(ca,k) - Block(ca,k+1))*Gain(1,k);
%             else
%                 Marginal_gain(k) = lambda/(mu+theta*P_m/k)*(Block(ca,k) - Block(ca,k+1))*Gain(1,k);
%                 %Marginal_gain(k) = (Block(ca,k) - Block(ca,k+1))*Gain(1,k);
%             end
% 
%             if A_d <= Marginal_gain(k)
%                opt_K(ca,lambda-L_min+1) = k;
%                Gain_max(ca,1) = Gain_max(ca,1) + (Marginal_gain(k)-A_d);
%             end
%         end
%         
%         check(lambda-L_min+1,2,:) = Main_copy2(lambda, opt_K(2,lambda-L_min+1));
%         Save (ca,lambda-L_min+1)= Gain_max(ca,1)/(P_m*g);

    
%%%%%%%%%%%%%%%%%%%%%%%%-----------Third Scheme, with removing cost-----------%%%%%%%%%%%%%   
 
        ca = 3;
        Marginal_gain_U = zeros (1,K);
        Marginal_gain_L = zeros (1,K);
        Gain_max = zeros(3,1);
        for k = 0:1:K
            Block(3,k+1) = Pr_Park(2,k+1)/sum(Pr_Park(2,1:1:k+1));% Block(k+1) is the blocking prob when having k spots
        end
% find the upper- and lower- bounds of optmimal k
        for k = 1:1:K
            Pr_ref = MMCK_pseudo( k, lambda, mu, theta, P_star, P_m );
            NeedPush = Block(1,k+1)-Pr_ref(1,k+1); % probability that bloking happens and not caused by discharging EVs
            workload = lambda*NeedPush*unPlug;
            if k<=c
                Marginal_gain_U(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-v) - workload;
                Marginal_gain_L(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-v)- workload;
            else
                Marginal_gain_U(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-(v_0+alpha*P_m/k))- workload;
                Marginal_gain_L(k) = lambda/(mu+theta*(P_m/k))*(Block(ca,k) - Block(ca,k+1))*(P_m/k)*(g-(g+v_0)/2)- workload;
%                 Marginal_gain_U(k) = lambda/(mu+theta*P_m/k)*(Block(ca,k) - Block(ca,k+1))*P_star*(g-v)- workload;
%                 Marginal_gain_L(k) = lambda/(mu+theta*P_star)*(Block(ca,k) - Block(ca,k+1))*(alpha*P_m^2*(1/(k-1) - 1/k))- workload;
            end
            
            if A_d <= Marginal_gain_U(k)
               opt_K_U(ca,idx) = k;
            end
            if A_d <= Marginal_gain_L(k)
               opt_K_L(ca,idx) = k;
               Gain_max(ca,1) = Gain_max(ca,1) + (Marginal_gain_L(k)-A_d);
            else
            end
        end
% find the optimal k between the upper- and lower- bounds        
        if opt_K_U(ca,idx)~=opt_K_L(ca,idx)
            opt_K_U_L_gap = opt_K_U(ca,idx)-opt_K_L(ca,idx)+1;
            RealCost = zeros(opt_K_U_L_gap,1);
            RealCost_min = P_m*g;
            for i = opt_K_L(ca,idx):1:opt_K_U(ca,idx);% 'i' is the candidate of the optimal 'k'
                Pr_Disch = MMCK_pseudo(i, lambda, mu, theta, P_star, P_m );
                for j = 1:1:i+1
                    RealCost(i-opt_K_L(ca,idx)+1,1) = RealCost(i-opt_K_L(ca,idx)+1,1) + Pr_Disch(j)*Cost( j-1, c, g, v_0, P_m, alpha);
                end
                RealCost(i-opt_K_L(ca,idx)+1,1)=RealCost(i-opt_K_L(ca,idx)+1,1)+i*A_d;
                if RealCost(i-opt_K_L(ca,idx)+1,1)<RealCost_min
                    RealCost_min = RealCost(i-opt_K_L(ca,idx)+1,1);
                    opt_K_R(ca,idx) = i;
                else
                end
            end
        else
            opt_K_R(ca,idx)=opt_K_U(ca,idx);
        end
        
        if opt_K_R(ca,idx)<=c
            Save (ca,idx) = Gain_max(3,1)/(P_m*g*0.85);
        else
            Save (ca,idx)= RealCost_min/(P_m*g*0.85);
            Save (ca,idx)=1-Save (ca,idx);
        end
end
cap = ['\caption{Cost saved as coming rate increases, when $\mu=$' num2str(mu) ', $\theta=$' num2str(theta) '}'];
string_parameters=['$\lambda$ goes from ' num2str(lambda_min) ' to ' num2str(lambda_max) ',$\mu=$' num2str(mu) ',$\theta = $' num2str(theta)];
axis_options = 'legend entries={Saving in Scheme 1,Saving in Scheme 2}';
x_label = 'Number of EVs coming per hour';
y_label = 'Saving in \%';
filename = 'Eco_lambda.tex';
fid=fopen(['Data_in_TEX/' filename],'w');
c=clock;
infos=['%Obtained from funtion Eco_theta with parameters of' string_parameters ', run on ' int2str(c(3)) '/' int2str(c(2)) '/' int2str(c(1)) ' at ' int2str(c(4)) ':' int2str(c(5)) ':' int2str(c(6))];
fprintf(fid,'%s\n',infos);
fprintf(fid,'%s\n',['{\footnotesize\begin{tikzpicture}\begin{axis}[legend style={at={(1,1.03)},anchor=south east},width=\figwidth,height=\figheight,cycle list name=\mylist,every axis legend/.append style={nodes={right}},xlabel=' x_label ',ylabel=' y_label ',' axis_options ']']);

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:nb_lambda
    fprintf(fid,'%s','(',num2str(lambda_vec(i)),',',num2str(100*Save(1,i)), ')');
end
fprintf(fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:nb_lambda
    fprintf(fid,'%s','(',num2str(lambda_vec(i)),',',num2str(100*Save(3,i)), ')');
end
fprintf(fid, '%s\n', '};');

fprintf(fid,'%s\n','\addplot coordinates{');
for i = 1:1:nb_lambda
    fprintf(fid,'%s','(',num2str(lambda_vec(i)),',',num2str(10), ')');
end
fprintf(fid, '%s\n', '};');

fprintf(fid,'%s\n','\end{axis}\end{tikzpicture}}');
fprintf(fid,'\n%s\n',cap);
fclose(fid);
end

%             B1 (k,lambda-29)= Pr(1,k+1)/sum(Pr(1,1:1:k+1));
%             B2 (k,lambda-29)= Pr(2,k+1)/sum(Pr(2,1:1:k+1));
%             Pr_test1 = MMKK(k, lambda, mu);
%             Block_test1(k,lambda-29) = Pr_test1(k+1);
%             Pr_test2 = MMCK_pseudo(k, lambda, mu, theta, P_star, P_m );
%             Block_test2(k,lambda-29) = Pr_test2(k+1);
