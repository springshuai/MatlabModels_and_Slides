function [ Result ] = Cost( n, c, g, v_0, P_m, epsilon)
%the cost of electricity when having 'n' cars discharging
%   Detailed explanation goes here
pEV = sqrt(v_0*g);
P_star = (1-sqrt(v_0/g))/epsilon;
if n<=c
    Result = g*P_m - n*P_star*(g-pEV);
else
    Result = P_m*v_0/(1-epsilon*(P_m/n));
end

end

