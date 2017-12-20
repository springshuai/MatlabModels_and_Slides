function [ Result ] = Cost( n, c, g, v_0, P_m, alpha)
%the cost of electricity when having 'n' cars discharging
%   Detailed explanation goes here

if n<=c
    Result = g*P_m-n*(g-v_0)^2/(4*alpha);
else
    Result = alpha*P_m^2/n + v_0*P_m;
end

end

