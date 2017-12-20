function [ Pr ] = MMCK_pseudo_Hyperopic( K, lambda, mu, theta, p_0_star, P_m )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
epsilon = 0.01;
p_star = p_0_star*(1-epsilon*p_0_star);
c = floor (P_m/p_star);

Pr = ones(1,K+1);

% Pr (i+1) is the probability of having i cars discharging
for i=1:1:K;
    if i<c+1
        Pr(1,i+1) = Pr(1,i)*lambda/(i*mu+i*theta*p_0_star);
    else
        Pr(1,i+1) = Pr(1,i)*lambda/(i*mu+i*theta*(1-sqrt(1-4*epsilon*P_m/i))/(2*epsilon));
    end
end
Pr = Pr/sum(Pr);

end

