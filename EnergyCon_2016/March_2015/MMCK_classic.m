function [ Pr ] = MMCK_classic( K, lambda, theta, P, P_m )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
c = floor(P_m/P);

Pr = ones(1,K+1);

% Pr (i+1) is the probability of having i cars discharging
for i=1:1:K;
    if i<c+1
        Pr(1,i+1) = Pr(1,i)*lambda/(i*theta*P);
    else
        Pr(1,i+1) = Pr(1,i)*lambda/(theta*P_m);
    end
end
Pr = Pr/sum(Pr);

end

