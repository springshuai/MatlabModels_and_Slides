function [ Pr ] = MMKK( K, lambda, mu)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Pr = ones(1,K+1);
rho = lambda/mu;
% Pr1(i+1) is the probability of having i cars discharging
for i=1:1:K
    Pr(i+1)= Pr(i)*rho/i;
end
Pr = Pr/sum(Pr);
end