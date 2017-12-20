function [ r_d_min, r_u_min ] = d_u_min( rho_d, rho_u, gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r_u_min=2-rho_u+gamma*rho_u^(-0.5)*(1-rho_u)^1.5;
r_d_min=1-(rho_d-gamma*sqrt(rho_d-(rho_d^2)));
end

