function [ Rx ] = test_diff( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

syms x rho_u rho_d rho_n guamma r_u r_d theta P_n P_d C_B t

Rx = theta*P_d*exp(-1-(1+((rho_u*r_u*x-rho_d*(1-r_d)*(1-x)-x)/(rho_d+rho_n*x)))*C_B*t/(P_d*theta*(1-(rho_d+rho_n*x-guamma*sqrt(rho_n*x^2+rho_d-(rho_d+rho_n*x)^2)))))...
     + theta*P_d*(rho_d+rho_n*x-guamma*sqrt(rho_n*x^2+rho_d-(rho_d+rho_n*x)^2))*(exp(-1+C_B*t*((rho_u*r_u*x-rho_d*(1-r_d)*(1-x)-x)/(rho_d+rho_n*x))/(P_d*theta*(rho_d+rho_n*x-guamma*sqrt(rho_n*x^2+rho_d-(rho_d+rho_n*x)^2))))...
                 +exp(-1-(1+((rho_u*r_u*x-rho_d*(1-r_d)*(1-x)-x)/(rho_d+rho_n*x)))*C_B*t/(P_d*theta*(1-(rho_d+rho_n*x-guamma*sqrt(rho_n*x^2+rho_d-(rho_d+rho_n*x)^2))))));





end

