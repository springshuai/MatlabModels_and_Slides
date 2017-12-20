function [ uber ] = User_welfare( tht_bar, c_b, p_d, p_a, t_c, t_r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

syms tht

if t_r/p_a<t_c/p_d

    uber=int((tht*p_a-t_r*c_b)*(exp(-tht/tht_bar)./tht_bar),tht,t_r*c_b/p_a,(t_c-t_r)*c_b/(p_d-p_a))...
        +int((tht*p_d-t_c*c_b)*(exp(-tht/tht_bar)./tht_bar),tht,(t_c-t_r)*c_b/(p_d-p_a),inf);
else
    uber=int((tht*p_d-t_c*c_b)*(exp(-tht/tht_bar)./tht_bar),tht,t_c*c_b/p_d,inf);
end    
uber=eval(uber);

end

