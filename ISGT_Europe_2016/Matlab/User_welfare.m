function [ uber ] = User_welfare( tht_bar, c_b, p_d, p_a, t_s, t_r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

syms tht

if (t_r/p_a)<(t_s/p_d)

    uber=int((tht*p_a-t_r*c_b)*(1/tht_bar)*exp(-tht/tht_bar),tht,t_r*c_b/p_a,(t_s-t_r)*c_b/(p_d-p_a))...
        +int((tht*p_d-t_s*c_b)*(1/tht_bar)*exp(-tht/tht_bar),tht,(t_s-t_r)*c_b/(p_d-p_a),inf);
else
    uber=int((tht*p_d-t_s*c_b)*(1/tht_bar)*exp(-tht/tht_bar),tht,t_s*c_b/p_d,inf);
end    
uber=eval(uber);

end

