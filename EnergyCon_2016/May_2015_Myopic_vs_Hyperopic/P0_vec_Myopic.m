function [ p_0, cost_ins, c ] = P0_vec_Myopic( k, g, v_0, epsilon, p_f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p_0 = zeros(1,k+1);
cost_ins = zeros(1,k+1);
p_0_star = (g-v_0)/(2*epsilon*g);
p_star = p_0_star*(1-epsilon*p_0_star);
c = p_f/p_star;

for index = 1:1:k+1 % index - 1 is the discharging EVs
    if index-1 <= c
        p_0(1,index) = p_0_star;
        cost_ins(1,index) = (index-1)*p_0(1,index)*v_0 + (p_f - (index-1)*p_star)*g;
    else
        p_0(1,index) = (1-sqrt(1-4*epsilon*p_f/(index-1)))/(2*epsilon);
        cost_ins(1,index) = (index-1)*p_0(1,index)*v_0;
    end
end

end

