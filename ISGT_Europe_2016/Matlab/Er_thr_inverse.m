function [c, thr1,thr2,ts_thr1,ts_thr2,ts_theo1,ts_theo2,ts_1,ts_2,y_min,y_max,y_thr1,y_thr2] = Er_thr_inverse( p_a )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

t = 0.03;
c_b = 50;
theta = 0.3;
p_d = 20;

ts_1 = t+(p_d-p_a)*theta/c_b;
ts_2 = t+p_d*theta/c_b;

ts_vec = linspace(ts_1,ts_2,100);

    for i=1:1:100
        ts = ts_vec(i);

        const = c_b*ts/(theta*(p_d-p_a));
        Numerator = theta*(p_d-p_a)*(1-exp(-const));
        Denominator = c_b*(p_d/p_a-1+exp(const));

        thr1(i) = Numerator/Denominator;

        thr2(i) = thr1(i)*(1+(p_d/p_a-1)*exp(const));

    end
    
[mthr1,ts_id1] = max(thr1);
[mthr2,ts_id2] = max(thr2);

ts_thr1 = ts_vec(ts_id1);
ts_thr2 = ts_vec(ts_id2);
    
y_max = exp(c_b*ts_2/(theta*(p_d-p_a)));
y_min = exp(c_b*ts_1/(theta*(p_d-p_a)));

y_thr1 = 1+sqrt(p_d/p_a);
ts_theo1 = (log(y_thr1))*theta*(p_d-p_a)/c_b;
c = p_d/p_a-1;
y_thr2 = (-1-(c+1)*sqrt(1-c))/(c^2+c-1);
ts_theo2 = (log(y_thr2))*theta*(p_d-p_a)/c_b;



end
