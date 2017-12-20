function [Pi, TM] = Steady_state_distribution_myopic( K, lambda, mu, theta, p_0_star, P_m )
%Transition_Matrix writes the transition matrix of a queueing system given the parameters 
%
epsilon = 0.01;
p_star = p_0_star*(1-epsilon*p_0_star);
C = floor (P_m/p_star);
MumberOfStates = (K+1)*(K+2)/2;  % This tells how many states are there given the parking lot size "K"
TM = zeros(MumberOfStates,MumberOfStates);  % This is the transition matrix
Pi = zeros (K+1,1);
State = zeros; 
Count = 0;
for m = 0:1:K
    for n = m:1:K
        Count = Count+1;
        State (m+1,n+1)= Count;
    end
end

function [index]=S(m,n)
    % This function maps the state of "m discharging,n parking" to its index in the TM matrix
    index = State(m+1,n+1);
end

for m = 0:1:K
    for n = m:1:K
        if K-n > 0
            TM (S(m,n),S(m+1,n+1))= lambda;  % new cars arrival
        end
        if n-m >= 1
            TM (S(m,n),S(m,n-1))= (n-m)*mu;   % a depleted car departs
        end
        if m > 0
            TM (S(m,n),S(m-1,n-1))= m*mu;   % a discharging car departs
        end
        if m > 0
            if C - m > 0
                TM (S(m,n),S(m-1,n))= theta*m*p_0_star;  % a discharging car depleted when discharging full rate
            else
                TM (S(m,n),S(m-1,n))= theta*m*(1-sqrt(1-4*epsilon*P_m/m))/(2*epsilon);  % a discharging car depleted when discharging partial rate
            end
        end
        TM (S(m,n),S(m,n))= - sum (TM(S(m,n),:));
    end
end

Pi_raw = null(transpose(TM));   % compute the kernal space
Pi_raw = Pi_raw/sum(Pi_raw);    % compute the probability distribution of all the states(m,n)
for m = 0:1:K
    for n = m:1:K
        Pi(m+1) = Pi(m+1) + Pi_raw (S(m,n));    % compute the probability distribution of states(m)
    end
end
Pi = transpose(Pi);
end

