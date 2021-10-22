%y = 0:10;
%y = y/10;
%x = [];
y = 0.5
% Variable coefficients for the quantum 0 and 1 states
P = [(0.5^(1/3)),(0.5^(1/3)),(0.5^(1/3)),0.5];
C = [0.5,0.5,0.5,0.5];

% Defining the states with i,j photons in each mode
state00 = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
state01 = [0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
state02 = [0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0];
state03 = [0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0];
state10 = [0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0];
state11 = [0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0];
state12 = [0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0];
state13 = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0];
state20 = [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0];
state21 = [0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0];
state22 = [0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0];
state23 = [0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0];
state30 = [0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0];
state31 = [0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0];
state32 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0];
state33 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1];

for r = y 

  %  phi00 = (P(4)/sqrt(6))*((r^(3/2))*state30 + 3*r*sqrt(1-r)*state21 + 3*(1-r)*sqrt(r)*state12 + ((1-r)^(3/2))*state03) + (P(3)/sqrt(2))*(r*state20 + 2*sqrt(r*(1-r))*state11 + (1-r)*state02) + P(2)*(sqrt(r)*state10 + sqrt(1-r)*state01) + P(1)*state00
  %  phi11 = (C(4)/sqrt(6))*(((1-r)^(3/2))*state30 + 3*(1-r)*sqrt(r)*state21 + 3*r*sqrt(1-r)*state12 + (r^(3/2))*state03) + (C(3)/sqrt(2))*((1-r)*state20 + 2*sqrt(r*(1-r))*state11 + r*state02) + C(2)*(sqrt(1-r)*state10 + sqrt(r)*state01) + C(1)*state00
    phi00 = (P(3)/sqrt(2))*(r*state20 + 2*sqrt(r*(1-r))*state11 + (1-r)*state02) + P(2)*(sqrt(r)*state10 + sqrt(1-r)*state01) + P(1)*state00;
    phi11 = (C(3)/sqrt(2))*((1-r)*state20 + 2*sqrt(r*(1-r))*state11 + r*state02) + C(2)*(sqrt(1-r)*state10 + sqrt(r)*state01) + C(1)*state00;

    outerphi00 = phi00 * phi00'
    disp(outerphi00*outerphi00)
    outerphi11 = phi11 * phi11';
    operator = (1/2) * (outerphi00 + outerphi11);
    
    % the SDP
    cvx_begin sdp quiet
        variable X(16,16) hermitian semidefinite
        maximize (trace(operator*X))
        trace(X) == 1
        X >= 0
    
    cvx_end
    
    %x(end + 1) = cvx_optval;
    
end

% now for finding the total cheating probability
%P_A = []
%for i = x
%    P_A(end+1) = (1/2)*i + (1/2);
%end

% disp(cvx_optval)

%disp(P_A)
%plot(y,P_A)