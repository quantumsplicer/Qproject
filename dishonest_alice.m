y = 0:100;
y = y/100;
x = [];
x_exp = [];
for c = y 
    phi00 = [sqrt(c);sqrt(1-c)];
    phi11 = [sqrt(1-c);sqrt(c)];
    outerphi00 = phi00 * phi00';
    outerphi11 = phi11 * phi11';
    operator = (1/2) * (outerphi00 + outerphi11);

    cvx_begin sdp quiet
        variable X(2,2) hermitian semidefinite
        maximize (trace(operator*X))
        trace(X) == 1
        X >= 0

    cvx_end
    
    x(end + 1) = cvx_optval;
    x_exp(end + 1) = (3/4) + (1/2) * sqrt(c * (1-c));
    
end

disp(x_exp)
plot(y,x_exp)
hold on
plot(y,x)
hold off