function [f, cons] = objective_fun(x)
% Number of objective is two, while it can have arbirtarly many decision
% variables within the range -5 and 5. Common number of variables is 3.

[m,n] = size(x);

for i = 1:m
    % Objective function one
    sum = -10*exp(-0.2*sqrt((x(i,1))^2 + (x(i,2))^2));
    % Decision variables are used to form the objective function.
    f(i,1) = -sum;

    % Objective function two
    sum = (abs(x(i,1))^0.8 + 5*(sin(x(i,1)))^3) + (abs(x(i,2))^0.8 + 5*(sin(x(i,2)))^3);
    % Decision variables are used to form the objective function.
    f(i,2) = -sum;
    if(x(i,1)-x(i,2)>0) 
        disp([num2str(x(i,1)),' - ',num2str(x(i,2)),' > 0']);
        cons(i,1) = x(i,1) - x(i,2);
    else
        cons(i,1) = 0;
    end
end

%cons = []; %pass empty if no constraints