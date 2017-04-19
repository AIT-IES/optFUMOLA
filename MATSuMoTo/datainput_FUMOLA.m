function Data=datainput_myproblem

obj = optFUMOLAuser();

Data.xlow = [1000 1000]; % variable lower bounds
Data.xup = [200000 200000];     % variable upper bounds
Data.dim = 2; %problem dimesnion
Data.integer = []; %indices of integer variables
Data.continuous = (1:2); %indices of continuous variables
%Data.objfunction=@(x)myfun(x); %handle to objective function
Data.parallel_eval = 1; %evaluate objective function in vectorized mode

Data.A_ineq = [ 1 -1 ];
Data.b_ineq = [ 0 ];
Data.mu_ineq = 0.0001; % f(x)+mu*(A_ineq*x-b_ineq)

end %function

%function y=myfun(x) %objective function
%
%y = objective_function(x); %objective function value
%
%end %myfun