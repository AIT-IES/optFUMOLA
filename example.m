%***************************************************************************
% This file is an example on how to run optFUMOLA once ObjFUNbase and 
% RunTASKSbase are implemented by the user.
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

options.LB = [1000 1000];		% lower bound
options.UB = [200000 200000];	% upper bound
options.A_eq = [];
options.b_eq = [];
options.A_ineq = [1 -1];		% left side of inequality constraint A_ineq*x <= b_ineq [1 -1; 2 -2; 3 -3]
options.b_ineq = [0];			% right side of inequality constraint A_ineq*x <= b_ineq [0; 0; 0]
options.maxfunevals = 300;		% maximum number of function evaluations (number of generations = ceil(maxfunevals/npop)
options.nvars = 2;				% number of variables
options.npop = 20;				% number of population members (for DE, PSO, PSwarm and NSGA-II) and number of new points selected at each iteration (for MATSuMoTo) 
options.nobj = 1;               % number of objective functions (only for NSGA-II)
algorithm = 'DE';				% algorithm used for optimization
runobj = RunTASKSsequential();   % create runTASKS object specifying task execution details (e.g., parallelization)
obj = ObjFUNfumola(runobj); 	% create objFUN object specifying simulation details (task descriptions, etc.) and objective function

result=optFUMOLA(obj, algorithm, options); % start optimization and simulation process

