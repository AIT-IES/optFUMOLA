%%%%%%%%%%%%%%%
%%% NSGA-II %%%
%%%%%%%%%%%%%%%
%{
rng('default');
addpath('NSGA-II')

options = nsgaopt();
options.popsize = 20;
options.lb = [1000 1000];
options.ub = [200000 200000];
%options.useParallel
options.useVectorization = 'yes';
options.numVar = 2;
options.numObj = 2;
options.vartype = [1 1];

options.numCons = 1;            
options.A_ineq = [1 -1];    	% left side of inequality constraint A_ineq*x <= b_ineq
options.b_ineq = [0];			% right side of inequality constraint A_ineq*x <= b_ineq

options.maxGen = 100;
%options.objfun = @objective_fun
%options.refPoints = [3.037 -10.91];
%options.refWeight = [1 1];
%options.refEpsilon = 0.001;
%options.refUseNormDistance = 'no'


obj = optFUMOLAusermo(); % create optFUMOLA object specifying FUMOLA simulation details (file paths, etc.)
%objective_function( obj, alg_pop, alg_iter )

res = nsga2(obj, options);


%}







%%%%%%%%%%
%%% DE %%%
%%%%%%%%%%
%
addpath('DE')

S_struct.F_VTR = 0; %"Value To Reach" (stop when ofunc < F_VTR)
S_struct.I_D = 2; %number of parameters of the objective function
S_struct.FVr_minbound = [ 1000 1000 ]; %vector of lower bounds of initial population *** note: these are no bound constraints!! ***
S_struct.FVr_maxbound = [ 200000 200000 ]; %vector of upper bounds of initial population *** note: these are no bound constraints!! ***
S_struct.I_bnd_constr = 1;  %1: use bounds as bound constraints, 0: no bound constraints
S_struct.I_NP = 10; %number of population members
S_struct.I_itermax = 15; %maximum number of iterations (generations)
S_struct.F_weight = 0.85; %DE-stepsize F_weight ex [0, 2]
S_struct.F_CR = 1.0; %crossover probabililty constant ex [0, 1]
S_struct.I_strategy = 2; %1:DE/rand/1, 2:DE/local-to-best/1, 3:DE/best/1 with jitter, 4:DE/rand/1 with per-vector-dither, 5:DE/rand/1 with per-generation-dither, 6:DE/rand/1 either-or-algorithm
S_struct.I_refresh = 1; %intermediate output will be produced after "I_refresh" iterations. No intermediate output will be produced if I_refresh is < 1
S_struct.I_plotting = 0; % I_plotting    Will use plotting if set to 1. Will skip plotting otherwise.
S_struct.I_funevalmax = 50;
S_struct.A = [1 -1];
S_struct.b = [0];

%objrun = RunTASKSvmcontrol(); 
objrun = RunTASKSsequential();
obj = ObjFUNfumola(objrun); % create optFUMOLA object specifying FUMOLA simulation details (file paths, etc.)

[FVr_x,S_y,I_nf] = DE( obj, S_struct );
%


%%%%%%%%%%%
%%% PSO %%%
%%%%%%%%%%%
%{
addpath('PSO')
% Problem definition
nvars = 2; %number of variables
A = [1 -1]; %non-linear constraints
b = [0];
LB = [1000 1000]; %boundary constraints
UB = [200000 200000];
options = psooptimset('ParticleInertia',0.3925,'CognitiveAttraction',2.5586,...
    'SocialAttraction',1.3358,'Display','iter','Generations',100,'PopulationSize',20,...
    'Vectorized','on');
% 'PlotFcns',@psoplotbestf,
%'VelocityLimit',abs(UB-LB) ???
obj = optFUMOLAuser();
[x, fval,exitflag,output,population,scores]=pso(obj,nvars,A,b,[],[],LB,UB,[],options);
%}



%%%%%%%%%%%%%%
%%% PSwarm %%%
%%%%%%%%%%%%%%
%{
addpath('PSwarm')

% Problem definition
Problem.A = [1 -1];
Problem.b = [0];
Problem.LB = [1000;1000];
Problem.UB = [200000;200000];

% Initial guess
%InitPop(1).x=[3; 0.5];
%InitPop(2).x=[1; 0.5];

% Algorithm options
Options.Size=20; %Population Size
Options.MaxIter=2000;
Options.MaxObj=300;
Options.SearchType=1; %make sure PSO is used in search step
%Options.Cache=1; %saves objective function
%Options.LoadCache=1;
%Options.SaveCache=1;
Options.Vectorized=1; %vectorized objective function

obj = optFUMOLAuser();

% Run the algorithm
[x,fx,RunData]=PSwarm(obj, Problem, [], Options);
nfo=RunData.ObjFunCounter
deg=RunData.Degenerate
x
fx

%}

%%%%%%%%%%%%%%%%%
%%% MATSuMoTo %%%
%%%%%%%%%%%%%%%%%
%{
addpath('MATSuMoTo');

rng('default');

maxeval = 200; %maximum number of allowed function evaluations
surogate_model = 'RBFcub'; %selected surrogate model type
sampling_technique = 'CANDloc'; %global randomized sampling strategy
initial_design = 'LHS'; %Matlab's lhsdesign.m as initial design
number_startpoints = 21; %15 points in the initial experimental design
starting_point = []; %no user-specified points to be added to the initial design
NumberNewSamples = 1; %1 new point is selected in each iteration

Problem.xlow = [1000 1000]; % variable lower bounds
Problem.xup = [200000 200000];     % variable upper bounds
Problem.dim = 2; %problem dimesnion
Problem.integer = []; %indices of integer variables
Problem.continuous = (1:2); %indices of continuous variables
Problem.parallel_eval = 1; %evaluate objective function in vectorized mode

Problem.A_ineq = [ 1 -1 ];
Problem.b_ineq = [ 0 ];
Problem.mu_ineq = 0.0001; % f(x)+mu*(A_ineq*x-b_ineq)

[x_best, f_best]=MATSuMoTo(Problem,maxeval,surogate_model,sampling_technique,...
    initial_design,number_startpoints,starting_point,NumberNewSamples);

%}
