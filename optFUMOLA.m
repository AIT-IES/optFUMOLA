function res=optFUMOLA(obj, algorithm, options)
%***************************************************************************
% This software is designed to work with FUMOLA
% https://sourceforge.net/projects/fumola/ 
% The conference paper ... should serve as documentation and further
% details are provided throughout the code as comments. 
% Any problems with optFUMOLA? Feel free to contact me
% benedikt.pesendorfer@ait.ac.at 
%
% optFUMOLA aims at providing an optimization library with suitable
% algorithms for simulation-based optimization. 
% Therefore different algorithms are include and interfaced with optFUMOLA
% in order to improve useability. 
% To get started you should provide all simulation details in objFUN<user>
% (subclass of the base class objFUNbase). 
%
% The optimization and simulation process is then started by calling the
% function optFUMOLA(obj, algorithm, options) with: 
%	- obj ... object with information provided in the subclass ObjFUN<user>
%	- see ObjFUN<user>.m and ObjFUNbase.m for more information; obj =
%	ObjFUN<user>(runobj);  
%             with: runobj ... object with information provided in the
%             subclass RunTASKS<user> - see RunTASKS<user>.m and
%             RunTASKSbase.m for more information; runobj =
%             RunTASKS<user>(runobj);   
%	- algorithm ... your choice of optimization algorithm; implemented
%	algorithms: 
% 		- 'DE' ... Differential Evolution
%		- 'PSO' ... Particle Swarm
%		- 'MATSuMoTo' ... Efficient Global Optimization with Surrogate
%		Model 
%		- 'PSwarm' ... Pattern Search with Particle Swarm Search Step
%		- 'NSGA-II'	--- Multi-Objective optimization --> specify more than
%		one objective function in ObjFUN<user>.m 
%	- options ... further option and definition of the optimization
%	problem; this struct can contain the following fields: 
%		- .LB ... lower bounds
%		- .UB ... upper bounds
%		- .A_eq ... left side of equality constraint A_eq*x <= b_eq; leave
%		empty if no constraints apply; (only applicalbe for PSO) 
%		- .b_eq ... right side of equality constraint A_eq*x <= b_eq; leave
%		empty if no constraints apply; (only applicalbe for PSO) 
%		- .A_ineq ... left side of inequality constraint A_ineq*x <=
%		b_ineq; leave empty if no constraints apply 
%		- .b_ineq ... right side of inequality constraint A_ineq*x <=
%		b_ineq; leave empty if no constraints apply 
%		- .maxfunevals ... maximum number of function evaluations
%		- .nvars ... number of design variables
%		- .npop ... number of population members (for DE, PSO, PSwarm and
%		NSGA-II) and number of new points selected at each iteration (for
%		MATSuMoTo)  
%       - .nobj ... number of objective functions (only used for
%       multi-objective optimization with NSGA-II) 
% To access additional options for the different algorithms modify settings
% below!! 
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

% check correct input of struct 'options'
    if ~isfield(options, 'LB')
        error('Lower bound not set: options.LB');
    end
    if ~isfield(options, 'UB')
        error('Upper bound not set: options.UB');
    end
    if ~isfield(options, 'maxfunevals')
        error('Maximum number of function evaluations not set');
    end
    if ~isfield(options, 'nvars')
        error('Number of variables not set: options.nvars');
    end
    if ~isfield(options, 'npop')
        error('Population size not set: options.npop');
    end
    if ~isfield(options, 'A_eq')
        A_eq = [];
    end
    if ~isfield(options, 'b_eq')
        b_eq = [];
    end
    if ~isfield(options, 'A_ineq')
        A_ineq = [];
    end
    if ~isfield(options, 'b_eq')
        b_eq = [];
    end

    switch algorithm
        case 'DE'
            % add folder to path
			addpath('DE')
            
            % algorithm specific parameters
			S_struct.F_VTR = -Inf; %"Value To Reach" (stop when ofunc < F_VTR)
			S_struct.I_bnd_constr = 1;  %1: use bounds as bound constraints, 0: no bound constraints
			S_struct.F_weight = 0.85; %DE-stepsize F_weight ex [0, 2]
			S_struct.F_CR = 1.0; %crossover probabililty constant ex [0, 1]
			S_struct.I_strategy = 2; %1:DE/rand/1, 2:DE/local-to-best/1, 3:DE/best/1 with jitter, 4:DE/rand/1 with per-vector-dither, 5:DE/rand/1 with per-generation-dither, 6:DE/rand/1 either-or-algorithm
			S_struct.I_refresh = 1; %intermediate output will be produced after "I_refresh" iterations. No intermediate output will be produced if I_refresh is < 1
            S_struct.I_plotting = 0;
                        
            % pass parameters
            S_struct.I_funevalmax = options.maxfunevals;
			S_struct.I_D = options.nvars; %number of parameters of the objective function
			S_struct.FVr_minbound = options.LB; %vector of lower bounds of initial population *** note: these are no bound constraints!! ***
			S_struct.FVr_maxbound = options.UB; %vector of upper bounds of initial population *** note: these are no bound constraints!! ***
			S_struct.I_NP = options.npop; %number of population members
			S_struct.I_itermax = round(options.maxfunevals/options.npop); %maximum number of iterations (generations)
            S_struct.A_ineq = options.A_ineq;
			S_struct.b_ineq = options.b_ineq;			
            if isfield(options, 'A_eq') | isfield(options, 'b_eq')
                if options.A_eq ~= [] & options.b_eq ~= []
                    error('DE does not support equality constraints')
                end
            end

            % start the optimization
			[res.x_best,res.f_best] = DE( obj, S_struct );
            
        case 'PSO'
            % add folder to path
			addpath('PSO')
            
            % algorithm specific parameters
            gen = ceil(options.maxfunevals/options.npop);
			psooptions = psooptimset('ParticleInertia',0.3925,'CognitiveAttraction',2.5586,'SocialAttraction',1.3358,'Display','iter','Generations',gen,'PopulationSize',options.npop,'Vectorized','on'); %'Vectorized' 'on' is crucial to work with input/output behavior of optFUMOLA
				% 'PlotFcns',@psoplotbestf) 
                
            % pass parameters
			nvars = options.nvars; %number of variables
			Aineq = options.A_ineq; %non-linear constraints
			bineq = options.b_ineq;
            Aeq = options.A_eq;
            beq = options.b_eq;
			LB = options.LB; %boundary constraints
			UB = options.UB;
                
            % start the optimization
            [res.x_best,res.f_best,res.exitflag,res.output,res.population,res.scores]=pso(obj,nvars,Aineq,bineq,Aeq,beq,LB,UB,[],psooptions)

        case 'MATSuMoTo'  
            % add folder to path
            addpath('MATSuMoTo');
                    
            % algorithm specific parameters
            surogate_model = 'RBFcub'; %selected surrogate model type
            sampling_technique = 'CANDloc'; %global randomized sampling strategy
            initial_design = 'LHS'; %Matlab's lhsdesign.m as initial design
            starting_point = []; %no user-specified points to be added to the initial design
            NumberNewSamples = 1; %1 new point is selected in each iteration
            Problem.parallel_eval = 1; %crucial to work with input/output behavior of optFUMOLA		
            Problem.integer = []; %indices of integer variables
            
                        
            % pass parameters
            Problem.dim = options.nvars; %problem dimesnion
            number_startpoints = options.npop; %number of points in the initial experimental design
            maxeval = options.maxfunevals; %maximum number of allowed function evaluations
            Problem.xlow =  options.LB; %lower bounds
            Problem.xup = options.UB; %upper bounds
            Problem.continuous = (1:options.nvars); %indices of continuous variables
            if isfield(options, 'A_eq') | isfield(options, 'b_eq')
                if options.A_eq ~= [] & options.b_eq ~= []
                    error('MATSuMoTo does not support equality constraints')
                end
            end

            % start the optimization
            [res.x_best, res.f_best]=MATSuMoTo(obj,Problem,maxeval,surogate_model,sampling_technique,initial_design,number_startpoints,starting_point,NumberNewSamples);
            
        case 'PSwarm'
            % add folder to path
			addpath('PSwarm')
            
            % algorithm specific parameters
			pswarmoptions.SearchType=1; %make sure PSO is used in search step
			pswarmoptions.Vectorized=1; %crucial to work with input/output behavior of optFUMOLA
            
            % pass parameters
			Problem.A = options.A_ineq;
			Problem.b = options.b_ineq;
            if isfield(options, 'A_eq') | isfield(options, 'b_eq')
                if options.A_eq ~= [] & options.b_eq ~= []
                    error('PSwarm does not support equality constraints')
                end
            end
			Problem.LB = transpose(options.LB);
			Problem.UB = transpose(options.UB);
			pswarmoptions.Size=options.npop; %Population Size
			pswarmoptions.MaxIter=ceil(options.maxfunevals/options.npop);
			pswarmoptions.MaxObj=options.maxfunevals;

			% start the optimization
			[res.x_best,res.f_best,res.RunData]=PSwarm(obj, Problem, [], pswarmoptions);
            
        case 'NSGA-II'
            % add folder to path
            addpath('NSGA-II');
                        
            % algorithm specific parameters
            if options.nobj < 2
                error('More than one objective function is needed for NSGA-II');
            end
            nsgaoptions = nsgaopt();
            nsgaoptions.useVectorization = 'yes'; %crucial to work with input/output behavior of optFUMOLA
            
            % pass parameters
            if isfield(options, 'A_eq') | isfield(options, 'b_eq')
                if options.A_eq ~= [] & options.b_eq ~= []
                    error('NSGA-II does not support equality constraints')
                end
            end
            nsgaoptions.numCons = size(options.A_ineq,1);
            nsgaoptions.A_ineq=options.A_ineq;
            nsgaoptions.b_ineq=options.b_ineq;
            nsgaoptions.popsize = options.npop;
            nsgaoptions.lb = options.LB;
            nsgaoptions.ub = options.UB;
            nsgaoptions.numVar = options.nvars;
            nsgaoptions.numObj = options.nobj;
            nsgaoptions.vartype = ones(1, nsgaoptions.numVar);
                        
            % start the optimization
            res = nsga2(obj, nsgaoptions);
            
        otherwise
            error('Algorithm name not recognized.\n');
            return
    end
    


end