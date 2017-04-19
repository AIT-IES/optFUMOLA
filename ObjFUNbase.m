%***************************************************************************
% This superclass handle consists of four methods:
%   * objective_function: specifies the workflow when called by the
%     optimization algorithm
%   * setup_all_tasks: all simulation tasks are defined here given the 
%     parameters from the optmization algorithm and using the abstract
%     method setup_single_task
%   * results_all_tasks: retrieves all results from all individual 
%     simulation tasks using the abstract method result_single_task
%   * check_result_file: checks if result file does exist and is ready to
%     be read, pauses for a short period of time if it does not exist;
%     should be used in the implementation of ObjFUNbase to check the 
%     result file exists before reading from it - to avoid errors
% and two abstract methods:
%   * setup_single_task: defines how a single simulation task has to look 
%     like and has to be defined by the user; see ObjFUNfumola.m as an 
%     example
%   * result_single_task: defines how results from a single simulation run
%     have to be read and how the objective function is computed out of 
%     these results; has to be defined by the user; see ObjFUNfumola.m as 
%     an example
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

classdef (Abstract) ObjFUNbase < handle
    
    properties
        
        objrun          % Class object from subclass of RunTASKSbase
        output_dir_base % All task outputs should be written in sub-directories (generated and named according
                        % to the respective iteration number of the algorithm) of this directory.
        n_retry_max     % Maximum number of times to check if a file exsists.
        retry_timeout   % Timeout (in seconds) before again trying to access a file.
        dry_run         % Set to 1 if result files already exist from previous run otherwise set to 0
        
    end % properties
    
    methods
        
        % Full constructor.
        function obj = ObjFUNbase( objrun, output_dir_base, n_retry_max, retry_timeout, dry_run );
            obj.objrun = objrun;
            obj.output_dir_base = output_dir_base;
            obj.n_retry_max = n_retry_max;
            obj.retry_timeout = retry_timeout;
            obj.dry_run = dry_run;
        end
        
        % Objective function as called by the optimization algorithms
        function alg_f = objective_function( obj, alg_pop, alg_iter )

            if alg_iter == 0 & obj.dry_run ~= 1
                % Clean-up previous result files at first iteration
                % IF ALGORITHM DOES NOT START AT ZERO MANIPULATE ITERATION COUNTER WHEN PASSING TO FUNCTION
                if exist(obj.output_dir_base, 'dir')  
                    try
                    	rmdir(obj.output_dir_base, 's');
                    catch                       
                    	pause(1);
                        rmdir(obj.output_dir_base, 's');
                    end
                end
                mkdir(obj.output_dir_base);
            end
            
            % Set output directory for current iteration.
            output_dir = [ obj.output_dir_base, '\iteration(', num2str( alg_iter ), ')' ];
            mkdir(output_dir);
            
            % Setup the co-simulation tasks.
            tasks = setup_all_tasks( obj, output_dir, alg_pop );
            
            % run all tasks
            if obj.dry_run ~= 1
                run_tasks(obj.objrun, tasks);
            end
            
            % Retrieve and combine results from all simulation runs
            alg_f = results_all_tasks(obj, output_dir, alg_pop);
            
        end
        
        % Setup all simulation tasks
        function tasks = setup_all_tasks( obj, output_dir, alg_pop )
            % Define simulation tasks for every population member.
            for task_iter = 1:size(alg_pop,1)
                tasks(task_iter) = setup_single_task( obj, output_dir, task_iter, alg_pop(task_iter,:) ); % Retrieve task description.
            end
        end
        
        % Retrieve all simulation results
        function alg_f = results_all_tasks(obj, output_dir, alg_pop)
            % Retrieve and combine results from all simulation runs
            for task_iter = 1:size(alg_pop,1)
                alg_f(task_iter,:) = result_single_task( obj, output_dir, task_iter);
            end
        end
        
        % Checks if file exists. Retry this check in case the file has not
        % been found. 
        function error = check_result_file( obj, file_name )
            error = 0;
            i_retry = 0;
            while exist( file_name, 'file' ) == 0 % Check if file exists.
                %disp( [ 'wait for "', file_name, '"' ] );
                i_retry = i_retry + 1; % Increment counter.
                if i_retry > obj.n_retry_max % Maximum numbers of retries.
                    disp( [ '[objFUNbase] File not found: "', file_name, '"' ] );
                    error = -1;
                    return
                end
                % Wait before retrying.
                pause( obj.retry_timeout );
            end
        end
        
    end % methods
    
    
    methods (Abstract)
        
        % Setup an individual simulation task. Implement this function in a subclass.
        sim_task = setup_single_task( obj, output_dir, task_iter, x  )
        
        % Retrieve results from simulation task. Implement this function in a sublass.
        sim_result = result_single_task( obj, output_dir, task_iter  )
        
    end % methods (Abstract)
    
end % classdef