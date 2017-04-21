%***************************************************************************
% Template for execution of list of tasks.
% This implementation of RunTASKSbase contains the actual realization of 
% the abstract method defined in RunTASKSbase.m
%   * run_tasks: receives a list of tasks from and executes these tasks.
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

%rename your class according to chosen file name
classdef RunTASKStemplate < RunTASKSbase
    
    methods
        
        function status = run_tasks( objrun, tasks )
		% Input:
		%	- objrun ... object of class RunTASKSbase
		%	- tasks ... cell array; list of tasks as defined in setup_single_task in ObjFUN<user>.m - single task is accessible through "tasks{i}"
		% Output:
		%	- status ... non-zero if error occured when executing tasks, zero otherwise
 
 
            % Run all tasks sequentially on local machine.
            for i=1:size(tasks,2)
                status = system( tasks{i} );
            end
			
			% Run all tasks using MATLAB's parallelization toolbox
			parfor i=1:size(tasks,2)
                status = system( tasks{i} );
            end
            
        end
        
    end % methods
    
end % classdef