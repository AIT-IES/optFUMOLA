%***************************************************************************
% Run tasks sequentially.
% This implementation of RunTASKSbase contains the actual realization of 
% the abstract method defined in RunTASKSbase.m
%   * run_tasks: receives a list of tasks from and executes these tasks 
%     sequentially.
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

classdef RunTASKSsequential < RunTASKSbase
    
    methods
        
        function status = run_tasks( objrun, tasks )
 
            % Run all tasks sequentially.
            for i=1:size(tasks,2)
                status = system( tasks{i} );
            end
            
        end
        
    end % methods
    
end % classdef