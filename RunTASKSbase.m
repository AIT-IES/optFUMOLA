%***************************************************************************
% This superclass handle only consists of one abstract method 
%   * run_tasks: receives a list of tasks from ObjFUNbase
% Subclass has to be implemented that containts the specific 
% parallelization handling.
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

classdef (Abstract) RunTASKSbase < handle
    methods (Abstract)
        
        % Actual parallelization of all received tasks. Implement this method in a subclass.
        run_tasks(objrun, tasks)
        
    end % methods (Abstract)
    
end % classdef