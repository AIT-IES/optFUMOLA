%***************************************************************************
% This is a template for an implementation of ObjFUNbase 
% The class below implements the following functionality:
%   * setup_single_task: specify individual simulation task given the
%     parameters from the optmization algorithm
%   * result_single_task: retrieve results from an individual simulation
%     task
%
% --------------------------------------------------------------
% Copyright (c) 2017, AIT Austrian Institute of Technology GmbH.
% All rights reserved.
% --------------------------------------------------------------
%***************************************************************************

%rename your class according to chosen file name
classdef ObjFUNtemplate < ObjFUNbase
    
    methods
        
        % Constructor.
		% rename constructor according to chosen file name
        function obj = ObjFUNtemplate(objrun)
		% Input:
		%	- objrun ... object of class RunTASKSbase
		% Output:
		%	- obj ... object of class ObjFUNbase
            
            output_dir_base = '<path>'; % directory base where results are saved in the form: output_dir_base\output_dir\result_file_name
            retry_timeout = 2; % timeout in seconds if file is not found or can not be read
            n_retry_max = 3; % number of times reading file is tried before error
            dry_run = 0; % 0 ... simulations are run; 1 ... simulations are not run and existing result files will be used (see manual!)
            obj@ObjFUNbase( objrun, output_dir_base, n_retry_max, retry_timeout, dry_run );
            
        end
        
        
        % Setup an individual simulation task.
        function sim_task = setup_single_task( obj, output_dir, task_iter, x )
		% Input:
		%	- obj ... object of class ObjFUNbase
		%	- output_dir ... string; directory where results have to be stored according to: output_dir_base\output_dir\result_file_name (also used by result_single_task)
		%	- task_iter ... integer; iteration number of task in this iteration to identify task and enable unique result filename
		%	- x ... array of numbers; design parameter values
		% Output:
		%	- sim_task ... cell array; task that will be invoked by RunTASKSbase.m 
            
            % name of the file were results should be saved by the simulation tool (number of this task for this iteration is passed to gurantee for unique result file name or folder)
			result_file_name = [ 'result', num2str(task_iter), '.csv' ];
		  
            %define command to run simulation
			%passing parameters through command line arguments
            sim_task = [ '<your_tool>.exe ', num2str(x(1)), num2str(x(2)), result_file_name]; % use synthax of your simulation tool and use as many input parameters x(n) as defined
            
			%passing parameters through files
			input_file_name = [ '<path>', 'input.csv']; %must be unique for each task otherwise file is just overwritten by next task
			fileID = fopen(input_file_name,'w'); % open file
			fprintf(fileID,'%g %g %g\n', x(1), x(2), x(3)); % write parameters to file (use synthax used by your simulation tool) and use as many input parameters x(n) as defined
			fclose(fileID); % close file
			sim_task = [ '<your_tool>.exe' ]; % tool must use the file input_file_name
			
			%convert sim_task (string array) to cell array	
			sim_task = cellstr(sim_task);
		           
        end
        
        % Retrieve results from simulation task and define objective function
        function sim_result = result_single_task( obj, output_dir, task_iter)
		% Input:
		%	- obj ... object of class ObjFUNbase
		%	- output_dir ... string; output directory for this algorithm iteration according to: output_dir_base\output_dir\result_file_name (also used by setup_single_task)
		%	- task_iter ... integer; iteration number of task in this iteration to identify task and enable unique result filename
        % Output:
		%	- sim_result ... numeric (integer or float); objective function used by optimization algorithm and constructed out of simulation results
		
            % Set result file name for current task (same as in setup_single_task)
            output_file_name = [ output_dir, '\result', num2str(task_iter), '.csv' ];
            
            % check if result file is ready
            error = check_result_file( obj, output_file_name );
            
            % read data from result file
            try
                data = dlmread( output_file_name ); % can also be another way of reading data with MATLAB
            catch % reading the data did not work as expected, make a short pause and then try to re-read
                pause( obj.retry_timeout );
                data = dlmread( output_file_name );
            end
            
			% definition of result values
            part1 = data( end, 1 ); 
            part2 = data( end, 2 );
            part3 = data( end, 3 );
			
			%define your objective function (for single objective)
            sim_result = 2*part1 + part2*part3;
            
            %define your objective function (for single objective)
            sim_result(1) = 2*part1; %minimize 2*part1
            sim_result(2) = -part2*part3;  %maximize part2*part3
        end
        
    end % methods
    
end % classdef