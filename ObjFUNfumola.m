%***************************************************************************
% For simulation with FUMOLA based on PtolemyII. 
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

classdef ObjFUNfumola < ObjFUNbase
    
    methods
        
        % Constructor.
        function obj = ObjFUNfumola(objrun)
            
            output_dir_base = 'E:\DeMat_HybridEnergySystem\_results';
            retry_timeout = 2;
            n_retry_max = 3;
            dry_run = 0;
            obj@ObjFUNbase( objrun, output_dir_base, n_retry_max, retry_timeout, dry_run );
            
        end
        
        
        % Setup an individual simulation task.
        function sim_task = setup_single_task( obj, output_dir, task_iter, x )
            
            % task_iter ... number of this population member (passed to gurantee for unique result file name or folder)
            
            % FUMOLA model name (PtolemyII xml-file)
            model_file = 'hybrid_energy_system.xml';
            % This directory (or its sub-directories) is supposed to contain all relevant
            % resources (co-simulation model, FMUs, additional data) for a single task
            % run with FUMOLA.
            resources_dir = 'E:\DeMat_HybridEnergySystem\HybridEnergySystem\fumola';
            % This directory is supposed to be used as working directory for a task and is
            % normaly kept empty.
            working_dir = 'X:\';
            % FUMOLA batch file name
            batch_file = 'E:\DeMat_HybridEnergySystem\runFumola.bat';
            
            % Set result file name for current task.
            task_result_file_name = [ 'result', num2str(task_iter), '.csv' ];
            
            el_heatpump_load_pon = x(1);
            el_store_c = x(2);
            
            % convert parameters to strings
            str_el_heatpump_load_pon = num2str( el_heatpump_load_pon );
            str_el_store_pmax_discharge = num2str( el_heatpump_load_pon );
            
            str_el_store_c = num2str( el_store_c );
            str_el_store_pmax_charge = num2str( 0.2 * el_store_c );
            
            str_pv_prod_scale = '1.5';
            
            %define command to run simulation
            task_optional_fumola_args = [ '-ResultsWriter.filename \"', task_result_file_name, '\" -sys.startValues \"thStoreM = 5000, elStoreC = ', str_el_store_c, ', elStoreP = 0, elHeatpumpLoadP = 0, pvProdScale = ', str_pv_prod_scale, '\" -ctrl.startValues \"thBoilerProdPon = 40000, elStorePmaxCharge = ', str_el_store_pmax_charge,', elStorePmaxDischarge = ', str_el_store_pmax_discharge,', elHeatpumpLoadPon = ', str_el_heatpump_load_pon, '\"' ];
            sim_task = [ 'python.exe E:\DeMat_HybridEnergySystem\runFumola.py ', working_dir, ' ', output_dir, ' ', resources_dir, ' ', batch_file, ' ', model_file, ' ', task_optional_fumola_args ];
            sim_task = cellstr(sim_task);
            
        end
        
        % Retrieve results from simulation task and define objective function
        function sim_result = result_single_task( obj, output_dir, task_iter)
            
            % Set result file name for current task (same as in setup_single_task)
            output_file_name = [ output_dir, '\result', num2str(task_iter), '.csv' ];
            
            % check if result file is ready
            error = check_result_file( obj, output_file_name );
            
            % read data from result file
            try
                data = dlmread( output_file_name );
            catch % reading the data did not work as expected, make a short pause and then try to re-read
                pause( obj.retry_timeout );
                data = dlmread( output_file_name );
            end
            
            th_total_e = data( end, 1 );
            th_boiler_total_e = data( end, 2 );
            el_store_utilization = data( end, 6 );
            sim_result = th_boiler_total_e / th_total_e / el_store_utilization;
            
            % For multiobjective optimization:
            % sim_result(1) = th_boiler_total_e / th_total_e;
            % sim_result(2) = -el_store_utilization;  %maximize
        end
        
    end % methods
    
end % classdef