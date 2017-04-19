function Solution= Optimization_continuous(Data)
%OPTIMIZATION_CONTINUOUS.m is the optimization phase for continuous 
%unconstrained problems: 
%iteratively select new sample site(s), do the expensive
%function evaluation, and update the chosen surrogate model; iterate until
%stopping criterion met
%The algorithm may restart from scratch if it determined being trapped in a
%local optimum
%
%--------------------------------------------------------------------------
%Copyright (c) 2013 by Juliane Mueller
%
%This file is part of MATSuMoTo.m - the MATLAB Surrogate Model Toolbox
%MATSuMoTo is free software: you can redistribute it and/or modify it under
%the terms of the GNU General Public License as published by the Free 
%Software Foundation, either version 3 of the License, or (at your option) 
%any later version.
%
%MATSuMoTo is distributed in the hope that it will be useful, but WITHOUT 
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
%more details.
%
%You should have received a copy of the GNU General Public License along 
%with MATSuMoTo.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%input:
%Data: structure array containing the optimization problem information
%output: 
%Solution -structure containing all information about the optimization
%after optimization stopped
%--------------------------------------------------------------------------

%initialize structure that collects all problem information from all
%restarts
Solution.S=[]; %sample points
Solution.Y=[]; %function values corresponding to sample points
Solution.fbest=[]; %best function value found in each trial
Solution.xbest=[]; %best point found in each trial
Solution.fevaltime=[]; %time needed for function evaluation
Solution.trial=0; %number of trials (within maxeval allowed evaluations)
Solution.EvalsEachTrial = []; %number of function evaluations in each trial


iteration_counter_FUMOLA = -1; %[PES]included to meet optFUMOLA requirements
funeval_counter = 0; %[PES] new function evaluation counter since old one did not work when more than one new trial run

% do until maximum number of allowed function evaluations has been reached
%while size(Solution.S,1)+(size(Data.starting_point,1)+Data.number_startpoints)< Data.maxeval 
while funeval_counter < Data.maxeval %[PES]
    %generate initial experimental design
    rankS=0;
    while rankS < Data.dim+1 %repeat until rank requirement met
        startPoint_cur=Data.starting_point; %user supplied starting points
        Data.S=StartingDesign(Data);  %create initial experimental design
        if isfield (Data,'starting_point') %user gives one or more starting points to add to the initital design
            %check if any given starting point is already contained in initial design
            ii = 1;
            while ii <= size(startPoint_cur,1)
                if any(sum((repmat(startPoint_cur(ii,:),size(Data.S,1),1)-Data.S).^2,2)==0)
                    startPoint_cur(ii,:)=[];
                else
                    ii=ii+1;
                end
            end
            %add remaining starting points to initial experimental design
            Data.S=[startPoint_cur; Data.S];
        end
        rankS = rank([ones(size(Data.S,1)), Data.S]); %compute matrix rank
    end

    if ~Data.parallel_eval
%%------------serial expensive function evaluations------------------------
        %do expensive function evaluations at points in initial experimental design
        Data.Y=zeros(size(Data.S,1),1); %initialize vector for function values
        for ii = 1:size(Data.S,1)
            fevalst=tic; %record function evaluation time 
            Data.Y(ii,1)=feval(Data.objfunction,Data.S(ii,:)); %expensive function evaluation
            Data.fevaltime(ii,1)=toc(fevalst);
        end 
%%------------end serial expensive function evaluations--------------------
    else
%%-------vectorized version of doing function evalutaions--------------------
        fvals=zeros(size(Data.S,1),1); %initialize vector for function values
        %object=Data.objfunction; %objective function handle
        Samples=Data.S;  %points to be evaluated in parallel
        Tfeval=zeros(size(Data.S,1),1);  %initialize vector for evaluation times
        fevalst=tic; %record function evaluation time 
        
        [i_m, i_n] = size(Samples);
        funeval_counter = funeval_counter + i_m;
        iteration_counter_FUMOLA=iteration_counter_FUMOLA+1;
        obj_f = objective_function(Data.obj, Samples, iteration_counter_FUMOLA); %[PES] evaluate all points at once
        obj_f = penalty_function(obj_f, Samples, Data); %[PES] integrate constraints
        
        for ii = 1:size(Data.S,1)  %parallel for
            fvals(ii,1)=obj_f(ii); %expensive function evaluations
            Tfeval(ii,1)=toc(fevalst);
        end  
        Data.Y=fvals;
        Data.fevaltime=Tfeval;
%%--------end of vectorized version------------------------------------------
    end
    
    [Data.fbest,idx]=min(Data.Y);  %best function value found so far
    Data.xbest=Data.S(idx,:); %best point found so far
    Data.Ymed=Data.Y; %vector where large function values are replaced by median
    MedY=median(Data.Y); %median of function values 
    Data.Ymed(Data.Y>MedY)=MedY; %replace large function values by median, used in calculation of surrogate model parameters
    [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data); % determine parameters of surrogate model

    %initialize parameter values
    xrange = Data.xup-Data.xlow; %range between lower and upper variable bounds
    minxrange = min(xrange);%smallest variable range
    % perturbation parameters
    sigma_stdev(1) = 0.2*minxrange;   
    sigma_stdev(2) = 0.1*minxrange;              
    sigma_stdev(3) = 0.05*minxrange;           
    NCandPoint=500*Data.dim;     %number of candidate points in each group    
    maxshrinkparam = 5;% maximal number of reduction of standard deviation for normal distribution when generating the candidate points
    failtolerance = max(5,Data.dim); %tolerance for number of allowed consecutive failed improvement trials
    succtolerance =3; %tolerance for number of allowed consecutive successful improvement trials
    shrinkctr = 0; % number of times perturbation radius has been decreased
    failctr = 0; % number of consecutive unsuccessful iterations
    succctr=0; % number of consecutive successful iterations
    % permutaion probability
    if Data.dim>5 %if problem dimension is greater than 5, perturb each variable with probability 5/dimension
        P=max(0.1,5/Data.dim);
    else %if problem dimension is <= 5, perturb all variables with probability 1
        P=1;
    end
    Data.tolerance = 0.001*minxrange;  % tolerance for deciding when two points coincide
    RScrit=1.25; %initialize weight for response surface criterion
    localminflag=0; %initialize indicator if local min has been encountered
    
    % iteration loop - until max number of function evaluations has been reached or local min encountered
    %while size(Solution.S,1)+ size(Data.S,1) < Data.maxeval && ~localminflag
    while funeval_counter < Data.maxeval && ~localminflag %[PES]
        %number of new points depends on how many evaluations are still in
        %the budget
        %Data.NumberNewSamples=min(Data.NumberNewSamples,Data.maxeval-(size(Solution.S,1)+ size(Data.S,1)));
        Data.NumberNewSamples=min(Data.NumberNewSamples,Data.maxeval-funeval_counter); %[PES]
        
        %sampling strategy
        if strcmp(Data.sampling_technique(1:4), 'CAND') %randomized strategy   
            weights=zeros(1,Data.NumberNewSamples); %initialize vecor with response surface criterion weights
            % determine weights for response surface criterion for each new sample site    
            for ii = 1:Data.NumberNewSamples
                RScrit=RScrit-.25;
                if RScrit < 0
                    RScrit = 1;
                end
                weights(ii) =RScrit;
            end
            CandPoint=Perturbation(Data, NCandPoint, sigma_stdev,P); %create potential new sample points by perturbation
            %select new sample sites
            xselected=SamplePointSelection(Data,CandPoint, weights,...
                lambda,gamma,mmodel,beta, w_m);
            clear CandPoint 
        else %surface minimum sampling
             xselected = SurfMin(Data,lambda,gamma, mmodel,beta, w_m);
        end
        no_xnew=size(xselected,1); %number of new sample sites   
        % perform expensive function evaluation at the selected point
        Fselected=zeros(no_xnew,1); %initialize vector for function values
        if no_xnew >1 && Data.parallel_eval
%%-------vectorized version of doing function evalutaions-------------------- 
            %object=Data.objfunction; %objective function handle
            Samples=xselected; %points to be evaluated
            Tfeval=zeros(no_xnew,1); %initialize vector for evaluation times
            fevalst=tic; %record function evaluation time 
            
            [i_m, i_n] = size(Samples);
            funeval_counter = funeval_counter + i_m;    
            iteration_counter_FUMOLA=iteration_counter_FUMOLA+1;
            obj_f = objective_function(Data.obj, Samples, iteration_counter_FUMOLA); %[PES] evaluate all points at once
            obj_f = penalty_function(obj_f, Samples, Data); %[PES] integrate constraints
            
            for ii = 1:no_xnew
                Fselected(ii,1)=obj_f(ii); %expensive function evaluation
                Tfeval(ii,1)=toc(fevalst);
            end  
            Data.fevaltime=[Data.fevaltime;Tfeval];
%%--------end of vectorized version------------------------------------------
        else         
%%-------------sequentially evaluating new points--------------------------
            for ii=1:no_xnew
                fevalst=tic; %record time for function evaluation
                
                [i_m, i_n] = size(xselected(ii,:));
                funeval_counter = funeval_counter + i_m;
                iteration_counter_FUMOLA=iteration_counter_FUMOLA+1;
                Fselected(ii,1) = objective_function(Data.obj, xselected(ii,:), iteration_counter_FUMOLA); %[PES] evaluate single point
                Fselected(ii,1) = penalty_function(Fselected(ii,1), xselected(ii,:), Data); %[PES] integrate constraints
                
                Data.fevaltime=[Data.fevaltime;toc(fevalst)];
            end
%%-------------end of sequentially evaluating points-----------------------
        end
        
        [Fselected_min, IDmin] = min(Fselected); %select best new point
        if Fselected_min < Data.fbest %new point is better than best point found so far
            if (Data.fbest - Fselected_min) > (1e-3)*abs(Data.fbest) %improvement must be better than 0.1% to be considered a success
                failctr = 0; %re-initialize fail-counter
                succctr=succctr+1; %update success counter
                shrinkctr = 0; %re-initialize shrink-counter
            else
                failctr = failctr+1; %update fail-counter
                succctr=0; %re-initialize success counter
            end     
            Data.xbest = xselected(IDmin,:); %update best point found so far
            Data.fbest = Fselected_min; %update best function value found so far  
        else %new points are no improvement
            failctr = failctr + 1; %update fail-counter
            succctr=0; %re-initialize success counter
        end
    
        % check if algorithm is in a local minimum
        shrinkflag = 1;      
        if failctr >= failtolerance %check how many consecutive failed improvement trials
            if shrinkctr >= maxshrinkparam %perturbation radius has been reduced as often as allowed
                shrinkflag = 0; %no further perturbation radius reduction allowed
            end
            failctr = 0; %re-initialize fail-counter
            if shrinkflag == 1 %reduce perturbation radius
                shrinkctr = shrinkctr + 1; %update shrink counter
                sigma_stdev = max(sigma_stdev/2, Data.tolerance); %reduce perturbation range
            else %local minimum detected / algorithm stuck
                localminflag = 1; %update indicator that local minimum encountered
                disp('Algorithm is probably in a local minimum! Restarting the algorithm from scratch.');
            end
        end

        if succctr>=succtolerance %check if number of consecutive improvements is large enough
            sigma_stdev=min(2*sigma_stdev, max(Data.xup-Data.xlow));%increase search radius
            succctr=0; %re-initialize success counter
            shrinkctr = 0; %re-initialize shrink-counter
        end

        %update data vectors
        Data.S = [Data.S; xselected]; %sample site matrix
        Data.Y = [Data.Y; Fselected]; %objective function values
        Data.Ymed = Data.Y; %data where large function values set to median, for calculation of surrogate model parameters
        % replace large function values by the median of all available function values
        MedY=median(Data.Y);
        Data.Ymed(Data.Y>MedY)=MedY;

        %check if any stopping criterion met
        %if size(Solution.S,1)+ size(Data.S,1)< Data.maxeval && localminflag == 0
        if funeval_counter < Data.maxeval && localminflag == 0 %[PES]
            %do not stop, update model parameters
            [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data);
            %save intermediate results
            %fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
            fprintf('Total number of function evaluations: %4.0f;\n',funeval_counter) %[PES]
            fprintf('Best feasible function value in this trial: %f\n',Data.fbest)
            if ~isempty(Solution.fbest)
                fprintf('Best feasible function value over all previous trials: %f\n',min(Solution.fbest))
            %else
                %fprintf('Best feasible function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution');
        %elseif localminflag == 1 && size(Solution.S,1) + size(Data.S,1) +...
        %        Data.number_startpoints +size(Data.starting_point,1)>= Data.maxeval 
        elseif localminflag == 1 && funeval_counter +...
                Data.number_startpoints +size(Data.starting_point,1)>= Data.maxeval  %[PES]
            %if local min encountered and not sufficiently many points in
            %budget to restart from scratch, do remaining steps as if no
            %local min detected
            %if size(Solution.S,1) + size(Data.S,1)  >= Data.maxeval
            if funeval_counter >= Data.maxeval %[PES]
                %budget of function evaluations exhausted
                Solution.S=[Solution.S;Data.S]; %update sample site matrix
                Solution.Y=[Solution.Y; Data.Y]; %update function value vector
                Solution.fbest=[Solution.fbest; Data.fbest]; %update best function value found
                Solution.xbest=[Solution.xbest; Data.xbest]; %update best point found
                Solution.fevaltime=[Solution.fevaltime; Data.fevaltime]; %update function evaluation time vector
                Solution.trial = Solution.trial + 1; %update number of restarts from scratch
                Solution.EvalsEachTrial = [Solution.EvalsEachTrial;length(Data.Y)]; %update number of evaluations in each trial
                Data.Ymed=[]; %re-initialize median function value vector
                Data.S=[]; %re-initialize sample site matrix 
                Data.Y=[]; %re-initialize function value vector
                Data.fevaltime=[]; %re-initialize function evaluation time vector
                Data.xbest=[]; %re-initialize best point found so far
                %fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
                fprintf('Total number of function evaluations: %4.0f;\n',funeval_counter) %[PES]
                fprintf('Best feasible function value in this trial: %f\n',Data.fbest)
                Data.fbest=[]; %re-initialize best function value found so far
            else           
                Data.localminflag =0; %re-initialize indicator that local min found 
                [lambda,gamma,mmodel,beta, w_m] = FitSurrogateModel(Data); %update surrogate model parameters
                %fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
                fprintf('Total number of function evaluations: %4.0f;\n',funeval_counter) %[PES]
                fprintf('Best feasible function value in this trial: %f\n',Data.fbest)
            end
            if ~isempty(Solution.fbest)
                fprintf('Best feasible function value over all previous trials: %f\n',min(Solution.fbest))
            %else
                %fprintf('Best feasible function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution'); %save intermediate results
        %elseif localminflag == 1 || size(Solution.S,1) + size(Data.S,1) >= Data.maxeval
        elseif localminflag == 1 || funeval_counter >= Data.maxeval %[PES]
            %local min detected or budget of function evaluations exhausted
            Solution.S=[Solution.S;Data.S]; %update sample site matrix
            Solution.Y=[Solution.Y; Data.Y]; %update function value vector
            Solution.fbest=[Solution.fbest; Data.fbest]; %update best function value found
            Solution.xbest=[Solution.xbest; Data.xbest]; %update best point found
            Solution.fevaltime=[Solution.fevaltime; Data.fevaltime]; %update function evaluation time vector
            Solution.trial = Solution.trial + 1; %update number of restarts from scratch
            Solution.EvalsEachTrial = [Solution.EvalsEachTrial;length(Data.Y)];  %update number of evaluations in each trial
            Data.S=[]; %re-initialize sample site matrix 
            Data.Y=[]; %re-initialize function value vector
            Data.fevaltime=[]; %re-initialize function evaluation time vector            
            Data.xbest=[]; %re-initialize best point found so far
            %fprintf('Total number of function evaluations: %4.0f;\n',size(Solution.S,1) +  size(Data.S,1))
            fprintf('Total number of function evaluations: %4.0f;\n',funeval_counter) %[PES]
            fprintf('Best feasible function value in this trial: %f\n',Data.fbest)
            Data.fbest=[]; %re-initialize best function value found so far
            if ~isempty(Solution.fbest)
                fprintf('Best feasible function value over all previous trials: %f\n',min(Solution.fbest))
            %else
                %fprintf('Best feasible function value over all previous trials: %f\n',nan)
            end
            save('Results','Data','Solution'); %save intermediate results
        end
    end
end %while for restart
end %function

% apply penalty function [PES]
function f = penalty_function(obj_f, pop, Data)
% loop over all population members
for i=1:size(pop,1)
    gx(i) = 0;
    % loop over all constraints
    for j=1:size(Data.b_ineq,1)
        gx(i) = gx(i)+Data.A_ineq(j,:)*pop(i,:)'-Data.b_ineq(j,:); %compute g(x)<=0
    end
    if gx(i) > 0 %constraints not fulfilled
        f(i) = obj_f(i) + Data.mu_ineq*gx(i);
    else %constraints are fulfilled
        f(i) = obj_f(i);
    end
end
end