function [pop, state] = evaluate(opt, pop, state, optFUMOLAobj, varargin)
% Function: [pop, state] = evaluate(opt, pop, state, varargin)
% Description: Evaluate the objective functions of each individual in the
%   population.
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

N = length(pop);
allTime = zeros(N, 1);  % allTime : use to calculate average evaluation times

%*************************************************************************
% Evaluate objective function in parallel
%*************************************************************************
if( strcmpi(opt.useParallel, 'yes') == 1 )
    curPoolsize = matlabpool('size');
    
    % There isn't opened worker process
    if(curPoolsize == 0)
        if(opt.poolsize == 0)
            matlabpool open local
        else
            matlabpool(opt.poolsize)
        end
        % Close and recreate worker process
    else
        if(opt.poolsize ~= curPoolsize)
            matlabpool close
            matlabpool(opt.poolsize)
        end
    end
    
    parfor i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d , Individual: %d / %d \n', state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end
    
    %*************************************************************************
    % Evaluate vectorized objective function [PES]
    %*************************************************************************
elseif( strcmpi(opt.useVectorization, 'yes') == 1 )
    %fprintf('\nEvaluating the objective function... Generation: %d / %d , Individual: %d / %d \n', state.currentGen, opt.maxGen, i, N);
    
    tStart = tic;
    
    for i = 1:N %create array with all population members for data handling purpose
        pop_temp(i,:) = pop(i).var;
    end
    obj_temp = objective_function(optFUMOLAobj, pop_temp, state.currentGen-1); %[PES]
    if( ~isempty(pop(i).cons) ) %check for and calculate constraints
        cons_temp = constraints(pop_temp, opt);
    end
    
    allTime = toc(tStart);
    for i = 1:N
        pop(i).obj(:) = obj_temp(i,:);  %pass objective to population members
        if( ~isempty(pop(i).cons) ) %pass constraints to population members
            pop(i).var(:)
            pop(i).cons(:) = cons_temp(i,:); %pass constraints to population members
            idx = find( cons_temp(i) );
            if( ~isempty(idx) )
                pop(i).nViol = length(idx);
                pop(i).violSum = sum( abs(cons_temp(i)) );
            else
                pop(i).nViol = 0;
                pop(i).violSum = 0;
            end
        end
    end
    
    
    %*************************************************************************
    % Evaluate objective function in serial
    %*************************************************************************
else
    for i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d , Individual: %d / %d \n', state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end
end

%*************************************************************************
% Statistics
%*************************************************************************
state.avgEvalTime   = sum(allTime) / length(allTime);
state.evaluateCount = state.evaluateCount + length(pop);

end


function [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Function: [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Description: Evaluate one objective function.
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-25
%*************************************************************************

tStart = tic;
y = objfun( indi.var, varargin{:} );
[y, cons] = objfun( indi.var, varargin{:} );
evalTime = toc(tStart);

% Save the objective values and constraint violations
indi.obj = y;
if( ~isempty(indi.cons) )
    idx = find( cons );
    if( ~isempty(idx) )
        indi.nViol = length(idx);
        indi.violSum = sum( abs(cons) );
    else
        indi.nViol = 0;
        indi.violSum = 0;
    end
end
end

% calculate constraint values for inequality constraints [PES]
function gx = constraints(pop, opt)
for i=1:size(pop,1) %loop over whole population
    for j=1:opt.numCons %loop over all constraints
        gx(i,j) = opt.A_ineq(j,:)*pop(i,:)'-opt.b_ineq(j,:);  %compute g(x)
        if gx(i,j) <= 0
            gx(i,j) = 0;
        end
    end
end
end