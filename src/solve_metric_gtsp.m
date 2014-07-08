function vertexSequence = solve_metric_gtsp(X, setMap, solver)


% verifying inputs
%--------------------------------------------------------------------------

% solver verification

if (strcmp(solver, 'LKH'))
    fprintf(1, 'Using LKH solver ... \n');
elseif(strcmp(solver,'Linkern'))
    fprintf(1, 'Using Linkern solver ... \n');   
else
    error('Invalid solver name. User either [Linkern] or [LKH]\n');
end

% verifying X and setMap

if(size(X,1) ~= size(setMap,1))
    error('X and setMap must be of the same size \n');   
end


% Sort node list according to sets
%--------------------------------------------------------------------------
numSets = max(setMap);
xSorted = [];
setSorted = [];

for i = 1:numSets
    currentSet = find(set == i);
 
    for j = 1:length(currentSet)
        setSorted = [setSorted; i];
        xSorted = [xSorted; X(currentSet(j),:)];
    end
end


% Create adjacency matrix (euclidean distances);
%--------------------------------------------------------------------------
numElements = length(setSorted);
A = zeros(numElements, numElements);

for i = 1:numElements
    
    for j = 1:numElements
        
        A(i,j) = norm(xSorted(i,:) - xSorted(j,:));
    end
end


% Transform GTSP to ATSP
%--------------------------------------------------------------------------
[atspMatrix infcost] = gtsp_to_atsp(A, setSorted); 


% Transform ATSP to TSP if needed (Needed only for Linkern)
%--------------------------------------------------------------------------
if (strcmp(solver,'LinKern'))
    symtsp = atsp_to_tsp(atspMatrix, infcost);
elseif (strcmp(solver,'LKH'))
    symtsp = atspMatrix;
else
    error('Invalid solver name. User either [Linkern] or [LKH]\n');
end

% Solver TSP using chosen solver
%--------------------------------------------------------------------------

% figure out getLinkernResult

if (strcmp(solver,'LinKern'))
    vertexSequenceOrdered = getLinkernResult(symtsp,setSorted);
elseif (strcmp(solver,'LKH'))   
     vertexSequenceOrdered = getLKHResult(symtsp, setSorted);
else
    disp(' Cannot solve using unknown solver ');
end



end