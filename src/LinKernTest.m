% LinKern Tests

%% Test 1
%--------------------------------------------------------------------------
% Test Noon-bean transformation on GTSP instance
% using various TSP solvers (LinKern, LKH)
%--------------------------------------------------------------------------

% Global Script Paramters
%--------------------------------------------------------------------------
SOLVER = 'LKH';         % [ LinKern | LKH ]

%SOLVER = 'LKH';

% Load Data file
%--------------------------------------------------------------------------
map = load('guitarnotes.mat'); % x and set
x = map.map.x;
set = map.map.set;

% Sort node list according to sets
%--------------------------------------------------------------------------
numCities = max(set);
xSorted = [];
setSorted = [];

for i = 1:numCities
    currentSet = find(set == i);
 
    for j = 1:length(currentSet)
        setSorted = [setSorted; i];
        xSorted = [xSorted; x(currentSet(j),:)];
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

% Transform ATSP to TSP if needed
%--------------------------------------------------------------------------

if (strcmp(SOLVER,'LinKern'))
    symtsp = atsp_to_tsp(atspMatrix, infcost);
elseif (strcmp(SOLVER,'LKH'))
    symtsp = atspMatrix;
    %symtsp = transformATSP2TSP(atspMatrix, infcost);
else
    disp('Unknown Solver');
end

% Solver TSP using SOLVER
vertexSequence = [];

if (strcmp(SOLVER,'LinKern'))
    vertexSequenceOrdered = get_linkern_result(symtsp,setSorted);
elseif (strcmp(SOLVER,'LKH'))   
     vertexSequenceOrdered = get_lkh_result(symtsp, setSorted);
else
    disp(' Cannot solve using unknown solver ');
end
    


%% Prune ATSP solution to the final GTSP solution

vertexSequenceFinal = [];
currentVertexIdx = 1;
setsFound = [];

%while(setsFound < numCities)
while(currentVertexIdx < length(setSorted))
    
    currentVertex = vertexSequenceOrdered(currentVertexIdx);
    vertexSequenceFinal = [vertexSequenceFinal currentVertex];
    thisSet = setSorted(currentVertex);
    setsFound = [setsFound; thisSet];
    indexes = find(setSorted == thisSet);
    currentVertexIdx = currentVertexIdx + (length(indexes));


end


%% Solver Analysis

% Did any sectors go unvisited
setsNotFound = [];

for i=1:numCities
    
    if (isempty(find(setsFound == i)))
        
        setsNotFound = [setsNotFound i];
    end
end

fprintf(1,'Sets not found : \n');
disp(setsNotFound');

%% Figure

%vertexSequenceFinal = vertexSequenceOrdered;

figure(1);
hold on;

% plotting
for i = 1:length(xSorted)
    
    scatter(xSorted(i,1), xSorted(i,2), 60, [(1-setSorted(i)/ size(setSorted,1)) 0.4 0.7], 'filled' );
    
end

% draw lines to mean
numsets = max(setSorted);

for i = 1:numsets
   vertices = find(setSorted == i);
   
   if(isempty(vertices) ~= 1)
      
       meanPt(1) = mean(xSorted(vertices,1));
       meanPt(2) = mean(xSorted(vertices,2));
       
       scatter(meanPt(1), meanPt(2), 0.1, 'b', 'filled');
       
       % lines to points
       for j = 1:length(vertices)
           
           line([meanPt(1) xSorted(vertices(j),1)], [ meanPt(2), xSorted(vertices(j),2)]);
       end
       
   end
    
    
end


for i = 1:length(vertexSequenceFinal)
    fret = xSorted(vertexSequenceFinal(i),1);
    string = xSorted(vertexSequenceFinal(i),2);
    
    fprintf(1, '(%d, %d) \n', fret, string);
    
end


totalcost = 0;
for i = 1:length(vertexSequenceFinal)-1
    
    X = [xSorted(vertexSequenceFinal(i),1) xSorted(vertexSequenceFinal(i+1),1)];
    Y = [xSorted(vertexSequenceFinal(i),2) xSorted(vertexSequenceFinal(i+1),2)];
    
    
    line(X,Y, 'Linewidth', 2, 'Color', 'r');
    
    totalcost = totalcost + norm(xSorted(vertexSequenceFinal(i+1)) - xSorted(vertexSequenceFinal(i)));
    
    
end

titleSentence = strcat('Hamiltonian Path:  Unique Notes on a 12 Fret Guitar,  Cost = ',num2str(totalcost));

title(titleSentence);

axis equal;
