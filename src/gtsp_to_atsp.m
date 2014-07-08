%##########################################################################
%
% Noon-Bean Transformation :
% Convert Generlized TSP to Asymmetric TSP
%
% Inputs: gtspAdjMatrix: Full GTSP Adjacency Matrix
%         setMap: node groupings in the same order as in gtspAdjMatrix
% 
% Outputs: atspAdjMatrix: Full transformed Adjacency Matrix for ATSP
%          infCost: Large numeric value assigned for infinite edge cost
%                   (This is used internally in the toolbox)
%
%##########################################################################

function [atspAdjMatrix infcost]  = gtsp_to_atsp(gtspAdjMatrix, setMap)


atspAdjMatrix = (gtspAdjMatrix*0) - 1; % marking edges
adjMatrixCostCalcs = gtspAdjMatrix*0;

% intra-cluster arcs
numSets = max(setMap);

for i=1:numSets
   
    % intra-cluster arcs
    
    indexes = find(setMap == i);
    
    if(length(indexes) > 1)
        
        for j = 1:length(indexes)-1
            atspAdjMatrix(indexes(j),indexes(j+1)) = 0;
        end

        % close the cluster loop
        atspAdjMatrix(indexes(end),indexes(1)) = 0;

    end
        %inter-cluster arc switching

        for j = 2:length(indexes)

            for k=1:length(setMap)

                if(setMap(indexes(j)) ~= setMap(k))

                    if(gtspAdjMatrix(indexes(j),k) ~= 0)

                        atspAdjMatrix(indexes(j-1),k) = gtspAdjMatrix(indexes(j),k);
                        adjMatrixCostCalcs(indexes(j-1),k) = gtspAdjMatrix(indexes(j),k);
                    end

                end

            end

        end


        for k=1:length(setMap)

            if(setMap(indexes(1)) ~= setMap(k))

                if(gtspAdjMatrix(indexes(1),k) ~= 0)

                    atspAdjMatrix(indexes(end),k) = gtspAdjMatrix(indexes(1),k);
                
                    adjMatrixCostCalcs(indexes(end),k) = gtspAdjMatrix(indexes(1),k);

                end

            end

        end
    
end



%totalcost = sum(sum(adjMatrixCostCalcs));
maxCost = max(max(adjMatrixCostCalcs));
beta = 2*maxCost;


% should i add the edges from the starting vertex out here?


for i=1:length(setMap)
    
    for j=1:length(setMap)
        
        if(atspAdjMatrix(i,j) ~= 0 && atspAdjMatrix(i,j) ~= -1)
            
            atspAdjMatrix(i,j) = (atspAdjMatrix(i,j) + beta) *1000; % random factor of 1000 added beware
            adjMatrixCostCalcs(i,j) = adjMatrixCostCalcs(i,j) + beta;
        end
    end
end

totalcost = sum(sum(adjMatrixCostCalcs));

beta2 = 2*totalcost;


for i=1:length(setMap)
    
    for j=1:length(setMap)
        
        if(atspAdjMatrix(i,j) == -1)
            
            atspAdjMatrix(i,j) = beta2;
        end
    end
end

infcost = beta2;

end