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

function [atspAdjMatrix]  = noon_bean_gtsp_to_atsp(gtspAdjMatrix, setMap)


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

% using max instead of sum. slight difference from original Noon-Bean

%totalcost = sum(sum(adjMatrixCostCalcs));
maxCost = max(max(adjMatrixCostCalcs));
beta = 2*maxCost;



for i=1:length(setMap)
    
    for j=1:length(setMap)
        
        if(atspAdjMatrix(i,j) ~= 0 && atspAdjMatrix(i,j) ~= -1)
            
            %------------ NOTE -------------------------------------
            % factor of a 1000 added. Ok since max is used earlier
            % this value can be changed around
            % potential bug: values become too large for TSP solvers. 
            %-------------------------------------------------------
            
            atspAdjMatrix(i,j) = (atspAdjMatrix(i,j) + beta) * 1000; 
            adjMatrixCostCalcs(i,j) = adjMatrixCostCalcs(i,j) + beta;
        end
    end
end

totalcost = sum(sum(adjMatrixCostCalcs));

beta2 = 5*totalcost;


for i=1:length(setMap)
    
    for j=1:length(setMap)
        
        if(atspAdjMatrix(i,j) == -1)
            
            atspAdjMatrix(i,j) = beta2;
        end
    end
end


end