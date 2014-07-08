%##########################################################################
%
% Transformation from ATSP to TSP :
%
% Inputs: aTSP: Asymmetric TSP adjacency matrix
%         infCost: large value assigned to infinite cost edges
% 
% Outputs: symTSP: symmetric TSP instance
%
%##########################################################################

function symTSP = atsp_to_tsp(aATSP, infcost)

n = size(aATSP,1);

asymTSP = aATSP;

symTSP = zeros(n*2, n*2);


for i=1:n
    for j=1:n
        if(i == j)
            asymTSP(i,j) = -100000;
        end
    end
end
       

symTSP(n+1:end,1:n) = asymTSP;
symTSP(1:n, n+1:end) = asymTSP';

symTSP(1:n,1:n) = infcost;
symTSP(n+1:end,n+1:end) = infcost;




end