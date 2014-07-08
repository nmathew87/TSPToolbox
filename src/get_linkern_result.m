function vertexSequenceOrdered = get_linkern_result(symtsp,setSorted)



%creating the input file for linkern
n = size(symtsp,1);

fid=fopen('superchargerTSP.txt', 'wt');

fprintf(fid,'TYPE: TSP\nDIMENSION:%7.0f\n',n);

fprintf(fid,'EDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT:FULL_MATRIX\nEDGE_WEIGHT_SECTION\n');


for i=1:n
    fprintf(fid, '%7.0f ', round(symtsp(i,:)*1000) );
    fprintf(fid, '\n');
end

fprintf(fid, 'EOF');

fclose(fid);

tic;
% run linkern solver
%!linkern -K 2 -R 2000 -r 10 -I 1 -t 5 -o superChargerTSPSolution.txt  superchargerTSP.txt> /dev/null
!linkern -o superChargerTSPSolution.txt  superchargerTSP.txt> /dev/null
endTime = toc;
disp('-------------------------------------------------------');
disp(endTime);
disp('-------------------------------------------------------');


%reading the output file of linkern
in = fopen('superChargerTSPSolution.txt');

%we read the number of points on the tsp
ntsp = fscanf(in,'%d',2);

ntsp = ntsp(1); %we get the first digit (the other one is the same)

%reading all the other information (the first line contains the sequence of points in the tsp)

ord = fscanf(in,'%d',[3,ntsp]);

fclose(in);

ord = ord';

ord(:,1:2) = ord(:,1:2)+1; % adding 1 to indexes (Likern starts at 0)

%atsp_ord = ord(1:2:end,:);

vertexSequenceTemp = ord(:,1); %selecting the first line

%disp('Edge Order');
%disp(atsp_ord');
% Pruning the Augmented Symmetric TSP Solution Down to the ATSP Solution

vertexSequence = vertexSequenceTemp(1:2:end);
%disp(' ');
%disp('Asymmetric TSP Vertex Sequence');
%disp(vertexSequence);


% sort to make sure sets are not dispersed

vertexSequenceOrdered = vertexSequence;

if ( setSorted(vertexSequence(1)) == setSorted(vertexSequence(end)))
   
    vertexSequenceOrdered = [];
    currentSet = setSorted(vertexSequence(1));
    idx = 1;
    goOn = 1;
    while( goOn == 1)
       
        if( setSorted(vertexSequence(idx)) ~= currentSet)
            goOn = 0;
        else
         idx = idx+1;
        end
    end
    
     vertexSequenceOrdered = [vertexSequence(idx:end); vertexSequence(1:idx-1)];
        
end





end