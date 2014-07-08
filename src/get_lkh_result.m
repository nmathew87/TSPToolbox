function vertexSequenceOrdered = get_lkh_result(symtsp,setSorted)

% creating params file for LKH

fidparam = fopen('params.txt','wt');
fprintf(fidparam, 'PROBLEM_FILE = lkhTSPinput.tsp\n');
fprintf(fidparam, 'OUTPUT_TOUR_FILE = lkhTSPoutput.out\n');
fprintf(fidparam, 'KICKS = 200\n');
fprintf(fidparam, 'KICK_TYPE = 6\n');
fprintf(fidparam, 'INITIAL_TOUR_ALGORITHM = NEAREST-NEIGHBOR\n');


%creating the input file for LKH

n = size(symtsp,1);

fid=fopen('lkhTSPinput.tsp', 'wt');

fprintf(fid,'TYPE: ATSP\nDIMENSION:%7.0f\n',n);

fprintf(fid,'EDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT:FULL_MATRIX\nEDGE_WEIGHT_SECTION\n');


for i=1:n
    fprintf(fid, '%7.0f ', round(symtsp(i,:)));
    fprintf(fid, '\n');
end

fprintf(fid, 'EOF');
fclose(fid);

tic;

% run LKH solver
!LKH params.txt

endTime = toc;
disp('-------------------------------------------------------');
disp(endTime);
disp('-------------------------------------------------------');


% parse LKH output

vertexSequenceOrdered = [];

%reading the output file of linkern
fin = fopen('lkhTSPoutput.out');

% reach tour section
cont = 1;
while( cont == 1)
    text = fscanf(fin, '%s',1);
    if(strcmp(text, 'TOUR_SECTION') == 1)
        cont = 0;
    end
end

%reading all the other information (the first line contains the sequence of points in the tsp)
sizeA = size(symtsp,1);

vertexSequenceTemp = fscanf(fin,'%d',[1,sizeA]);

fclose(fin);


%disp('Edge Order');
%disp(atsp_ord');
% Pruning the Augmented Symmetric TSP Solution Down to the ATSP Solution

%vertexSequence = vertexSequenceTemp(1:2:end);
vertexSequence = vertexSequenceTemp';


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