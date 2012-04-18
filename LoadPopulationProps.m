function LoadPopulationProps(j)

global modelAtm %NOTHING ON WALL

%Pop = modelAtm.Pop;

    PopString = int2str(j);
     
   %For now, assume that all populations have properties of SOA
    eval(['modelAtm.Pop' PopString '.lambda = modelAtm.SOA.lambda;']);
    eval(['modelAtm.Pop' PopString '.alpha = modelAtm.SOA.alpha;']);
    eval(['modelAtm.Pop' PopString '.Diff = modelAtm.SOA.Diff;']);
    eval(['modelAtm.Pop' PopString '.va = modelAtm.SOA.va;']);
    eval(['modelAtm.Pop' PopString '.rho = modelAtm.SOA.rho;']);
    %eval(['modelAtm.Pop' PopString '.coverage = 0;']);
  
%    if j==1 
        eval(['modelAtm.Pop' PopString '.shapeFact = 1;']);
        eval(['modelAtm.Pop' PopString '.wallYN = 0;']);
%     else
%         eval(['modelAtm.Pop' PopString '.shapeFact = 1;']);
%        % eval(['modelAtm.Pop' PopString '.shapeFact = 2^(-1/3);']);
%         %eval(['modelAtm.Pop' PopString '.wallYN = 1;']);
%         eval(['modelAtm.Pop' PopString '.wallYN = 0;']);
%     
%     
%     
%         eval(['NumConc = modelAtm.Pop' PopString '.NumConc0;']);
%         eval(['Dp = modelAtm.Pop' PopString '.Dp0;']);
%         eval(['shapeFact = modelAtm.Pop' PopString '.shapeFact;']);
%         n_a = NumConc/modelAtm.StoV; %#/m2
%         coverage = n_a*pi*Dp^2*shapeFact;
%         eval(['modelAtm.Pop' PopString '.coverage = coverage']);
%     end

    
    