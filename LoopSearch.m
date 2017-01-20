function []=LoopSearch(periodID)
% Loop on CandidateAnalysis_v2 for al runs in period periodId
% OMH 18/06/2013

%% Build mrun matrix of all periods
periods = zeros(10,2);

periods(1,:) = [2538 2585];  % 
periods(2,:) = [2685 2890];  %
periods(3,:) = [3000 3086];  % 
periods(4,:) = [3157 3256];  %
periods(5,:) = [3336 3371];  % 
periods(6,:) = [3562 3733];  % New DAQ
periods(7,:) = [3835 3999];  % 
periods(8,:) = [4183 4389];  % 
periods(9,:) = [4444 4844];  % NS polar
periods(10,:) = [4845 5070]; % NS polar + DAQ upgrade
periods(11,:) = [5071 5227]; % Same (to be merged with period 10)

%% Loop on runs of period
runstart = periods(periodID,1);
runstop = periods(periodID,2);
disp(sprintf('Scanning period %d: R%d-%d',periodID,runstart,runstop))

for i=runstart:runstop
    CandidateAnalysis_v2(i)
end


