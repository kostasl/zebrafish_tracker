%% Import CSV Files

%%FOLDER NAMES SHOULD BE OF THE FORMAT: EXP_6-7_20151123_5sec
%5sec is the timelapse period in seconds

scriptPath = which('processFiles.m');
[framePeriod,VialAge,ExpIDs,ExpN ] = importCSVtoCell( '*V*_N','EXP*' );
save('LarvaCountData.mat','ExpN');


%% Plot Results for Active Larva Per Vial
%Give 3 days data points 1 sec each.
% Genotypes are 3 organized in this order : 1st WT (oregonR), 2nd Genetic Control, 3rd AlfaBeta Mutant
ConditionIndex = 1; %Experimental Condition ID Say 1 OR NF 
% The videos have 2 rows of 9 vials - Vials 1-10 have identical conditions so they go in PAIRS
VialPairsPerCondition = [[1,10];[2,11];[3,12];[4,13];[5,14];[6,15];[7,16];[8,17];[9,18]]; %OR Normal Food
timePoints = max(VialAge) + 24*3*3600;%Total Time points in seconds over which to analyse data
%FramePeriod sampled at each timelapse Experiment -

%%THESE DO NOT CORRESPOND TO TIMES OF VIDEOS
%framePeriod = [20;5;5;2;2;2;2;2];
display('The loaded videos had frame Periods designed in folder name:')
display(framePeriod);

%N = length(timePoints);
%datV{iVial} = zeros([size(ExpN,1),timePoints]);
%Copies all Exp data from the imported in the Cells - to a matrix of fixed
%time length with the data points place at the right timepoints

for (ConditionIndex=1:9)
    VialPairs = VialPairsPerCondition(ConditionIndex,:);
    [ResultsSourceRefIndex,datNV] = collectResultsInTimeVector( ExpN,VialPairs,VialAge,framePeriod,timePoints );
    %Make Cell Array of Data organized in Structures
    %Could Start Mean from min(VialAge) so to exclude lots of empty indexes
    meanNLarva = mean(datV{ConditionIndex}(:,1:timePoints));
    %STD Error - (Std dev normalized by sample size
    stdNLarva = std(datV{ConditionIndex}(:,1:timePoints),1,1) ;
    
    %Make structure for each experimental condition such that it contains
    %raw & processed data but also a table we can refer to Which Experiment
    %and Vial the data was sourced
    vialtrackerResults{ConditionIndex} = struct('ActiveCount',datNV,'MeanActiveCount',meanNLarva,'MeanErrCount',stdNLarva,'DataRowSourceIndex',ResultsSourceRefIndex);
end



%% PLOT MEAN Active LARVAE
h=figure('Name','Mean Number of Larva Moving NF');
plot((1:timePoints)/3600,vialtrackerResults{1}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{2}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{3}.MeanActiveCount);
xlabel('hours');
ylabel('N active larva');
legend('NF OR','NF CT','NF AB');
saveas(h,'figures/NFMeanNActive.pdf');

h=figure('Name','Mean Number of Larva Moving DMSO 0.5%');
plot((1:timePoints)/3600,vialtrackerResults{4}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{5}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{6}.MeanActiveCount);
xlabel('hours');
ylabel('N active larva');
legend('0.5% OR','0.5% CT','0.5% AB');
saveas(h,'figures/DMSO05NActive.pdf');

h=figure('Name','Mean Number of Larva Moving DMSO 1.0%');
plot((1:timePoints)/3600,vialtrackerResults{7}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{8}.MeanActiveCount,(1:timePoints)/3600,vialtrackerResults{9}.MeanActiveCount);
xlabel('hours');
ylabel('N active larva');
legend('1% OR','1% CT','1% AB');
saveas(h,'figures/DMSO10NActive.pdf');

%% PLOT STD Error Of MEAN
h= figure('Name','STD Error in Mean Number of Larva Moving NF');;
plot((1:timePoints)/3600,vialtrackerResults{1}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{2}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{3}.MeanErrCount);
xlabel('hours');
ylabel('N active larva');
legend('NF OR','NF CT','NF AB');
saveas(h,'figures/NFSTDErrorMeanNActive.pdf');


h=figure('Name','STD Error in Mean Number of Larva Moving DMSO 0.5%');;
plot((1:timePoints)/3600,vialtrackerResults{4}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{5}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{6}.MeanErrCount);
xlabel('hours');
ylabel('N active larva');
legend('0.5% OR','0.5% CT','0.5% AB');
saveas(h,'figures/DMSO05STDErrorMeanNActive.pdf');

h=figure('Name','STD Error in Mean Number of Larva Moving DMSO 1.0%');;
plot((1:timePoints)/3600,vialtrackerResults{7}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{8}.MeanErrCount,(1:timePoints)/3600,vialtrackerResults{9}.MeanErrCount);
xlabel('hours');
ylabel('N active larva');
legend('1% OR','1% CT','1% AB');
saveas(h,'figures/DMSO10STDErrorMeanNActive.pdf');

%% Plot Control Vials In High DMSO 1%
ConditionIndex = 8;


span = 2000; %10K points is 2.5 hours
h=figure('Name',strcat('Individual Vials CT with 1%DMSO Smoothed ',num2str(2000/3600), ' Hours'));
plot((1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(1,1:timePoints),span),...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(2,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(3,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(4,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(5,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(6,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(7,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(8,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(vialtrackerResults{ConditionIndex}.ActiveCount(9,1:timePoints),span),'*' ...
);
xlabel('hours');
ylabel('N active larva');
%%CT Vial with 1%DMSO Legends
legend(strcat(num2str(vialtrackerResults{ConditionIndex}.DataRowSourceIndex(1,:)),'V'),strcat(ExpIDs{2},'V'),strcat(ExpIDs{3},'V'),strcat(ExpIDs{4},'V'),strcat(ExpIDs{5},'V')); %,strcat(ExpIDs{6},'V'),strcat(ExpIDs{7},'V'),strcat(ExpIDs{8},'V'),strcat(ExpIDs{9},'V')

saveas(h,'figures/CTIndividualVialsNSmoothed.pdf');

% ALL TOGETHER
%% PLOT MEAN Active LARVAE
figure('Name','Mean Number of Larva Moving ');
plot((1:timePoints)/3600,smooth(meanNLarva{1},2000),'.', ...
     (1:timePoints)/3600,smooth(meanNLarva{2},2000),'-', ...
     (1:timePoints)/3600,smooth(meanNLarva{3},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{4},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{5},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{6},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{7},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{8},2000), ...
     (1:timePoints)/3600,smooth(meanNLarva{9},2000) ...
 ) ;
xlabel('hours');
ylabel('N active larva');
legend('NF OR','NF CT','NF AB');

%figure('Name','Error Bars');
%errorbar((1:timePoints)/3600,meanNLarva{ConditionIndex},stdNLarva{ConditionIndex},'r')

% 
% %%%% NOTES/SCRAPBOOK %%%
% 
% hold on;
% plot((1:timePoints)/3600,datV{ConditionIndex}([7;8],1:timePoints))
% plot((1:timePoints)/3600,datV{iVial}(3,1:timePoints))
% plot((1:timePoints)/3600,datV{iVial}(4,1:timePoints))
% plot((1:timePoints)/3600,datV{iVial}(5,1:timePoints))
% 
% datV{iVial}()
% % 
% 
% medfilt1(ExpN{1, 1}(1:N,3),100);
% timeE{1} = ExpN{1, 1}(1:N,1)/(3600/5);
% datV1{2} = medfilt1(ExpN{2, 1}(1:N,3),100);
% timeE{2} = ExpN{1, 1}(1:N,1)/(3600/5);
% datV1{3} = medfilt1(ExpN{3, 1}(:,3),100);
% timeE{3} = ExpN{1, 1}(:,1)/(3600/2);
% 
% plot(timeE{1},datV1{1},'.','DisplayName','Exp2 V1');
% hold on;
% plot(timeE{2},datV1{2},'x','DisplayName','Exp3 V1');
% %plot(timeE{3},datV1{3},':','DisplayName','Exp4 V1');
