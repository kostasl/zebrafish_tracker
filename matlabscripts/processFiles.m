%% Import CSV Files


scriptPath = which('processFiles.m');
ExpN = importCSVtoCell( '*V*_N','EXP*' );
save('LarvaCountData.mat','ExpN');


%% Plot Results for Active Larva Per Vial
%Give 3 days data points 1 sec each.
% Genotypes are 3 organized in this order : 1st WT (oregonR), 2nd Genetic Control, 3rd AlfaBeta Mutant
ConditionIndex = 1; %Experimental Condition ID Say 1 OR NF 
% The videos have 2 rows of 9 vials - Vials 1-10 have identical conditions so they go in PAIRS
VialPairsPerCondition = [[1,10];[2,11];[3,12];[4,13];[5,14];[6,15];[7,16];[8,17];[9,18]]; %OR Normal Food
timePoints = 24*3*3600;%Total Time points in seconds over which to analyse data
%FramePeriod sampled at each timelapse Experiment
framePeriod = [20;5;5;2;2;2;2;2];

%N = length(timePoints);
%datV{iVial} = zeros([size(ExpN,1),timePoints]);
%Copies all Exp data from the imported in the Cells - to a matrix of fixed
%time length with the data points place at the right timepoints

for (ConditionIndex=1:9)
    VialPairs = VialPairsPerCondition(ConditionIndex,:);
    datV{ConditionIndex} = collectResultsInTimeVector( ExpN,VialPairs,framePeriod,timePoints );
    meanNLarva{ConditionIndex} = mean(datV{ConditionIndex}(:,1:timePoints));
    %STD Error - (Std dev normalized by sample size
    stdNLarva{ConditionIndex} = std(datV{ConditionIndex}(:,1:timePoints),1,1) ;
end

%% PLOT MEAN Active LARVAE
figure('Name','Mean Number of Larva Moving ');
plot((1:timePoints)/3600,meanNLarva{1},(1:timePoints)/3600,meanNLarva{2},(1:timePoints)/3600,meanNLarva{3});
xlabel('hours');
ylabel('N active larva');
legend('NF OR','NF CT','NF AB');

figure('Name','Mean Number of Larva Moving ');
plot((1:timePoints)/3600,meanNLarva{4},(1:timePoints)/3600,meanNLarva{5},(1:timePoints)/3600,meanNLarva{6});
xlabel('hours');
ylabel('N active larva');
legend('0.5% OR','0.5% CT','0.5% AB');


figure('Name','Mean Number of Larva Moving ');
plot((1:timePoints)/3600,meanNLarva{7},(1:timePoints)/3600,meanNLarva{8},(1:timePoints)/3600,meanNLarva{9});
xlabel('hours');
ylabel('N active larva');
legend('1% OR','1% CT','1% AB');


%% Plot Control Vials In High DMSO 1%
ConditionIndex = 8;

h=figure('Name','Individual Vials CT with 1%DMSO Smoothed 2.5Hours');
span = 20000; %10K points is 2.5 hours
plot((1:timePoints)/3600,smooth(datV{ConditionIndex}(1,1:timePoints),span),...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(2,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(3,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(4,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(5,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(6,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(7,1:timePoints),span), ...
        (1:timePoints)/3600,smooth(datV{ConditionIndex}(8,1:timePoints),span), ...
    (1:timePoints)/3600,smooth(datV{ConditionIndex}(9,1:timePoints),span) ...
);
xlabel('hours');
ylabel('N active larva');
%%CT Vial with 1%DMSO Legends
legend('Exp1 V8','Exp2 V8','Exp2 V17','Exp3 V8','Exp3 V17','Exp4 V8','Exp4 V17','Exp5 V8','Exp5 V17');

saveas(h,'CTIndividualVialsNSmoothed.pdf');

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

figure('Name','STD Error in Mean Number of Larva Moving ');;
plot((1:timePoints)/3600,stdNLarva{1},(1:timePoints)/3600,stdNLarva{2},(1:timePoints)/3600,stdNLarva{3});
xlabel('hours');
ylabel('N active larva');
legend('NF OR','NF CT','NF AB');


figure('Name','STD Error in Mean Number of Larva Moving ');;
plot((1:timePoints)/3600,stdNLarva{4},(1:timePoints)/3600,stdNLarva{5},(1:timePoints)/3600,stdNLarva{6});
xlabel('hours');
ylabel('N active larva');
legend('0.5% OR','0.5% CT','0.5% AB');


figure('Name','STD Error in Mean Number of Larva Moving ');;
plot((1:timePoints)/3600,stdNLarva{7},(1:timePoints)/3600,stdNLarva{8},(1:timePoints)/3600,stdNLarva{9});
xlabel('hours');
ylabel('N active larva');
legend('1% OR','1% CT','1% AB');



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
