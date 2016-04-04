%Produces plot of Mean Speed / Activity based on Number of tracklet samples
%For EXp Set 1 - 9 conditions :clo
%ExpCondTitles = {' OR',' GC',' AB',' OR',' GC',' AB',' OR',' GC',' AB'};
%For EXp Set 2 - 9 conditions :
ExpCondTitles = {' ATTP40',' BWD47',' BWD48',' ATTP2',' 48.2',' 34'};
ExpCondFood = {'0.0% D','0.0% D','0.0% D','0.0% D','0.0% D','0.0% D','0.0% D','0.0% D','0.0% D'};
%ExpCondFood = {'0.0% DMSO','0.0% DMSO','0.0% DMSO','0.5% DMSO','0.5% DMSO','0.5% DMSO','1.0% DMSO','1.0% DMSO','1.0% DMSO'};

nbins = 100;
ylimits = 7;
ylimitsTracklets =  2500;
clear meanConditionSpeeds;
clear mu;
clear n;
clear stdd;

%cd /media/kostasl/FlashDrive/PilotVialTrack/ExpSet2_201603/DataOut %Home
%load(strcat('LarvaTrackData',strOutputTag,'.mat'));
%% Do Mean Speed Per Condition Per Time Window
meanConditionSpeeds  = {};
n  = zeros(length(ExpTrackResultsInTime),9);
mu  = zeros(length(ExpTrackResultsInTime),9);
stdd  = zeros(length(ExpTrackResultsInTime),9);
ConditionIndex = 1;


for t=1:length(ExpTrackResultsInTime)
    for (ConditionIndex=1:ConditionIndexMax)
            ExpTrackResults         = ExpTrackResultsInTime{t};
            ResSet                  = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});

            if isempty(ResSet)
                continue;
            end
            meanConditionSpeeds{ConditionIndex}   = vertcat(ResSet.MeanSpeed);
            n(t,ConditionIndex)                   = length(meanConditionSpeeds{ConditionIndex});
            mu(t,ConditionIndex)                  = mean(meanConditionSpeeds{ConditionIndex});
            stdd(t,ConditionIndex)                = std(meanConditionSpeeds{ConditionIndex});
    end
    
end

ylimitsTracklets = ceil(max(n(:))/100)*100;
ylimits = ceil(max(mu(:)));

strtitle = sprintf('Mean Activity -Sliding Window %d hours',TimeFrameWidth/3600);
% Plot Results - In Groups or - In sets of 3-genotypes For each Food Condition
CondGrouping = ConditionIndexMax; %When =ConditionIndexMax then plot them all together
for (ConditionIndex=1:CondGrouping:ConditionIndexMax)

    t = length(ExpTrackResultsInTime);
    Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;
    hf = figure('Name',strcat(ExpCondFood{ConditionIndex},strtitle));
    
    subplot(3,1,1)
    plot(Exptime,mu(:,ConditionIndex:(ConditionIndex+CondGrouping-1))   );
    title('Mean Speed in px/sec');
    ylim([0 ylimits]);
    subplot(3,1,2)
    plot(Exptime,n(:,ConditionIndex:(ConditionIndex+CondGrouping-1)));
    title('Number of samples');
    hh= subplot(3,1,3)
    %get(hh,'position')


    plot(Exptime,stdd(:,ConditionIndex:(ConditionIndex+CondGrouping-1)));
    title('STD Dev ');
    xlabel('Hour');
    ylim([0 ylimits]);
    
    %Make Legend
    strLegend = '';
    j = 0;
    for i=ConditionIndex:(ConditionIndex+CondGrouping-1)
        j=j+1; 
        strLegend{j} = strcat(ExpCondFood{i},ExpCondTitles{i});
    end
    legend( strLegend,'Location','southoutside','Orientation','horizontal')
    set(hh,'position',[0.13 0.2 0.77 0.12]); %Fix Last plot after adding legends
    saveas(hf,strcat('figures/meanVial',strOutputTag,'SpeedSlidingWindow',ExpCondFood{ConditionIndex},'.png'));
end
%,strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})

ConditionIndex =1;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},strtitle));
subplot(3,1,1)
    plot(Exptime,n(:,ConditionIndex+0),Exptime,n(:,ConditionIndex+1),Exptime,n(:,ConditionIndex+2));
    
    title(ExpCondFood{ConditionIndex});
    ylim([0 ylimitsTracklets ]);
    
if ConditionIndexMax > 3
    ConditionIndex =4;
    subplot(3,1,2)
    plot(Exptime,n(:,ConditionIndex+0),Exptime,n(:,ConditionIndex+1),Exptime,n(:,ConditionIndex+2));
    ylim([0 ylimitsTracklets ]);
    legend(strcat(ExpCondTitles{ConditionIndex+0}),strcat(ExpCondTitles{ConditionIndex+1}),strcat(ExpCondTitles{ConditionIndex+2}))
    title(ExpCondFood{ConditionIndex});
end

if ConditionIndexMax > 6
    ConditionIndex =7;
    subplot(3,1,3)
    plot(Exptime,n(:,ConditionIndex+0),Exptime,n(:,ConditionIndex+1),Exptime,n(:,ConditionIndex+2));
    ylim([0 ylimitsTracklets ]);
    legend(strcat(ExpCondTitles{ConditionIndex+0}),strcat(ExpCondTitles{ConditionIndex+1}),strcat(ExpCondTitles{ConditionIndex+2}))
    title(ExpCondFood{ConditionIndex});
end

legend(strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2}),'Location','southoutside','Orientation','horizontal')
saveas(hf,strcat('figures/meanVial',strOutputTag,'ActivityFor3FoodConditions.png'));


return;

%% Plot Histogram Of Speed Within A chosen time Window

t= (goToHour*3600 - VialAge(1))/timeAdvance;
ExpTrackResults = ExpTrackResultsInTime{t};

hold off;
ConditionIndex = 1;
strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'MEAN SPEED'));

subplot(3,1,1);
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
[cnt,bin]                            = hist(meanConditionSpeeds{ConditionIndex},nbins);
hist(meanConditionSpeeds{ConditionIndex},nbins);
n       = length(meanConditionSpeeds{ConditionIndex});
mu      = mean(meanConditionSpeeds{ConditionIndex});
stdd    = std(meanConditionSpeeds{ConditionIndex});
title(strTitle);
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','blue');
% set(h,'EdgeColor','w');

ylimits =  2*ceil(max(cnt)/10)*10;
xlimits = 10;
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 2;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
n       = length(meanConditionSpeeds{ConditionIndex});
mu      = mean(meanConditionSpeeds{ConditionIndex});
stdd    = std(meanConditionSpeeds{ConditionIndex});
strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );
title(strTitle);
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 3;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
n       = length(meanConditionSpeeds{ConditionIndex});
mu      = mean(meanConditionSpeeds{ConditionIndex});
stdd    = std(meanConditionSpeeds{ConditionIndex});
strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );
title(strTitle);
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/NFTracklet%sSpeedHist-%dHour.png',strOutputTag,goToHour))


%%Box Plot Of Mean Speeds per tracklet
ConditionIndex = 1;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Crawl-Run SPEEDs'));
groups = [ zeros( length(meanConditionSpeeds{1}) ,1); ones(length(meanConditionSpeeds{2}),1); 2*ones(length(meanConditionSpeeds{3}),1) ];
boxplot([meanConditionSpeeds{1};meanConditionSpeeds{2};meanConditionSpeeds{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
ylim([0 15]);
saveas(hf,sprintf('figures/NFTracklet%sSpeedBoxPlot-%dHour.png',strOutputTag,goToHour));


%% DMSO 0.5
%Check If Condition Exists

    
hold off;
ConditionIndex = 4;

if size(ExpTrackResults,2) >= VialPairsPerCondition(ConditionIndex,: )
    hf = figure('Name','0.5% DMSO MEAN SPEED');
    subplot(3,1,1);
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
    hist(meanConditionSpeeds{ConditionIndex},nbins);
    title(strcat('OR 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);
else
    
    display('Exp. Conditions Missing - Stopping Plots')
    return;
end

subplot(3,1,2);
ConditionIndex = 5;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('CT 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 6;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/DMSO05Tracklet%sSpeedHist.pdf',strOutputTag))


%% DMSO 1%
hold off;
hf = figure('Name','1% DMSO MEAN SPEED');
title('1% DMSO MEAN SPEED');
subplot(3,1,1);
ConditionIndex = 7;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('OR 1% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 8;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('CT 1% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 9;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB 1% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,strcat('figures/DMSO10Tracklet',strOutputTag,'SpeedHist.pdf'));
