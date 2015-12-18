
ExpCondTitles = {' OR',' GC',' AB',' OR',' GC',' AB',' OR',' GC',' AB'};
ExpCondFood = {'0.0% DMSO','0.0% DMSO','0.0% DMSO','0.5% DMSO','0.5% DMSO','0.5% DMSO','1.0% DMSO','1.0% DMSO','1.0% DMSO'};

nbins = 100;

clear meanConditionSpeeds;
clear mu;
clear n;
clear stdd;
%% Do Mean Speed Per Condition Per Time Window

meanConditionSpeeds  = {};
n  = zeros(length(ExpTrackResultsInTime),9);
mu  = zeros(length(ExpTrackResultsInTime),9);
stdd  = zeros(length(ExpTrackResultsInTime),9);
ConditionIndex = 1;
for t=1:length(ExpTrackResultsInTime)
    
    for (ConditionIndex=1:3)
        ExpTrackResults = ExpTrackResultsInTime{t};
        ResSet                   = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
        
        if isempty(ResSet)
            continue;
        end
        meanConditionSpeeds{ConditionIndex}   = vertcat(ResSet.MeanSpeed);
        n(t,ConditionIndex)                   = length(meanConditionSpeeds{ConditionIndex});
        mu(t,ConditionIndex)                  = mean(meanConditionSpeeds{ConditionIndex});
        stdd(t,ConditionIndex)                = std(meanConditionSpeeds{ConditionIndex});
    end
end

title('Mean Activity in hour frames');
plot((1:t)*timeAdvance/3600,mu(:,1));
legend(strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3}))
xlabel('Hour');
%% Plot Histogram Of Speed Within A time Window
hold off;
ConditionIndex = 1;
strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'MEAN SPEED'));
subplot(3,1,1);

ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
n       = length(meanConditionSpeeds{ConditionIndex});
mu      = mean(meanConditionSpeeds{ConditionIndex});
stdd    = std(meanConditionSpeeds{ConditionIndex});
strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );
title(strTitle);
 ylim([0 ylimits]);
 xlim([0 xlimits]);
saveas(hf,'figures/NFTrackletSpeedHist.pdf')


%%Box Plot Of Mean Speeds per tracklet
ConditionIndex = 1;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Crawl-Run SPEEDs'));
groups = [ zeros( length(meanConditionSpeeds{1}) ,1); ones(length(meanConditionSpeeds{2}),1); 2*ones(length(meanConditionSpeeds{3}),1) ];
boxplot([meanConditionSpeeds{1};meanConditionSpeeds{2};meanConditionSpeeds{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
ylim([0 15]);
saveas(hf,'figures/NFTrackLetSpeedBoxPlot.pdf');
%% DMSO 0.5
%Check If Condition Exists

    
hold off;
ConditionIndex = 4;

if size(ExpTrackResults,2) >= VialPairsPerCondition(ConditionIndex )
    hf = figure('Name','0.5% DMSO MEAN SPEED');
    subplot(3,1,1);
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
    meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
    hist(meanConditionSpeeds{ConditionIndex},nbins);
    title(strcat('OR 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);
else
    
    error('Exp. Conditions Missing - Stopping Plots')
    return;
end

subplot(3,1,2);
ConditionIndex = 5;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('CT 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 6;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB 0.5% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO05TrackletSpeedHist.pdf')


%% DMSO 1%
hold off;
hf = figure('Name','1% DMSO MEAN SPEED');
title('1% DMSO MEAN SPEED');
subplot(3,1,1);
ConditionIndex = 7;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB 1% DMSO \mu:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO10TrackletSpeedHist.pdf')
