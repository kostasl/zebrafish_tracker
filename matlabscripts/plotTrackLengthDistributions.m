
ExpCondTitles = {' OR',' GC',' AB',' OR',' GC',' AB',' OR',' GC',' AB'};
ExpCondFood = {'0.0% DMSO','0.0% DMSO','0.0% DMSO','0.5% DMSO','0.5% DMSO','0.5% DMSO','1.0% DMSO','1.0% DMSO','1.0% DMSO'};

%% Box plot of track lengths

ConditionIndex = 1;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);


hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths '));
groups = [ zeros( length(meanConditionLength{1}) ,1); ones(length(meanConditionLength{2}),1); 2*ones(length(meanConditionLength{3}),1) ];
boxplot([meanConditionLength{1};meanConditionLength{2};meanConditionLength{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
ylim([0 250]);
saveas(hf,sprintf('figures/NFTrackletLengthBoxPlot-%dHour.png',goToHour));


%% Plot Path Distance Distribution

nbins = 100;

hf = figure('Name','NF Path Length Distribution in pixels');

subplot(3,1,1);
ConditionIndex = 1;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
[cnt,bin] = hist(meanConditionLength{ConditionIndex},nbins);
hist(meanConditionLength{ConditionIndex},nbins)
title(strcat('OR NF \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
%xlabel('px distance');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
xlimits = 400;
ylimits =  2*ceil(max(cnt)/10)*10;
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 2;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT NF \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
%xlabel('px distance');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 3;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('AB NF \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
%xlabel('px distance');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/NFTrackLengthHist-%dHour.png',goToHour));

%%
hf = figure('Name','DMSO 0.5% Path Length Distribution ');
subplot(3,1,1);
ConditionIndex = 4;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('OR DMSO 0.5% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 5;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT DMSO 0.5% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 6;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('AB DMSO 0.5% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/DMSO05TrackLengthHist-%dHour.png',goToHour));

%%
hf = figure('Name','DMSO 1% Path Length Distribution ');
subplot(3,1,1);
ConditionIndex = 7;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('OR DMSO 1% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 8;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT DMSO 1% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 9;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('AB DMSO 1% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/DMSO10TrackLengthHist-%dHour.png',goToHour));

