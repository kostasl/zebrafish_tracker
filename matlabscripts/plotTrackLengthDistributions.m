
ExpCondTitles = {' OR',' GC',' AB',' OR',' GC',' AB',' OR',' GC',' AB'};
ExpCondFood = {'0.0% DMSO','0.0% DMSO','0.0% DMSO','0.5% DMSO','0.5% DMSO','0.5% DMSO','1.0% DMSO','1.0% DMSO','1.0% DMSO'};

%% Box plot of track lengths at specific Time
goToHour        = 110; %90 is arround the usual recording time ~ VialAge(1)
t               = (goToHour*3600- VialAge(1))/timeAdvance; %(goToHour*3600 - VialAge(1))/timeAdvance;
ExpTrackResults = ExpTrackResultsInTime{t};

ConditionIndex = 1;
for ConditionIndex=1:9
    ExpTrackResults         = ExpTrackResultsInTime{t};
    ResSet                  = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    %ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    
end

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)));
groups = [ zeros( length(meanConditionLength{1}) ,1); ones(length(meanConditionLength{2}),1); 2*ones(length(meanConditionLength{3}),1) ];
boxplot([meanConditionLength{1};meanConditionLength{2};meanConditionLength{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
title(strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)))
ylim([0 250]);
saveas(hf,sprintf('figures/NFTrackletLengthBoxPlot-%dHour.png',goToHour));

if ConditionIndexMax < 4  break;

ConditionIndex = 4;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)));
groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
title(strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)))
ylim([0 250]);
saveas(hf,sprintf('figures/DMSO05TrackletLengthBoxPlot-%dHour.png',goToHour));

if ConditionIndexMax < 7  break;
ConditionIndex = 7;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)));
groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
title(strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths hour@',num2str(goToHour)))
ylim([0 250]);
saveas(hf,sprintf('figures/DMSO10TrackletLengthBoxPlot-%dHour.png',goToHour));

%Across all time
ExpTrackResultsAllTime = ExtractFilteredTrackData(ExpTrack,ExpIDs,framePeriod,MinLifetime, MaxLifetime, MinDistance, MaxStepLength, TimeFrameWidth ,maxRecordingTime, MinStepLength ,bVerbose);
for ConditionIndex=1:9
    ResSet                  = vertcat(ExpTrackResultsAllTime{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    %ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    
end

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths All time'));
groups = [ zeros( length(meanConditionLength{1}) ,1); ones(length(meanConditionLength{2}),1); 2*ones(length(meanConditionLength{3}),1) ];
boxplot([meanConditionLength{1};meanConditionLength{2};meanConditionLength{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
title('Track length across all time')
ylim([0 250]);
saveas(hf,sprintf('figures/NFTrackletLengthBoxPlot-Allt.png',goToHour));

if ConditionIndexMax < 4  break;
ConditionIndex = 4;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths All time'));
groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
title('Track length across all time')
ylim([0 250]);
saveas(hf,sprintf('figures/DMSO05TrackletLengthBoxPlot-Allt.png',goToHour));

if ConditionIndexMax < 7  break;
ConditionIndex = 7;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'Run  tracklet lengths All time'));
groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
title('Track length across all time')
ylim([0 250]);
saveas(hf,sprintf('figures/DMSO10TrackletLengthBoxPlot-Allt.png',goToHour));


%% Plot Path Distance Distribution

nbins = 100;

hf = figure('Name','NF Path Length Distribution in pixels');

subplot(3,1,1);
ConditionIndex = 1;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT NF \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
%xlabel('px distance');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 3;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT DMSO 0.5% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 6;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
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
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('CT DMSO 1% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 9;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
hist(meanConditionLength{ConditionIndex},nbins);
title(strcat('AB DMSO 1% \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,sprintf('figures/DMSO10TrackLengthHist-%dHour.png',goToHour));

