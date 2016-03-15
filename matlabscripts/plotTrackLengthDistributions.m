%For Exp Set 1
%ExpCondTitles = {' OR',' GC',' AB',' OR',' GC',' AB',' OR',' GC',' AB'};
%For exp Set 2
ExpCondTitles = {' ATTP40',' B47',' B48'};
ExpCondFood = {'0.0% DMSO','0.0% DMSO','0.0% DMSO','0.5% DMSO','0.5% DMSO','0.5% DMSO','1.0% DMSO','1.0% DMSO','1.0% DMSO'};

%% Box plot of track lengths at specific Time
goToHour        = 85; %90 is arround the usual recording time ~ VialAge(1)
t               = (goToHour*3600- VialAge(1))/timeAdvance; %(goToHour*3600 - VialAge(1))/timeAdvance;
ExpTrackResults = ExpTrackResultsInTime{t};

ConditionIndex = 1;
for ConditionIndex=1:ConditionIndexMax
    ExpTrackResults         = ExpTrackResultsInTime{t};
    ResSet                  = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    %ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    
end

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)));
groups = [ zeros( length(meanConditionLength{1}) ,1); ones(length(meanConditionLength{2}),1); 2*ones(length(meanConditionLength{3}),1) ];
boxplot([meanConditionLength{1};meanConditionLength{2};meanConditionLength{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
title(strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)))
ylim([0 250]);
saveas(hf,sprintf('figures/NFTracklet%sLengthBoxPlot-%dHour.png',strOutputTag,goToHour));

if ConditionIndexMax > 3  
    ConditionIndex = 4;
    hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)));
    groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
    boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
    title(strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)))
    ylim([0 250]);
    saveas(hf,sprintf('figures/DMSO05Tracklet%sLengthBoxPlot-%dHour.png',strOutputTag,goToHour));
end
 
if (ConditionIndexMax >= 7)
    ConditionIndex = 7;
    hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)));
    groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
    boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
    title(strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths hour@',num2str(goToHour)))
    ylim([0 250]);
    saveas(hf,sprintf('figures/DMSO10Tracklet%sLengthBoxPlot-%dHour.png',strOutputTag,goToHour));
end

%% Across all time
FromFrame = 0;
TimeWindow = max(ExpTrack{1}(:,1));
ExpTrackResultsAllTime = ExtractFilteredTrackData(ExpTrack,ExpIDs,framePeriod,MinLifetime, MaxLifetime, MinDistance, MaxStepLength, FromFrame ,TimeWindow, MinStepLength ,bVerbose);
for ConditionIndex=1:ConditionIndexMax
    ResSet                  = vertcat(ExpTrackResultsAllTime{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    %ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    
end

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths All time'));
groups = [ zeros( length(meanConditionLength{1}) ,1); ones(length(meanConditionLength{2}),1); 2*ones(length(meanConditionLength{3}),1) ];
boxplot([meanConditionLength{1};meanConditionLength{2};meanConditionLength{3}],groups,'labels',{strcat(ExpCondFood{1},ExpCondTitles{1}),strcat(ExpCondFood{2},ExpCondTitles{2}),strcat(ExpCondFood{3},ExpCondTitles{3})})
title('Track length across all time')
ylim([0 250]);
saveas(hf,sprintf('figures/NFTracklet%sLengthBoxPlot-from%dtoEnd.png',strOutputTag,goToHour));

if ConditionIndexMax > 3  
   
    ConditionIndex = 4;
    hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths All time'));
    groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
    boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
    title('Track length across all time')
    ylim([0 250]);
    saveas(hf,sprintf('figures/DMSO05Tracklet%sLengthBoxPlot-Allt.png',strOutputTag,goToHour));
end
    
if (ConditionIndexMax > 6)
    
    ConditionIndex = 7;
    hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Run  tracklet lengths All time'));
    groups = [ zeros( length(meanConditionLength{ConditionIndex+0}) ,1); ones(length(meanConditionLength{ConditionIndex+1}),1); 2*ones(length(meanConditionLength{ConditionIndex+2}),1) ];
    boxplot([meanConditionLength{ConditionIndex+0};meanConditionLength{ConditionIndex+1};meanConditionLength{ConditionIndex+2}],groups,'labels',{strcat(ExpCondFood{ConditionIndex+0},ExpCondTitles{ConditionIndex+0}),strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})})
    title('Track length across all time')
    ylim([0 250]);
    saveas(hf,sprintf('figures/DMSO10Tracklet%sLengthBoxPlot-Allt.png',strOutputTag,goToHour));
end
%% Plot Path Distance Distribution

nbins = 100;

ConditionIndex = 1;
hf = figure('Name',strcat(ExpCondFood{ConditionIndex}, ' Path Length Distribution in pixels'));

subplot(3,1,1);

ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
[cnt,bin] = hist(meanConditionLength{ConditionIndex},nbins);
hist(meanConditionLength{ConditionIndex},nbins)
title(strcat(ExpCondTitles{ConditionIndex},'  \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
%xlabel('px distance');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
xlimits = max(ceil(bin/100))*100;
ylimits =  ceil(max(cnt)/10)*10;
ylim([0 ylimits]);
xlim([-1 xlimits]);

if ConditionIndexMax > 1 
    subplot(3,1,2);
    ConditionIndex = 2;
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    %xlabel('px distance');
    ylim([0 ylimits]);
    xlim([0 xlimits]);
end

if ConditionIndexMax > 2 
    subplot(3,1,3);
    ConditionIndex = 3;
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex}, '  \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    %xlabel('px distance');
    ylim([0 ylimits]);
    xlim([0 xlimits]);
    saveas(hf,sprintf('figures/NFTrackLength%sHist-%dHour.png',strOutputTag,goToHour));
end
%%
if ConditionIndexMax > 3 
    ConditionIndex = 4;
    hf = figure('Name',ExpCondFood{ConditionIndex},' Path Length Distribution ');
    subplot(3,1,1);
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
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
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);

    subplot(3,1,3);
    ConditionIndex = 6;
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);
    saveas(hf,sprintf('figures/DMSO05TrackLength%sHist-%dHour.png',strOutputTag,goToHour));
end
%%
if ConditionIndexMax > 6 
    ConditionIndex = 7;
    hf = figure('Name',ExpCondFood{ConditionIndex}, ' Path Length Distribution ');
    subplot(3,1,1);
    
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
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
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);

    subplot(3,1,3);
    ConditionIndex = 9;
    ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
    meanConditionLength{ConditionIndex}  = vertcat(ResSet.Length);
    hist(meanConditionLength{ConditionIndex},nbins);
    title(strcat(ExpCondTitles{ConditionIndex},ExpCondFood{ConditionIndex},' \mu:',num2str(mean(meanConditionLength{ConditionIndex}))));
    ylim([0 ylimits]);
    xlim([0 xlimits]);
    saveas(hf,sprintf('figures/DMSO10TrackLength%sHist-%dHour.png',strOutputTag,goToHour));
end
