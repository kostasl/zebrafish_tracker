ylimits = 500;
xlimits = 2;
nbins = 1000;

hf = figure('Name','Normal Food MEAN SPEED');
subplot(3,1,1);
ConditionIndex = 1;
ResSet                      = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hold off;
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('OR NF μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));h = findobj(gca,'Type','patch');
set(h,'FaceColor','blue');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 2;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('CT NF μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','blue');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 3;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB NF μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/NFTrackletSpeedHist.pdf')


hold off;
hf = figure('Name','0.5% DMSO MEAN SPEED');
subplot(3,1,1);
ConditionIndex = 4;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('OR 0.5% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 5;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('CT 0.5% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,3);
ConditionIndex = 6;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('AB 0.5% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO05TrackletSpeedHist.pdf')


%DMSO 1%
hold off;
hf = figure('Name','1% DMSO MEAN SPEED');
subplot(3,1,1);
ConditionIndex = 7;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title(strcat('OR 1% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
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
title(strcat('CT 1% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
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
title(strcat('AB 1% DMSO μ:',num2str(mean(meanConditionSpeeds{ConditionIndex}))));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO10TrackletSpeedHist.pdf')
