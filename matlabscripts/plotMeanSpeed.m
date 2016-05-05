%Produces plot of Mean Speed / Activity based on Number of tracklet samples
%Add Centroids To Plots And Save Centroid_TAG variable

%Conditions Label Are in  process FileTracks
nbins = 100;
ylimits = 7;
ylimitsTracklets =  2500;
clear meanConditionSpeeds;
clear mu;
clear n;
clear stdd;

plotcoloursPerVial = [1,0,0; ...
                      1,0,0; ...
                      1,0,0; ...
                      0.5,0.3,0.1; ...
                      0.5,0.3,0.1; ...
                      0.5,0.3,0.1; ...
                      0,0,1; ...
                      0,0,1; ...
                      0,0,1; ...
                      0,0,0; ...
                      0,0,0; ...
                      0,0,0; ...
                      1,0,1; ...
                      1,0,1; ...
                      1,0,1; ...
                      0.3,1,0; ...
                      0.3,1,0; ...
                      0.3,1,0; ...
                      ];

plotcoloursPerCondition = [1,0,0; ...
                      0.5,0.3,0.1; ...
                      0,0,1; ...
                      0,0,0; ...
                      1,0,1; ...
                      0.3,1,0; ...
                      ];

                  %cd /media/kostasl/FlashDrive/PilotVialTrack/ExpSet2_201603/DataOut %Home
%load(strcat('LarvaTrackData',strOutputTag,'.mat'));

%%Collect Centroids
%centroid = [centroid_R6_(1,:); centroid_R7_(1,:);centroid_R8_(1,:); centroid_R9_(1,:); centroid_R10_(1,:);]

%% Calc Data Per Vial Independently %%
meanConditionSpeedsV  = {};
maxVialCount = 18;
nV  = zeros(length(ExpTrackResultsInTime),maxVialCount);
muV  = zeros(length(ExpTrackResultsInTime),maxVialCount);
stddV  = zeros(length(ExpTrackResultsInTime),maxVialCount);
ConditionIndex = 1;

for (VialIndex=1:1:maxVialCount)
    for t=1:length(ExpTrackResultsInTime)
    
            ExpTrackResults         = ExpTrackResultsInTime{t};
            %TODO: Add Filters Here
            %ExpTrackResults =             
            
            ResSet                  = vertcat(ExpTrackResults{:,VialIndex});

            if isempty(ResSet)
                continue;
            end
            meanConditionSpeedsV{VialIndex}   = vertcat(ResSet.MeanSpeed);
            nV(t,VialIndex)                   = length(meanConditionSpeedsV{VialIndex});
            muV(t,VialIndex)                  = mean(meanConditionSpeedsV{VialIndex});
            stddV(t,VialIndex)                = std(meanConditionSpeedsV{VialIndex});

    end
           %Calc Centroids
           centrTime = sum((1:t)'.*nV(:,VialIndex))/sum( nV(:,VialIndex)); %Time Cntr X
           if ~isnan(centrTime)
                tcV(VialIndex) =  Exptime( max( round(centrTime ),1 )) ;
           else
                tcV(VialIndex) = NaN;
           end
           
           ncV(VialIndex) = mean(nV(:,VialIndex));
end
%Make Output Var Of Centroids - Append to file
eval(strcat('centroid',strOutputTag,'= [tcV; ncV]'));
save('ActivityCentroids.mat',strcat('centroid',strOutputTag),'-append')

ylimitsNTracklets = ceil(max(nV(:))/100)*100;
ylimitsmu = ceil(max(muV(:)));

%% Plot All Vials Independently %%%
        strtitle = sprintf('Mean Activity Of Each Vial-Sliding Window %d hours',TimeFrameWidth/3600);
        t = length(ExpTrackResultsInTime);
        Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;
        hf = figure('Name',strcat(ExpCondFood{ConditionIndex},strtitle));

        % Get the initial set of default plot colors.
        initialColorOrder = get(gca,'ColorOrder') % Initial
        
        % Apply the new default colors to the current axes.
        set(gca, 'ColorOrder', plotcoloursPerVial, 'NextPlot', 'replacechildren');
        
        subplot(3,1,1)
        set(gca, 'ColorOrder', plotcoloursPerVial, 'NextPlot', 'replacechildren');
        plot(Exptime,muV(:,:));
        title('Mean Speed in px/sec');
        ylim([0 ylimitsmu]);

        hh= subplot(3,1,2);
        %get(hh,'position')
         % Apply the new default colors to the current axes.
        set(gca, 'ColorOrder', plotcoloursPerVial, 'NextPlot', 'replacechildren');
        plot(Exptime,stddV(:,:));
        title('STD Dev ');
        
        
        
        subplot(3,1,3)
        title('Number of samples');
        % Apply the new default colors to the current axes.
                
        hold on;
   

        set(gca, 'ColorOrder', plotcoloursPerVial, 'NextPlot', 'replacechildren');
        plot(Exptime,nV,'LineWidth',1.6);
        
        hold on;
        %Plot Centroids
        for (VialIndex=1:1:maxVialCount)
            plot(tcV( VialIndex) ,ncV(VialIndex),'.','markers',22,'MarkerEdgeColor',plotcoloursPerVial(VialIndex,:)) ;
            pltC = plot(  tcV( VialIndex) ,ncV(VialIndex),'o','markers',22,'MarkerEdgeColor','k');
        end
      

        xlabel('Hour');
        ylim([0 ylimitsNTracklets]);
        hold off;
        
        %Make Legend
        strLegend = '';
        j = 0;
        for k=1:maxVialCount
            [cond,vial] =find(VialPairsPerCondition==k);
            
            j=j+1; 
            strLegend{j} = strcat(ExpCondTitles{cond},'-V',num2str(vial));
        end
        legend( strLegend,'Location','southoutside','Orientation','vertical','Position',[0.84 0.45 0.124 0.43])
        %set(hh,'position',[0.13 0.2 0.77 0.12]); %Fix Last plot after adding legends
        saveas(hf,strcat('figures/meanALLVialIndy-',strOutputTag,'SpeedSlidingWindow',ExpCondFood{cond},'.png'));
%    end


%% Do Mean Speed Per Condition Per Time Window
meanConditionSpeeds  = {};
n  = zeros(length(ExpTrackResultsInTime),9);
mu  = zeros(length(ExpTrackResultsInTime),9);
stdd  = zeros(length(ExpTrackResultsInTime),9);
ConditionIndex = 1;


for (ConditionIndex=1:(ConditionIndexMax))
    for t=1:length(ExpTrackResultsInTime)
            ExpTrackResults         = ExpTrackResultsInTime{t};
            %TODO: Add Filters Here
            %ExpTrackResults =             
            
            ResSet                  = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});

            if isempty(ResSet)
                continue;
            end
            meanConditionSpeeds{ConditionIndex}   = vertcat(ResSet.MeanSpeed);
            n(t,ConditionIndex)                   = length(meanConditionSpeeds{ConditionIndex});
            mu(t,ConditionIndex)                  = mean(meanConditionSpeeds{ConditionIndex});
            stdd(t,ConditionIndex)                = std(meanConditionSpeeds{ConditionIndex});
    end
    
         %Calc Central Moment of Inertia
         tc(ConditionIndex) = sum((1:t)'.*n(:,ConditionIndex))/sum( n(:,ConditionIndex));
         nc(ConditionIndex) =   mean(n(:,ConditionIndex));
    
end



%% Plot Vials Grouped per condition %%
ylimitsTracklets = ceil(max(n(:))/100)*100;
ylimits = ceil(max(mu(:)));

strtitle = sprintf('Mean Activity -Sliding Window %d hours',TimeFrameWidth/3600);
% Plot Results - In Groups or - In sets of 3-genotypes For each Food Condition
CondGrouping = ConditionIndexMax; %When =ConditionIndexMax then plot them all together

for i=1:length(ConditionGroups) %plot Per Condition Groups 
    CondIndexes = ConditionGroups{i};
    
%    for (ConditionIndex=1:CondGrouping:ConditionIndexMax)

        t = length(ExpTrackResultsInTime);
        Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;
        hf = figure('Name',strcat(ExpCondFood{CondIndexes(1)},strtitle));

        subplot(3,1,1)
        %plot(Exptime,mu(:,ConditionIndex:(ConditionIndex+CondGrouping-1))   );
        plot(Exptime,mu(:,CondIndexes)   );
        title('Mean Speed in px/sec');
        ylim([0 ylimits]);

        subplot(3,1,2)
        plot(Exptime,n(:,CondIndexes));


        hold on;
        title('Number of samples');
        %Plot Centroids
        plot(Exptime(max(round(tc),1)),nc,'x')
        
        hold off;
        
        hh= subplot(3,1,3)
        %get(hh,'position')

        plot(Exptime,stdd(:,CondIndexes));
        title('STD Dev ');
        xlabel('Hour');
        ylim([0 ylimits]);

        %Make Legend
        strLegend = '';
        j = 0;
        for k=1:length(CondIndexes)
            j=j+1; 
            strLegend{j} = strcat(ExpCondFood{CondIndexes(k)},ExpCondTitles{CondIndexes(k)});
        end
        legend( strLegend,'Location','southoutside','Orientation','vertical','Position',[0.84 0.01 0.124 0.43])
        %set(hh,'position',[0.13 0.2 0.77 0.12]); %Fix Last plot after adding legends
        saveas(hf,strcat('figures/meanVial',strOutputTag,'SpeedSlidingWindow',ExpCondFood{CondIndexes(1)},'-G',num2str(i),'.png'));
%    end
end %end Condition Groups
%,strcat(ExpCondFood{ConditionIndex+1},ExpCondTitles{ConditionIndex+1}),strcat(ExpCondFood{ConditionIndex+2},ExpCondTitles{ConditionIndex+2})


%% PLOT ALL GENOTYPES %%%
for (ConditionIndex=1:CondGrouping:ConditionIndexMax)
        CondIndexes = ConditionIndex:(ConditionIndex+CondGrouping-1);
        t = length(ExpTrackResultsInTime);
        Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;
        hf = figure('Name',strcat(ExpCondFood{ConditionIndex},strtitle));

        subplot(3,1,1)
        plot(Exptime,mu(:,ConditionIndex:(ConditionIndex+CondGrouping-1))   );
        title('Mean Speed in px/sec');
        ylim([0 ylimits]);
                
        hh= subplot(3,1,2)
        %get(hh,'position')
        set(gca, 'ColorOrder', plotcoloursPerCondition, 'NextPlot', 'replacechildren');
        plot(Exptime,stdd(:,CondIndexes));
        title('STD Dev ');
        xlabel('Hour');
        ylim([0 ylimits]);

        
        subplot(3,1,3)
        set(gca, 'ColorOrder', plotcoloursPerCondition, 'NextPlot', 'replacechildren');
        plot(Exptime,n(:,CondIndexes),'LineWidth',2.1);
        title('Number of samples');
        hold on;
        %Plot Centroids
        for (k=1:1:ConditionIndexMax)
            plot(Exptime(max(round(tc(k)),1)),nc(k),'.','markers',22,'MarkerEdgeColor',plotcoloursPerCondition(k,:));
            pltC = plot(Exptime(max(round( tc(k) ),1) ),nc(k),'o','markers',22,'MarkerEdgeColor','k');
        end
      
        
        %Make Legend
        strLegend = '';
        j = 0;
        for k=1:length(CondIndexes)
            j=j+1; 
            strLegend{j} = strcat(ExpCondFood{CondIndexes(k)},ExpCondTitles{CondIndexes(k)});
        end
        legend( strLegend,'Location','southoutside','Orientation','vertical','Position',[0.84 0.01 0.124 0.43])
        %set(hh,'position',[0.13 0.2 0.77 0.12]); %Fix Last plot after adding legends
        saveas(hf,strcat('figures/meanALLVial',strOutputTag,'SpeedSlidingWindow',ExpCondFood{CondIndexes(1)},'.png'));
%    end
end %end Condition Groups

%% Plot Histogram Of Speed Within A chosen time Window
goToHour = 110;
t= round((goToHour*3600 - VialAge(1))/timeAdvance);
ExpTrackResults = ExpTrackResultsInTime{t};

strtitle = sprintf('Speed Histogram -@t:%d for %d hours',goToHour,TimeFrameWidth/3600);

for i=1:length(ConditionGroups) %plot Per Condition Groups 
    CondIndexes = ConditionGroups{i};

    hf = figure('Name',strcat(ExpCondFood{CondIndexes(1)},strtitle));
    for (k=1:length(CondIndexes))
        ConditionIndex = CondIndexes(k);
        t = length(ExpTrackResultsInTime);
        Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;

        subplot(length(CondIndexes),1,k);
        ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
        
        if (length(ResSet) == 0) 
            continue;
        end
        meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
        [cnt,bin]                            = hist(meanConditionSpeeds{ConditionIndex},nbins);
        hist(meanConditionSpeeds{ConditionIndex},nbins);
        n       = length(meanConditionSpeeds{ConditionIndex});
        mu      = mean(meanConditionSpeeds{ConditionIndex});
        stdd    = std(meanConditionSpeeds{ConditionIndex});
        strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );
        title(strTitle);

        %Make Legend
        strLegend = '';
        j = 0;
        for k=1:length(CondIndexes)
            j=j+1; 
            strLegend{j} = strcat(ExpCondFood{CondIndexes(k)},ExpCondTitles{CondIndexes(k)});
        end
        %legend( strLegend,'Location','southoutside','Orientation','vertical','Position',[0.84 0.01 0.124 0.43])
        %set(hh,'position',[0.13 0.2 0.77 0.12]); %Fix Last plot after adding legends
        saveas(hf,strcat('figures/SpeedHistTracklet',strOutputTag,'tHour',num2str(goToHour),'_',ExpCondFood{CondIndexes(1)},'-G',num2str(i),'.png'));
    end
end %end Condition Groups




%% Plot Histogram Of Speed Within A chosen time Window
t= round((goToHour*3600 - VialAge(1))/timeAdvance);
ExpTrackResults = ExpTrackResultsInTime{t};

hold off;
ConditionIndex = 1;


hf = figure('Name',strcat(ExpCondFood{ConditionIndex},'MEAN SPEED'));

subplot(3,1,1);
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
if length(ResSet) > 0 
    meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
    [cnt,bin]                            = hist(meanConditionSpeeds{ConditionIndex},nbins);
    hist(meanConditionSpeeds{ConditionIndex},nbins);
    n       = length(meanConditionSpeeds{ConditionIndex});
    mu      = mean(meanConditionSpeeds{ConditionIndex});
    stdd    = std(meanConditionSpeeds{ConditionIndex});
    strTitle = sprintf('%s %s mean: %0.3f std:%0.3f n:%d',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex},mu,stdd, n );
    title(strTitle);
end

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


%% Do Scatter Plots Of Speeds
hold off;

hf = figure('Name',strcat(ExpCondFood{ConditionIndex},' Tracklet mean speeds across time'));
startT               = round((goToHour*3600 - VialAge(1))/timeAdvance);

for ConditionIndex=1:ConditionIndexMax
    display(ConditionIndex);
    
    cnt = 0;
    speeds = [];
    ttime = [];
    for t=1:length(ExpTrackResultsInTime)
        cnt = cnt+1;
        ExpTrackResults = ExpTrackResultsInTime{t};
        ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex,: )});
        if isempty(ResSet) 
            continue;
        end
        meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
        Exptime = (VialAge(1)+(1:t)*timeAdvance)/3600;

        mxy = max(meanConditionSpeeds{ConditionIndex});
        if ylimits < mxy 
            ylimits = mxy
        end

        n       = length(meanConditionSpeeds{ConditionIndex});
        %mu      = mean(meanConditionSpeeds{ConditionIndex});
        %stdd    = std(meanConditionSpeeds{ConditionIndex});
        %title(strTitle);
        for i=1:length(meanConditionSpeeds{ConditionIndex})
            speeds(cnt) = meanConditionSpeeds{ConditionIndex}(i);
            ttime(cnt)  = (t*timeAdvance + VialAge(1))/3600;
        end

    end

    hold off;
    subplot(ConditionIndexMax,1,ConditionIndex)
    strTitle = sprintf('%s %s ',ExpCondFood{ConditionIndex},ExpCondTitles{ConditionIndex} );
    scatter(ttime,speeds,'.');
    ylim([0 ylimits]);
    xlim([80 (maxRecordingTime + VialAge(1))/3600]);

    title(strTitle);

end
saveas(hf,sprintf('figures/TrackletMeanSpeedScatter-%s-%dHour.png',strOutputTag,goToHour))

