%%Import Track Data into A form that retains Experiment/vial structure -
% Process data giving Euclidean distance in Px / normalize by framerate
% Note: Lifetimes In track data are not reliable!

%% Import CSV Files
%%FOLDER NAMES SHOULD BE OF THE FORMAT: EXP_6-7_20151123_5sec
%5sec is the timelapse period in seconds

scriptPath = which('processFilesTracks.m');
%frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive
[framePeriod,VialAge,ExpIDs,ExpTrack ] = importCSVtoCell( '*V*_tracks','EXP*' );
save('LarvaTrackData.mat','ExpTrack');


%% Organize and Process Imported data
%Give 3 days data points 1 sec each.
% Genotypes are 3 organized in this order : 1st WT (oregonR), 2nd Genetic Control, 3rd AlfaBeta Mutant
ConditionIndex = 1; %Experimental Condition ID : Food/Genetype Combinations
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

%%How to process The tracks - 
%Calc:
% * Sort/filter Tracks by lifetime 
% * mean Tracklet speed 
% * mean Tracklet In Time
%Holds Track id That fit the spec of lifetime
FilteredTracks = {}; 
ExpTrackResults = {};
%Filter Each Experiments Data set
MinLifetime = 5;
for (e=1:size(ExpTrack,1))
    display(char(ExpIDs(e)));
    for (v=1:size(ExpTrack,2))
        %Check If Empty cell - Data for Vial N/A
       if isempty(ExpTrack{e,v}) 
            continue;
       end
       %Lifetimes In track data are not reliable!
       % Get Track Ids that are longer than MinLifetime
       FilteredTrackIDs = unique(ExpTrack{e,v}( find(ExpTrack{e,v}(:,6)>MinLifetime) ,2 ));
       %Unfiltered Unique Ids
       FilteredTrackIDs = unique(ExpTrack{e,v}(:,2));
       
       %Go through Each ID - Get Mean Speed
       ii = 0;
       for i = 1:length(FilteredTrackIDs)
           
           trkID = FilteredTrackIDs(i);
           % Find Positions /Sorted By Frame Number / Get distance
           % travelled
           trackData = sort(ExpTrack{e,v}(find(ExpTrack{e,v}(:,2)==trkID),[1,4,5]),1);
           %Check Track Lifetime Again - Filter If Less than Required Size
           if (length(trackData) < MinLifetime)
               continue; %Go to Next
           end
           ii = ii + 1;
           %Calc distance moved in px at each frame -frameN is on col 1, X(col 2) Y (col 3) Take Diff in sqrt((Xn-Xn+1)^2+(Yn-Yn+1)^2) -
           %Normalize By FrameRate
           meanspeed(ii) = mean(sqrt(diff(trackData(:,2)).^2+diff(trackData(:,3)).^2)) / framePeriod(e);
           stdspeed(ii) = std(sqrt(diff(trackData(:,2)).^2+diff(trackData(:,3)).^2)) / framePeriod(e);
           FilteredTracks{ii} = struct('TrackID',trkID,'Length',length(trackData),'Positions',trackData,'MeanSpeed',meanspeed(ii),'StdDevSpeed', stdspeed(ii));
           
           assert(~( isnan(meanspeed(ii)) || isnan(stdspeed(ii))  ), ...
           sprintf('NaN Encountered in mean track speeds, %s  V:%d, TrackID: %d',char(ExpIDs(e)),v,trkID ));
       
       end %Each TrackID
       
       meanVialSpeed = mean(meanspeed);
       meanVialStd   = mean(stdspeed);
       assert(~(isnan(meanVialSpeed) || isnan(meanVialStd)),'NaN Encountered in mean track speeds');
       
       %%%,'VialTrackCount',length(FilteredTracks),'VialMeanSpeed',meanVialSpeed,'VialMeanStdDevSpeeds',meanVialStd
       ExpTrackResults{e,v} = struct('Tracks',FilteredTracks);
       ExpTrackResults{e,v} = vertcat(ExpTrackResults{e,v}.Tracks);
       %Empty Buffer Cell Array
       FilteredTracks = {}; 
    end %Each Vial

end %Each Experiment

clear FilteredTracks
clear meanspeed
clear stdspeed
ylimits = 500;
xlimits = 2;
nbins = 1000;
%%Plot Indicative results - Distribution of mean Tracklet Speeds
hf = figure('Name','Normal Food');
subplot(3,1,1);
ConditionIndex = 1;
ResSet                      = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hold off;
hist(meanConditionSpeeds{ConditionIndex},nbins);
title('OR Normal Food');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','blue');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);

subplot(3,1,2);
ConditionIndex = 2;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title('CT Normal Food');
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
title('AB Normal Food');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/NFTrackletSpeedHist.pdf')



hold off;
hf = figure('Name','0.5% DMSO');
subplot(3,1,1);
ConditionIndex = 4;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title('OR 0.5% DMSO ');
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
title('CT 0.5% DMSO ');
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
title('AB 0.5% DMSO ');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO05TrackletSpeedHist.pdf')


%DMSO 1%
hold off;
hf = figure('Name','1% DMSO');
subplot(3,1,1);
ConditionIndex = 7;
ResSet                               = vertcat(ExpTrackResults{:,VialPairsPerCondition(ConditionIndex )});
meanConditionSpeeds{ConditionIndex}  = vertcat(ResSet.MeanSpeed);
hist(meanConditionSpeeds{ConditionIndex},nbins);
title('OR 1% DMSO ');
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
title('CT 1% DMSO ');
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
title('AB 1% DMSO ');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
set(h,'EdgeColor','w');
ylim([0 ylimits]);
xlim([0 xlimits]);
saveas(hf,'figures/DMSO10TrackletSpeedHist.pdf')





%% 
% %Combine the vial Indexes
% for (ConditionIndex=1:9)
%     for (vi = 1:length(VialPairsPerCondition)) 
% 
%         %Go through Each Experiments Data set
%         for (e=1:size(ExpTrack,1))
%             cl = VialAge(e); %Col -timepoint - Srtart from Vial Age
%             rescolIndex = (vi-1)*size(ExpTrack,1) + e;
%             %display(rescolIndex);
%             %Go through each data point and Place it in the correct bin in the timeData vector
% 
%             %Check If Experiment Had this Vial Number - Replication 
%             if (VialPairsPerCondition(vi) <= size(ExpTrack,2))
%                 %%Collect the data from Vial to The matrix of results - Limit
%                 %%To the timeframe set by timepoints
%                 for (i=1:length(ExpTrack{e, VialPairsPerCondition(vi)}))
% 
%                     %Check If Data longer than required Time Vector
%                     if (cl > timePoints)
%                         break;
%                     end
%                     %Fill The gap between each sample with the same value
%                     dataMatrix(rescolIndex,cl:(cl+framePeriod(e))) = ExpTrack{e, VialPairsPerCondition(vi)}(i,3);
%                     cl = cl + framePeriod(e); %INCREMENT TO NEXT REAL TIME SAMPLED
%                 end
%             end
%             %Filter
%             %dataMatrix(e,:) = medfilt1(dataMatrix(e,:),2*40*framePeriod(e));
%             %display(strcat('Median Filter applied :',num2str(40*framePeriod(e))) );
%         end
%     end
% end

%vialtrackerResults{ConditionIndex} = struct('ActiveCount',datNV,'MeanActiveCount',meanNLarva,'MeanErrCount',stdNLarva);
