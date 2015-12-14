%%Import Track Data into A form that retains Experiment/vial structure -
% Process data giving Euclidean distance in Px / normalize by framerate
% Note: Lifetimes In track data are not reliable!

%% Import CSV Files
%%FOLDER NAMES SHOULD BE OF THE FORMAT: EXP_6-7_20151123_5sec
%5sec is the timelapse period in seconds

%Put the script dir in the path
addpath(fileparts(which('processFilesTracks.m')))
%Change dir to where the data files are
%frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive
cd /home/klagogia/Videos/LarvaTrackPilot/DataOut
%[framePeriod,VialAge,ExpIDs,ExpTrack ] = importCSVtoCell( '*V*_tracks','EXP*' );
%Transform - Y Inversion
%ExpTrack{:,:}(:,5) = 768 - ExpTrack{:,:}(:,5)

%save('LarvaTrackData.mat','ExpTrack');


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
MinLifetime = 10;
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
       %FilteredTrackIDs = unique(ExpTrack{e,v}(:,2));
       
       %Go through Each ID - Get Mean Speed
        ii = 0;
        npathSteps   = zeros(1,1);
        pathdistance = zeros(1,1);
        meanspeed    = zeros(1,1);
        stdspeed     = zeros(1,1);
        
         for (i=1:length(FilteredTrackIDs))
            trkID = FilteredTrackIDs(i);           
           %trkID = FilteredTrackIDs(i);
           % Find Positions /Sorted By Frame Number / Get distance
           % travelled
           trackData = ExpTrack{e,v}(find(ExpTrack{e,v}(:,2)==trkID),[1,4,5]);
           
           %Calc distance moved in px at each frame -frameN is on col 1, X(col 2) Y (col 3) Take Diff in sqrt((Xn-Xn+1)^2+(Yn-Yn+1)^2) -
           
           %Fix TrackData remove Spurious Jumps in frame
           datbreakpoint  = find(diff(trackData(:,1))>1);
           if ( isempty(datbreakpoint) == 0)
                trackData = trackData(1:datbreakpoint(1),:);
           end
           %Check Track Lifetime Again - Filter If Less than Required Size
           if (length(trackData) < MinLifetime)
               continue; %Go to Next
           end
          
           %Normalize By FrameRate
           pathSteps          = sqrt(diff(trackData(:,2)).^2+diff(trackData(:,3)).^2);
           
           %Check When Track Stops And Truncate
           datbreakpoint  = find(pathSteps(:,1)<2);
            if (length(datbreakpoint) > 1)
                pathSteps = pathSteps(1:datbreakpoint(1),:);
            end
           if (length(pathSteps) < MinLifetime)
               continue; %Go to Next
           end
           ii = ii + 1;     
            
           npathSteps(ii)     = length(pathSteps);
           pathdistance(ii)   = sum(pathSteps);
           meanspeed(ii)      = pathdistance(ii)/(npathSteps(ii)*framePeriod(e));
           stdspeed(ii)     = std(pathSteps) / framePeriod(e);
           
           FilteredTracks{ii} = struct('TrackID',trkID,'PointCount',length(trackData),'Positions',trackData,'Length',pathdistance(ii),'MeanSpeed',meanspeed(ii),'StdDevSpeed', stdspeed(ii));
           
           assert(~( isnan(meanspeed(ii)) || isnan(stdspeed(ii))  ), ...
           sprintf('NaN Encountered in mean track speeds, %s  V:%d, TrackID: %d',char(ExpIDs(e)),v,trkID ));
       
       end %Each TrackID
       
       meanVialSpeed    = mean(meanspeed);
       meanVialStd      = mean(stdspeed);
       meanVialDistance = mean(pathdistance);
       meanTrackSteps   = mean(npathSteps);
       assert(~(isnan(meanVialSpeed) || isnan(meanVialStd)),'NaN Encountered in mean track speeds');
       
       sprintf('Mean Vial Speed: %0.3f std:%0.4f mean path Distance: %0.3f Steps:%d',meanVialSpeed,meanVialStd,meanVialDistance,meanTrackSteps)
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

%% Plot Indicative results - Distribution of mean Tracklet Speeds
plotMeanSpeed;
%% Plot Track Length
plotTrackLengthDistributions;



%% Plot Example Tracks
colour = ['r','m','k','c','b','g','k'];
hf = figure('Name','Tracks');
hold on;

e = 1;
MinLifetime = 0;
for (v=1:1)
    FilteredTrackIDs = unique(ExpTrack{e,v}( find(ExpTrack{e,v}(:,6)>MinLifetime),2 ));
    
    for (i=1:length(FilteredTrackIDs))
        trkID = FilteredTrackIDs(i);
        trackData = ExpTrack{e,v}(find(ExpTrack{e,v}(:,2)==trkID),[1,4,5]);
        
        plot(trackData(:,2),trackData(:,3),'Color',colour(randi(7)));
        scatter(trackData(1,2),trackData(1,3),'x');
        l = length(trackData);
        lRec(i) = l;
        scatter(trackData(l,2),trackData(l,3),'.')
        
    end
end
title('Plot Sample Track');
xlim([0 1024]);
ylim([0 768]);
saveas(hf,'figures/VialTrackLets.png')




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
