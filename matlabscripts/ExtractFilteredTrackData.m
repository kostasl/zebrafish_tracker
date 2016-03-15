function [ExpTrackResults] = ExtractFilteredTrackData(ExpTrack,ExpIDs,framePeriod,MinLifetime, MaxLifetime, MinDistance, MaxpxSpeed, FromTime,TimeWindow, MinpxSpeed,bVerb ) 
% Summary: Returns Cell Array of Structures holding each track selected vial the given filters
% Parameters :
% ExpTrack - The imported track Data from the CSV files organized in cell
% array of per experiment / vial
% ExpIDs Holds Names of Experiements 
% MinLifetime Minimum Number of Path Steps
%MaxLifetime Maximum Number of Path Steps
%MinDistance   Minimum Track length to consider
%MaxStepLength  Between two frames rejects steps larger than this
% FromTime : Filter Tracklets that are X sec after beginning of video
%TimeWindow  : Set Time after FromTime to extract data points from
FilteredTracks = {}; 
%MinpxSpeed = 2;
% bVerb = 0; Flag to Set if function should be verbose

for (e=1:size(ExpTrack,1))
    if (bVerb)
        disp(char(ExpIDs(e)));
    end
    
    for (v=1:size(ExpTrack,2))
        %Check If Empty cell - Data for Vial N/A
       if isempty(ExpTrack{e,v}) 
            continue;
       end
       %Lifetimes In track data are not reliable!
       % Get Track Ids that are longer than MinLifetime
       FiltIndexes = find( ExpTrack{e,v}(:,6)>MinLifetime & ...
                           ExpTrack{e,v}(:,6)  < MaxLifetime & ...           
                           ExpTrack{e,v}(:,1)*framePeriod(e) > FromTime & ...
                           ExpTrack{e,v}(:,1)*framePeriod(e) < (FromTime + TimeWindow) );
       
       FilteredTrackIDs = unique(ExpTrack{e,v}(FiltIndexes, 2));
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
                %For some Reason Tracker Reuses  trackIDs
                %error('Found Broken Lifetime');
           end
           %Check Track Lifetime Again - Filter If Less than Required Size
           if (length(trackData(:,2)) < MinLifetime || length(trackData(:,2)) < 2)
               continue; %Go to Next
           end
          
           %Normalize By FrameRate
           pathSteps          = sqrt(diff(trackData(:,2)).^2+diff(trackData(:,3)).^2);
           
           
           %Find 1st point When Track Stops or goes too fast And Truncate
           datbreakpoint  = find(pathSteps(:,1) < MinpxSpeed | pathSteps(:,1) > MaxpxSpeed );
           %datbreakpoint = join(datbreakpoint, );
           
            if (length(datbreakpoint) > 1)
                pathSteps = pathSteps(1:datbreakpoint(1),:);
            end
            
            stepsCount = length(pathSteps);
            trackDist = sum(pathSteps);
            %Apply Step Count and MinDistance Filter 
           if ((stepsCount < MinLifetime) || (stepsCount > MaxLifetime) || (trackDist < MinDistance))
               continue; %Go to Next
           end
           
           ii = ii + 1;     
            
           npathSteps(ii)     = stepsCount;
           pathdistance(ii)   = trackDist;
           meanspeed(ii)      = trackDist/(stepsCount*framePeriod(e));
           stdspeed(ii)       = std(pathSteps) / framePeriod(e);
           
           FilteredTracks{ii} = struct('TrackID',trkID,'PointCount',length(trackData),'Positions',trackData,'Length',pathdistance(ii),'MeanSpeed',meanspeed(ii),'StdDevSpeed', stdspeed(ii));
           
           assert(~( isnan(meanspeed(ii)) || isnan(stdspeed(ii))  ), ...
           sprintf('NaN Encountered in mean track speeds, %s  V:%d, TrackID: %d',char(ExpIDs(e)),v,trkID ));
       
       end %Each TrackID
       
       meanVialSpeed    = mean(meanspeed);
       meanVialStd      = mean(stdspeed);
       meanVialDistance = mean(pathdistance);
       meanTrackSteps   = mean(npathSteps);
       assert(~(isnan(meanVialSpeed) || isnan(meanVialStd)),'NaN Encountered in mean track speeds');
       
       if (bVerb)
            disp( sprintf('V:%d Mean Vial Speed: %0.3f std:%0.4f mean path Distance: %0.3f Steps:%d',v,meanVialSpeed,meanVialStd,meanVialDistance,meanTrackSteps));
       end
       
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

end
