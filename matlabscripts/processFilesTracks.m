%Import Track Data into A form that retains Experiment/vial structure -
% Process data giving Euclidean distance in Px / normalize by framerate
% Note: Lifetimes In track data are not reliable!

%% Import CSV Files
%%FOLDER NAMES SHOULD BE OF THE FORMAT: EXP_6-7_20151123_5sec
%5sec is the timelapse period in seconds

%Put the script dir in the path
addpath(fileparts(which('processFilesTracks.m')))
%Change dir to where the data files are
%frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive
%cd /home/klagogia/Videos/LarvaTrackPilot/DataOut %Office
cd /media/kostasl/FlashDrive/PilotVialTrack/DataOut %Home
%%Import FROM CSV FILES
%VialAge : Age of the vials for an experiment j - from embryo to the beginning of timelapse Recording
[framePeriod,VialAge,ExpIDs,ExpTrack ] = importCSVtoCell( '*V*_tracks','EXP_A*' );



%Transform - Y Inversion
%ExpTrack{:,:}(:,5) = 768 - ExpTrack{:,:}(:,5)



%% Organize and Process Imported data
%Put the script dir in the path

%Give 3 days data points 1 sec each.
% Genotypes are 3 organized in this order : 1st WT (oregonR), 2nd Genetic Control, 3rd AlfaBeta Mutant
ConditionIndex      = 1; %Experimental Condition ID : Food(Condition)/Genetype Combinations
ConditionIndexMax   = 3; %Defines max exp. configuration being replicated - ex. 1= Food1/Gen1 1= Food2/Gen1. Combos

% The videos have 2 rows of 9 vials - Vials 1-10 have identical conditions so they go in PAIRS
%VialPairsPerCondition = [[1,10];[2,11];[3,12];[4,13];[5,14];[6,15];[7,16];[8,17];[9,18]]; %OR Normal Food


%For new 2016/03 Setup We have 1 row - 3 conditions - 3 reps Each
VialPairsPerCondition = [[1,2,3];[4,5,6];[7,8,9];]; %OR Normal Food
timePoints = max(VialAge) + 24*3*3600;%Total Time points in seconds over which to analyse data
%FramePeriod sampled at each timelapse Experiment -

%%THESE DO NOT CORRESPOND TO TIMES OF VIDEOS
%framePeriod = [20;5;5;2;2;2;2;2];
display('The loaded videos had frame Periods designated in folder name:')
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

ExpTrackResultsInTime = {};
bVerbose=0;
%Filter Each Experiments Data set

MinLifetime     = 2; %Minimum Number of Path Steps
MaxLifetime     = 20000; %Maximum Number of Path Steps
MinDistance     = 5; %Minimum Track length to consider def 10
MinStepLength   = 1; %%Cut Tracklet when 2-frame displacement drops below value 
MaxStepLength   = 55; %MaxpxSpeed -->Between two frames rejects steps larger than this
TimeFrameWidth  = 1*3600; %Frame Sliding Window in sec Overwhich results are averaged

% Organize data in a Sliding Window
InitTime = 0*3600; %Start processing Data from InitTime / Default 0
wi = 0;
%maxRecordingTime = 3*24*3600; %3 Days
%Estimate Max FrameN from 1st Experiment
e = 1;
maxRecordingTime = max(vertcat([ExpTrack{e,1}(:,1)]))*framePeriod(e);
timeAdvance = 5*60; %Fwd Time Step in secs


for StartTime=(InitTime + TimeFrameWidth):timeAdvance:(maxRecordingTime)
   wi = wi+1;
   %                                                    (ExpTrack,ExpIDs,framePeriod,MinLifetime, MaxLifetime, MinDistance, MaxpxSpeed, FromTime,TimeWindow, MinpxSpeed,bVerb ) 
   ExpTrackResultsInTime{wi} = ExtractFilteredTrackData(ExpTrack,ExpIDs,framePeriod,MinLifetime, MaxLifetime, MinDistance, MaxStepLength, StartTime,TimeFrameWidth, MinStepLength ,bVerbose);
   disp(StartTime/maxRecordingTime);%%Sho Fraction of Calculation Completed
end

save('LarvaTrackData_B.mat');
%% Plot Indicative results - Distribution of mean Tracklet Speeds
plotMeanSpeed;

%% Plot Track Length
plotTrackLengthDistributions;



%% Plot Example Tracks
goToHour =85; %Exp Hour - with 0 being Embryo Placement

t= (goToHour*3600 - VialAge(1))/timeAdvance; %
ExpTrackResults = ExpTrackResultsInTime{t};

% NOTE: Y values are inverted since 0 Point in plot as at the bottom
imgSize = [1024,768];
colour = {'r','m','k','c','b','g','k','--r','--m','--k','--c','--b','--g','--k','-r','-k','-g'};
hf = figure('Name',sprintf('Tracks at t=%d-%d hours',goToHour,goToHour+TimeFrameWidth/3600));
xlim([0 imgSize(1)]);
ylim([0 imgSize(2)]);

hold on;

e = 1;
MinLifetime = 0;
for e=1:size(ExpTrackResults,1)
    for (v=1:size(ExpTrackResults,2))
        
        %FilteredTrackIDs = unique(ExpTrackResults{e,v}( find(ExpTrackResults{e,v}(:,6)>MinLifetime),2 ));
        %find(vertcat(ExpTrackResults{e,v}.PointCount) > 6)
        %trackN = length(FilteredTrackIDs);
        trackN = length(ExpTrackResults{e,v});
        for (i=1:1:trackN)
            %trkID = FilteredTrackIDs(i);
            trackData = ExpTrackResults{e,v}(i).Positions;
            trackData(:,3) = imgSize(2) - trackData(:,3); 
            mark = colour{randi(17)};
            plot(trackData(:,2),trackData(:,3),mark);
            scatter(trackData(1,2),trackData(1,3),'.'); %Start
            l = size(trackData,1); %Count Records

            %lRec(i) = l;
            scatter(trackData(l,2),trackData(l,3),'X') %Finsih
        end
    end
end
title('Plot Sample Track');
saveas(hf,sprintf('figures/VialTracklets_t%d-%dHours.png',goToHour,goToHour+TimeFrameWidth/3600));




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
