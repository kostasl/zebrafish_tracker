%% Import CSV Files

%%FOLDER NAMES SHOULD BE OF THE FORMAT: EXP_6-7_20151123_5sec
%5sec is the timelapse period in seconds

scriptPath = which('processFilesTracks.m');
%frameN,TrackID,TrackBlobLabel,Centroid_X,Centroid_Y,Lifetime,Active,Inactive
[framePeriod,VialAge,ExpIDs,ExpTrack ] = importCSVtoCell( '*V*_tracks','EXP*' );
save('LarvaTrackData.mat','ExpTrack');


%% Plot Results for Active Larva Per Vial
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
%%
%Holds Track id That fit the spec of lifetime
FilteredTrackIDs = {}; 
%Filter Each Experiments Data set
MinLifetime = 50;
for (e=1:size(ExpTrack,1))
    for (v=1:size(ExpTrack,2))
        %Check If Empty cell - Data for Vial N/A
       if isempty(ExpTrack{e,v}) 
            continue;
       end
       % Get Track Ids that are longer than MinLifetime
       FilteredTrackIDs{e,v} = [unique(ExpTrack{e,v}(find(ExpTrack{e,v}(:,6)>MinLifetime),2))];
       
       %Go through Each ID - Get Mean Speed
       for i = 1:length(FilteredTrackIDs{e,v})
           trkID = FilteredTrackIDs{e,v}(i);
           
           % Find Positions /Sorted By Frame Number / Get distance
           % travelled
           
       end
           
    end
end


%Combine the vial Indexes
for (ConditionIndex=1:9)
    for (vi = 1:length(VialPairsPerCondition)) 

        %Go through Each Experiments Data set
        for (e=1:size(ExpTrack,1))
            cl = VialAge(e); %Col -timepoint - Srtart from Vial Age
            rescolIndex = (vi-1)*size(ExpTrack,1) + e;
            %display(rescolIndex);
            %Go through each data point and Place it in the correct bin in the timeData vector

            %Check If Experiment Had this Vial Number - Replication 
            if (VialPairsPerCondition(vi) <= size(ExpTrack,2))
                %%Collect the data from Vial to The matrix of results - Limit
                %%To the timeframe set by timepoints
                for (i=1:length(ExpTrack{e, VialPairsPerCondition(vi)}))

                    %Check If Data longer than required Time Vector
                    if (cl > timePoints)
                        break;
                    end
                    %Fill The gap between each sample with the same value
                    dataMatrix(rescolIndex,cl:(cl+framePeriod(e))) = ExpTrack{e, VialPairsPerCondition(vi)}(i,3);
                    cl = cl + framePeriod(e); %INCREMENT TO NEXT REAL TIME SAMPLED
                end
            end
            %Filter
            %dataMatrix(e,:) = medfilt1(dataMatrix(e,:),2*40*framePeriod(e));
            %display(strcat('Median Filter applied :',num2str(40*framePeriod(e))) );
        end
    end
end

%vialtrackerResults{ConditionIndex} = struct('ActiveCount',datNV,'MeanActiveCount',meanNLarva,'MeanErrCount',stdNLarva);
