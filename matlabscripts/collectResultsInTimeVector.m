function [ dataMatrix ] = collectResultsInTimeVector( ExpDataCell,vialIndexes,framePeriod,timePoints )
%UNTITLED Collects results per vial across conditions from the Imported ExpData cell to a  Matrix
%   The data are collected in a data/second vector that syncs the data in
%   time according to the given framerate ie. The period of the timelapse experiment  
% ex: framePeriod = [5;5;2;2;2];
% Can combine data across vial indexes given as a vector in vialIndex
% Performes MedianFiltering to reduce the gaps/noise 

dataMatrix = zeros([length(vialIndexes)*size(ExpDataCell,1),timePoints+max(framePeriod)]);


%Combine the vial Indexes
for (vi = 1:length(vialIndexes)) 

    %Go through Each Experiments Data set
    for (e=1:size(ExpDataCell,1))
        cl = 1; %Col -timepoint
        rescolIndex = (vi-1)*size(ExpDataCell,1) + e;
        %display(rescolIndex);
        %Go through each data point and Place it in the correct bin in the timeData vector
        
        %Check If Experiment Had this Vial Number - Replication 
        if (vialIndexes(vi) <= size(ExpDataCell,2))
            %%Collect the data from Vial to The matrix of results - Limit
            %%To the timeframe set by timepoints
            for (i=1:length(ExpDataCell{e, vialIndexes(vi)}))

                %Check If Data longer than required Time Vector
                if (cl > timePoints)
                    break;
                end
                %Fill The gap between each sample with the same value
                dataMatrix(rescolIndex,cl:(cl+framePeriod(e))) = ExpDataCell{e, vialIndexes(vi)}(i,3);
                cl = cl + framePeriod(e); %INCREMENT TO NEXT REAL TIME SAMPLED
            end
        end
        %Filter
        dataMatrix(e,:) = medfilt1(dataMatrix(e,:),2*40*framePeriod(e));
        display(strcat('Median Filter applied :',num2str(40*framePeriod(e))) );
    end
end

%return dataMatrix;
end

