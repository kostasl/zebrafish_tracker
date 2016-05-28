function [ frameRates,ExpVialAge,ExpID ,ResultsCell] = importCSVtoCell(filePatt,dirPatt )
%IMPORTCSVTOCELL Imports experiment results from CSV to CELLArray
%   scans all sub dirs that start wit h Exp, use filePatt like '*V*_N'
dataFileCount = 0;

files = dir(fullfile(dirPatt));
Dirs = files(find(vertcat(files.isdir)));

assert(length(Dirs) > 0,'No source directories found');
ExpID = {'','','','',''};

for d=1:length(Dirs)

    importDir = Dirs(d);
    %if (isempty(strfind(importDir.name, 'Exp') ))
%        continue;
%    else
     display(strcat('Importing from Dir:',importDir.name));
     if exist(fullfile(importDir.name,'timing.csv'),'file') == 0 
         error('Experiment folder is missing timing.csv file. This file needs two data rows 1st the embryo collected date-time,2nd and video recording date-time');
     end
     ExpID{d} = importDir.name; %Save Directory Name To Identify Experiment
%    end
    %%Get Time Data stored in File timing.csv
    
    timings = csvread(fullfile(importDir.name,'timing.csv'),1);
    ExpVialAge(d) = etime(timings(2,:),timings(1,:));
    sprintf('VialAge in days %0.2f',((ExpVialAge(d)/3600)/24.0))
    files = dir(fullfile(importDir.name,strcat(filePatt,'.csv')));
    if (exist('ResultsCell') == 0)
        ResultsCell = cell(size(Dirs,1)-2,size(files,1) );
    end
    
    strtimelapsePeriod= regexp(importDir.name,'[-+]?([0-9]*\.[0-9]+|[0-9]+)sec','match');
    frameRates(d) = str2double(regexp(strtimelapsePeriod{1},'[-+]?([0-9]*\.[0-9]+|[0-9]+)','match'));
   
    for i=1:length(files);
        %Get Vial Number from File Name
        VNo = regexp(files(i).name,'V\d+','match');
        VNo = str2double(regexp(VNo{1},'\d+','match'));
        display(strcat(fullfile(importDir.name,files(i).name),' -Vial Number:',num2str(VNo) ));
       
        %Exp{d,i} = ;
        %Put Result on Index representing 
        ResultsCell{d,VNo} = csvread(fullfile(importDir.name,files(i).name),1);
        display(strcat('* CSV imported :',num2str(size(ResultsCell{d,VNo},1)),  ' lines') );
        dataFileCount = dataFileCount + 1;
    end

    display(strcat('    Total Imported :',num2str(i), ' *******' ));
end

display(strcat('    CSV total Imported:',num2str(dataFileCount)) );


end


