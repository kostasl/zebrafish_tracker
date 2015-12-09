function [ ResultsCell ] = importCSVtoCell(filePatt,dirPatt )
%IMPORTCSVTOCELL Imports experiment results from CSV to CELLArray
%   scans all sub dirs that start with Exp, use filePatt like '*V*_N'
dataFileCount = 0;

files = dir(fullfile(dirPatt));
Dirs = files(find(vertcat(files.isdir)));

for d=1:length(Dirs)

    importDir = Dirs(d);
    %if (isempty(strfind(importDir.name, 'Exp') ))
%        continue;
%    else
        display(strcat('Importing from Dir:',importDir.name));
%    end
        
    files = dir(fullfile(importDir.name,strcat(filePatt,'.csv')));
    if (exist('ResultsCell') == 0)
        ResultsCell = cell(size(Dirs,1)-2,size(files,1) );
    end
    
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


