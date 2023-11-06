function [S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name)    
file_name='ieee14cdf.txt';
fid=fopen(file_name); %open data file

%% bus data
l_str=fgetl(fid);% read the first line of data
%Title=[string(l_str(2:9)), string(l_str(11:30)), string(l_str(32:37)), string(l_str(39:42)), string(l_str(44)), string(l_str(46:end))]; 
%S_Base=str2num(Title(3));
S_Base=str2num(l_str(32:37));
l_str=fgetl(fid);%read the second line of data
Bus_data=[];
No_of_Buses=0;
while ischar(l_str) %loop until the end of bus data
    l_str=fgetl(fid);
    if(strcmp(l_str(1:4),'-999')==1)
        break;
    end
    index=19; l_str_num=l_str(index:end);
   l_num= str2num(l_str_num);
    No_of_Buses=No_of_Buses+1;
    Bus_data=[Bus_data; [No_of_Buses  l_num]];
end
%% line data
l_str=fgetl(fid);
Line_data=[];No_of_Lines=0;
while ischar(l_str)
    l_str=fgetl(fid);
    if(strcmp(l_str(1:4),'-999')==1)
        break;
    end
    index=1; l_str_num=l_str(index:end);
    l_num= str2num(l_str_num);
    No_of_Lines=No_of_Lines+1;
    Line_data=[Line_data;[No_of_Lines  l_num]];
end
end

    
    
    
    

