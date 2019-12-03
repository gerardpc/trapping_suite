%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
%   time = time array of the timetrace
%   x    = a matrix dim(N_Ch,tt_time*SampRate) where each row contains
%           the timetrace of each axis. Marc$Fran x = [z;x;y;x_out-of-loop]
%   
%%%%%%%%% FUNCTION INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  iname   = file name including extension [string]
%  ipath   = File path [string]
%  tt_time = Duration of the time trace [s]
%  N_Ch    = Number of chanels to read  
%  SampRate       = Sampling rate [samples/s]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  example 4CH : [time,x] = extract_tt_NCh(iname,ipath,40,4,625000)

function [time,x_mode] = extract_trace(file_name, file_path, rec_time)

ifile   =   fullfile(file_path, file_name);
i_fid   =   fopen(ifile,'r');

SampRate = 625000;
N_Ch = 4;

t_sample = 1/SampRate;
N_tt = rec_time*SampRate;

%%% Reading all data points and puting them in an array such as tt = [x1,y1,z1,...,n_ch1,x2,y2,z2,...,n_ch2,...]   
n_To_Throw = 150000;
Data_To_Throw  = fread(i_fid,N_Ch*n_To_Throw,'int16',0,'ieee-le.l64');    % Reads single string containing all the dim.   

tt_all = fread(i_fid,'int16',0,'ieee-le.l64');
tt_all = tt_all';

%%% Separating channels

x = zeros(N_Ch,N_tt);

for j=1:N_Ch
   x(j,:) = tt_all(j:N_Ch:N_tt*N_Ch);
end 

%%% Time array

time = (t_sample:t_sample:rec_time);
x_mode = x;
fclose(i_fid);
end

