%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert analog number to red pitaya resolution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dig_out] = rp_resolution(an_in)
% redpitaya resolution
res = 2^13;
dig_out = an_in; % preallocate

%% round
%round up
up = ceil(an_in*res)/res;
%round down
down = floor(an_in*res)/res;
avg_round = mean([up, down]);

for i = 1:length(an_in)
    if(avg_round(i) > an_in(i))
        dig_out(i) = down(i);
    else 
        dig_out(i) = up(i);
    end
end
end