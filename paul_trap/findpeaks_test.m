v_length = 500;
x = linspace(0, 10, v_length);
y = sin(x) + normrnd(0, 0.01, 1, v_length);
%plot(x,y);
min_dist = 2;
y = smooth(y);
[pks, locs] = findpeaks(y, x, 'MinPeakDistance',min_dist, 'SortStr','descend', 'Npeaks',2);
max(locs);