function [date_marker, time_marker, log_time] = get_clock()
c = clock;
date_marker = sprintf('%4.0f%02.0f%02.0f',c(1:3));
time_marker = sprintf('%02.0f%02.0f%02.0f',c(4:6));
log_time = sprintf('%02.0f:%02.0f:%02.0f',c(4:6));
return
