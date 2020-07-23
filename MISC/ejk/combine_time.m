function [ time ] = combine_time( func,time1,time2 )
%either add or subtract times
%input:
%func: @plus or @minus
%time1/time2: [hour,minute,second] in 24 hour time
%if there is subtraction, it would be time2-time1
%output: 1x3 vector same as time1 and time2
time1_in_seconds = time1(1)*3600+time1(2)*60+time1(3);
time2_in_seconds = time2(1)*3600+time2(2)*60+time2(3);

time_in_seconds = feval(func,time2_in_seconds,time1_in_seconds);


if time_in_seconds > 0 && time_in_seconds < 24*60*60
    time(1) = floor(time_in_seconds/3600);
    time(2) = floor(mod(time_in_seconds/60,60));
    time(3) = mod(time_in_seconds-(time(1)*3600),60);
elseif time_in_seconds == 0 || time_in_seconds == (24*60*60)
    time = [0,0,0];
elseif time_in_seconds > 24*60*60
    time_in_seconds = time_in_seconds - (24*60*60);
    time(1) = floor(time_in_seconds/3600);
    time(2) = floor(mod(time_in_seconds/60,60));
    time(3) = mod(time_in_seconds-(time(1)*3600),60);
else
    time_in_seconds = 24*60*60 + time_in_seconds;
    time(1) = floor(time_in_seconds/3600);
    time(2) = floor(mod(time_in_seconds/60,60));
    time(3) = mod(time_in_seconds-(time(1)*3600),60);

end
end

