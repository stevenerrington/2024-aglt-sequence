function [mltime,mltime_sec,ml_datestr] = unix2matlab_time(unix_time)

% Convert unix time (time since 1/1/9170) in micro- milli- or unit seconds 
% to matlab time format (days since 1/1/1

unix0 = datenum('1/1/1970'); %Beginning of time according to unix


if unix_time > 1000*365*24*3600 % Assume it is not seconds in this case
    
    unix_time = double(unix_time)/1000;
    
    if unix_time > 1000*365*24 % Must  be microseconds!
        
          unix_time = double(unix_time)/1000;
    end
end

mltime = unix_time/3600/24 + unix0;

if nargout > 1
    mltime_sec = mltime*24*3600;
end

if nargout > 2
    ml_datestr = datestr(mltime);
end