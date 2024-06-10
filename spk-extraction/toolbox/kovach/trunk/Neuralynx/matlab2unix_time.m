function unix_time = matlab2unix_time(matlab_time)

% Convert matlab time (days since 1/1/1) to unix time (seconds since 1/1/9170) 

unix0 = datenum('1/1/1970'); %Beginning of time according to unix


if ischar(matlab_time)
    matlab_time = datenum(matlab_time);
end

unix_time = (matlab_time-unix0)*3600*24;