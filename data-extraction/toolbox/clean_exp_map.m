function ephysLog_out = clean_exp_map(ephysLog)

% Remove sessions from which I couldn't find the raw data dir
ephysLog_out = ephysLog(1:2:end,:);

end