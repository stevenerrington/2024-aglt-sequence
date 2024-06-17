function [event] = ft_read_event_BA(filename)
%%

[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV(filename, [1 1 1 1 1], 1, 1, [] );

%%
temp = find( cellfun(@length, EventStrings) <30); %find the strings with less than 30 characters. those will not have codes but just starting or ending recording
EventStrings(temp) = {'TTL Input on AcqSystem0_0 board 0 port 6 value (0x1388).'}; %the number will be 5000 when converted from hex to decimal
%%

ttls = num2cell(TTLs);
timestamp = num2cell(Timestamps);
eventID = num2cell(EventIDs);

%%
port = extractAfter(EventStrings,'port ');
port = extractBefore(port,' value');
% %%
% port = extractBetween(EventStrings,39,40); % extracting port number (0 for audio/visual and 2 for TTL codes)
% port = num2cell(hex2dec(port));
port = num2cell(hex2dec(port));
%%

% value = extractBetween(EventStrings,51,54); %extracting the codes from the string cell array
value = extractAfter(EventStrings,'(0x');
value = extractBefore(value,')');
value = num2cell(hex2dec(value));

event = struct('ttls',ttls,...
    'timestamp',timestamp,...
    'eventID',eventID,...
    'port',port',...
    'value',value',...
    'EventStrings', EventStrings');
%%

end
