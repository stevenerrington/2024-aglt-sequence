function pw = getpw(prompt)

% get a password without echoing the result
if nargin < 1
    prompt='Enter password: ';
end
f = figure('visible','off','position',[-100 -100 0.1 0.1],'windowstyle','modal');
drawnow
pw='';
currchar='x';
figure(f)
fprintf('\n%s',prompt)
while ~ismember(double(currchar),[double(sprintf('\n')),10,13])
    %figure(f)
    waitforbuttonpress
    drawnow
    currchar= get(f,'CurrentCharacter');
    if ismember(currchar,[8 double(sprintf('\n'))]) && ~isempty(pw)
        pw(end)=[];
        fprintf('\b');
    else
        
        pw(end+1) =currchar;
        fprintf('*')
    end

end
fprintf('\b\n')
close(f)
pw=deblank(pw);
