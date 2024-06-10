function bl = scan_protocol_blocks(startrecn,endrecn)


% function bl = scan_protocol_blocks(startrecn,endrecn)
%
% Loads block data from the protocol server for all recent records between startrecn and endrecn.
% Record numbers are located to the right on the protocol server search page.
%
% See the description in parse_block_url.m for details on the structure bl. 
%
% See also PARSE_BLOCK_URL PULLDATA


% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/scan_protocol_blocks.m $
% $Revision: 213 $
% $Date: 2013-06-04 19:48:04 -0500 (Tue, 04 Jun 2013) $
% $Author: ckovach $
% ------------------------------------------------


url = 'http://cruprotocol.neurosurgery.uiowa.edu/append.aspx?RecordID=%i'

if nargin >1
	recns = startrecn:endrecn;
else
	recns = startrecn;
end


%%
for i = 1:length(recns)
    
    recs(i) = parse_block_url(sprintf(url,recns(i)));
    fprintf('\nscanning %i',recns(i));
end

[uqr,~,rn] = unique({recs.identifier});

for i = 1:length(recs)    
    bl(rn(i)) = recs(i);
end

