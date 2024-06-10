function B = chooseperm(N, k, item, rows)

% function B = chooseperm( N, k )
% 
%    Returns indices for every unordered combination of k items from a
%    population of N using a clever recursive algorithm.
%   
%    B is a C x k matrix where C is the binomial coefficient (N,k)
% 
%   Example:
% 
%                 >> B = chooseperm(4,3)
% 
%                 B =
%                      1     2     3
%                      1     2     4
%                      1     3     4
%                      2     3     4
% 
% 
%   See also NCHOOSEK

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/misctools/chooseperm.m $
% $Revision: 36 $
% $Date: 2011-06-04 23:10:51 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------


%
% C. Kovach 2007
%

if nargin < 3
    rows = 1;
    item = 1;
end

B = [];

if item >= k
    
    B = (rows:N)';
    
else
    
    for strow = rows:N

       b = chooseperm( N, k, item+1, strow+1);
    
       if ~isempty(b)
           
           B = cat(1,B,cat(2,ones(size(b,1),1)*strow,b));
           
       else
           
           return
           
       end

    end
    
end

