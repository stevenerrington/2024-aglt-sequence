function rth = vonmisesrand(m,k,N)

%  r = vonmisesrand(m,k,N)
%
% Generates a vector of N random draws from a von-mises distribution
% using Best and Fisher's algorithm (Best and Fisher, 1979)
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2011

if nargin < 3
    N = 1;
end

if nargin < 2
    k = 1;
end


tau = @(k) 1 + sqrt(1 + 4*k^2);

rho = @(k) (tau(k) - sqrt(2*tau(k)) )/(2*k);

%%% Step 0 
r = ( 1 + rho(k).^2 ) / (2*rho(k));

rth = zeros(1,N);
nr = 0;

while nr < N
    
    %%% generate 3 uniform rand nums.
    u = rand(1,3);
    
    
    %%% Step 1 
    z = cos( pi * u(1) );
    f = (1 + r*z)/(r+z);
    c = k*(r-f);
    
    %%% Step 3 
    if c*(2-c) - u(2) > 0  ||  log( c / u(2) ) + 1 - c >= 0
       
        %%% Step 4
        nr = nr+1;
        
        rth(nr) = sign( u(3) - .5 ) * acos(f); 
        
    end
end

rth = atan2( sin(rth + m), cos( rth + m ) );


        