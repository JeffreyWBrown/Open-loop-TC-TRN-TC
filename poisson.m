%function poisson[r,t0,t1]
%   r is constant rate for spike train
%   t0 is initial time
%   t1 is final time
%  this function will generate a train of spikes using algorithme
%    t(i+1) = t(i) - log(xrand)/r
%
function [nout,ti] = poisson(r,t0,t1)
nout = 0;       %% number of out
t = -log(rand(1))/r + t0;  %% first spike
   
if t < t1
    tt(1) = t;              %% assign first spike
    nout = 1;
    while t < t1            %% continue assigning new spike times
        t = t + -log(rand(1))/r;
        tt(nout +1) = t;
        nout = nout +1;
    end
else    nout = 0 ;
end

    if nout > 0
        ti(1:nout-1) = tt(1:nout-1);  %% remove last element
        nout = nout -1;
    else
        ti = 0;
    end
    
end
