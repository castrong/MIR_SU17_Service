function [ D, E ] = TH_DTW_C_to_DE( C, Parameter )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Thomas Helten 
%Date: 2008/11/21
%
%This function computes the accumulated cost matrix D and the step index
%matrix E.
%
%   The first argument must be a NxM matrix of type double containing the
%   costs.
%
%   The second argument is a struct defining optional parameter
%       dn : 1xS integer array defining valid steps (N direction of C).
%            Default is [1 1 0].
%       dm : 1xS integer array defining valid steps (M direction of C).
%            Default is [1 0 1].
%       dw : 1xS double array defining the weight of each step.
%            Default is [1 1 1].
%       SubSequence : a boolean value defining wether this function shall
%            do a subsequence match (true) or a global match (false=default). 
%
%The returned arguments are the D matrix of type double and size NxM and
%the E matrix of type uint8 and size NxM. Matrix D hold the accumulated
%costs according to the allowed steps. E(n,m) holds the index of the step
%take to determine the value of D(n,m). If E(n,m) is zero, no valid step
%was possible. NaNs in the cost matrix are preserved, invalid fields in the
%cost matrix are NaNs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('You must at least provide a cost matrix.');
end

if ~isnumeric(C) || ~isreal(C)
    error('The cost matrix must be numeric but not complex');
end

if ~exist('Parameter','var')
    Parameter = struct();
end

if isfield(Parameter,'dn')
    dn = Parameter.dn;
else
    dn = int32([1 1 0]);
end

if isfield(Parameter,'dm')
    dm = Parameter.dm;
else
    dm = int32([1 0 1]);
end

if isfield(Parameter,'dw')
    dw = Parameter.dw;
else
    dw = [1 1 1];
end

if ~isinteger(dn) || any(dn < 0) || all(dn == 0)
    error('The dn parameter must contain only integers greater or equal zero and at least one integer greater zero');
end

if ~isinteger(dm) || any(dm < 0) || all(dm == 0)
    error('The dm parameter must contain only integers greater or equal zero and at least one integer greater zero');
end

if ~isnumeric(dw) || any(dw <= 0)
    error('The dw parameter must contain only numbers greater than zero');
end

N = int32(size(C,1));
M = int32(size(C,2));
S = int32(size(dn,2));

if S ~= size(dm,2) || S ~=  size(dw,2)
    error('The parameters dn,dm, and dw must be of equal length.');
end

if ~isfield(Parameter,'SubSequence')
   Parameter.SubSequence = false;
else
    if ~islogical(Parameter.SubSequence)
        error('Parameter SubSequence is not of type logical');
    end
end

%% calc bounding box size of steps
sbbn = max(dn);
sbbm = max(dm);

%% initialize E
E = uint8(zeros(N,M));

%% initialize extended D matrix
D = ones(sbbn+N,sbbm+M)*nan;

if Parameter.SubSequence
    for m=1:M
        D(sbbn+1,sbbm+m) = C(1,m);
    end
else
    D(sbbn+1,sbbm+1) = C(1,1);
end


%% accumulate
for m=(1:M)+sbbm
    for n=(1:N)+sbbn
        for s=1:S
            cost = D(n-dn(s),m-dm(s))+C(n-sbbn,m-sbbm)*dw(s);
            [D(n,m),Idx] = min([D(n,m) cost]);
            if Idx == 2
               E(n-sbbn,m-sbbm) = s; 
            end
        end
    end
end

D = D((1:N)+sbbn,(1:M)+sbbm);

end
