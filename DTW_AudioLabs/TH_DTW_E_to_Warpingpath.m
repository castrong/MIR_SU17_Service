function [ WarpingPath ] = TH_DTW_E_to_Warpingpath( E, Parameter )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Thomas Helten 
%Date: 2008/11/21
%Version: 1.1
%
%This function computes a warping path based on the provided matrix E and
%the allowed steps.
%
%   The first argument must be a NxM matrix of type in32 containing the
%   indices of the steps.
%
%   The second argument is a struct defining optional parameter
%       dn : 1xS integer array defining valid steps (N direction of C).
%            Default is [1 1 0].
%       dm : 1xS integer array defining valid steps (M direction of C).
%            Default is [1 0 1].
%       EndIndex : In case of subsequence DTW.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('You must at least provide a cost matrix.');
end

if ~isinteger(E)
    error('The step matrix must be of type integer');
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



if ~isinteger(dn) || any(dn < 0) || all(dn == 0)
    error('The dn parameter must contain only integers greater or equal zero and at least one integer greater zero');
end

if ~isinteger(dm) || any(dm < 0) || all(dm == 0)
    error('The dm parameter must contain only integers greater or equal zero and at least one integer greater zero');
end

N = int32(size(E,1));
M = int32(size(E,2));
S = int32(size(dn,2));

if S ~= size(dm,2)
    error('The parameters dn and dm must be of equal length.');
end

if ~isfield(Parameter,'EndIndex')
   Parameter.EndIndex = int32(M);  
end

if ~isfield(Parameter,'SubSequence')
   Parameter.SubSequence = false;
else
    if ~islogical(Parameter.SubSequence) || numel(Parameter.SubSequence) ~= 1
        error('Parameter field SubSequence is not a logical scalar!');
    end
end

if Parameter.EndIndex ~= int32(M) && ~Parameter.SubSequence
    error('DTW_Toolbox:Parameter:MustBeSet', ...
            ['The optional parameter field SubSequence must be set ', ...
             'to logical true if sub sequence DTW should be used!']);
    %%Parameter.SubSequence = true;
end

if ~isinteger(Parameter.EndIndex) || Parameter.EndIndex > M || Parameter.EndIndex < 1
    error('EndIndex is not an integer or out of range [1 M].');
end

m = Parameter.EndIndex;
n = N;

if E(n,m) == 0
    error(['A warping path ending in [' num2str(n) ', ' num2str(m) '] is not possible.']);
end

WarpingPath = zeros(2,n+m);

Index = 1;

if Parameter.SubSequence
    while n > 1 && E(n,m) > 0
        WarpingPath(:,Index) = [n m]';
        StepIndex = E(n,m);
        m = m-dm(StepIndex);
        n = n-dn(StepIndex);
        Index = Index+1;
    end
else
    while E(n,m) > 0 && (m > 1 || n > 1)
        WarpingPath(:,Index) = [n m]';
        StepIndex = E(n,m);
        m = m-dm(StepIndex);
        n = n-dn(StepIndex);
        Index = Index+1;
    end
end

WarpingPath(:,Index) = [n m]';
WarpingPath = WarpingPath(:,Index:-1:1);

end
