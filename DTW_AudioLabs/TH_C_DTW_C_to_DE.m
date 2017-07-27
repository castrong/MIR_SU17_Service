function [D,E] = TH_C_DTW_C_to_DE(C,Parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Thomas Helten 
%Date: 2008/11/21
%
%This function computes the accumulated cost matrix D and the step index
%matrix E.
%
%   The first argument must be a NxM matrix of type double containing the
%   costs (no NaNs allowed but Infs).
%
%   The second argument is a struct defining optional parameter
%       dn : 1xS integer array defining valid steps (n direction of C).
%            Default is [1 1 0].
%       dm : 1xS integer array defining valid steps (m direction of C).
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
%was possible.
%
%  !!! NaNs in the Cost Matrix are not preserved, they become inf      !!!
%  !!! Invalid fields in the accumulated cost matrix have the value inf!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_DTW_Compilation_Intructions;