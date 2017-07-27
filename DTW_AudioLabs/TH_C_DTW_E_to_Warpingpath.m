function [ WarpingPath ] = TH_C_DTW_E_to_Warpingpath( E, Parameter )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Thomas Helten 
%Date: 2008/11/21
%
%This function computes a warping path based on the provided matrix E and
%the allowed steps.
%
%   The first argument must be a NxM matrix of type int32 containing the
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
show_DTW_Compilation_Intructions;