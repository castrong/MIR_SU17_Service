%% Test subsequence DTW implementation

% Set up cost matrix
addpath ../DTW_AudioLabs;
clear all;
close all;

% Create cost matrix test cases
C = {};

% Test #2
load('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat');
V = f_CENS';
Q = V(:,112:180);
load('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat');
R = f_CENS';
C{2} = 1-Q'*R;


% Run AudioLabs subsequence DTW code as reference
parameter.dn = int32([1 1 2 1]); % allowable steps
parameter.dm = int32([1 2 1 3]);
parameter.dw = [1.0 1.0 2.0 3.0];
parameter.SubSequence = true;
[S1,B1] = TH_DTW_C_to_DE(C{2},parameter);
[optCost1,optOffset1] = min(S1(end,:));
parameter.EndIndex = int32(optOffset1); % must specify endpoint
P1 = TH_DTW_E_to_Warpingpath(B1,parameter);

% load in things
pythonResults = load('testCythonDTW_subseq_python.mat')
S2 = pythonResults.accumCost;
B2 = pythonResults.stepMatrix;
optCost2 = pythonResults.endCost;
optOffset2 = pythonResults.offset + 1; % convert to matlab indexing
P2 = pythonResults.path + 1; % convert to matlab indexing

failedTests = 0;
if max(max(abs(S1-S2))) > 1e-10
    disp(sprintf('Cumulative cost matrix is incorrect'));
    failedTests = failedTests + 1;
end
if max(max(abs(uint8(B2) - B1))) ~= 0
    disp(sprintf('Backtrace matrix is incorrect'));
    failedTests = failedTests + 1;
end
if abs(optCost2 - optCost1) > 1e-10
    disp(sprintf('Optimal cost incorrect: %f',optCost2));
    failedTests = failedTests + 1;
end
if optOffset1 ~= optOffset2
    disp(sprintf('Optimal offset incorrect: %d',optOffset2));
    failedTests = failedTests + 1;
end
if ~isequal(P1,P2)
    disp(sprintf('Backtrace path incorrect'));
    failedTests = failedTests + 1;
end
if failedTests == 0
    disp(sprintf('Subseq DTW implementation passes subseq test'));
end

% Run AudioLabs non subsequence DTW code as reference
parameter.dn = int32([1 1 2 1]); % allowable steps
parameter.dm = int32([1 2 1 3]);
parameter.dw = [1.0 1.0 2.0 3.0];
parameter.SubSequence = false;
[S1,B1] = TH_DTW_C_to_DE(C{2},parameter);
[optCost1,optOffset1] = min(S1(end,:));
optCost1 = S1(end,end);
optOffset1 = length(S1(end,:));

parameter.EndIndex = int32(optOffset1); % must specify endpoint
P1 = TH_DTW_E_to_Warpingpath(B1,parameter);

% load in things
pythonResults = load('testCythonDTW_nonSubseq_python.mat');
S2 = pythonResults.accumCost;
B2 = pythonResults.stepMatrix;
optCost2 = pythonResults.endCost;
optOffset2 = pythonResults.offset + 1; % convert to matlab indexing
P2 = pythonResults.path + 1; % convert to matlab indexing

failedTests = 0;
if max(max(abs(S1-S2))) > 1e-10
    disp(sprintf('Cumulative cost matrix is incorrect'));
    failedTests = failedTests + 1;
end
if max(max(abs(uint8(B2) - B1))) ~= 0
    disp(sprintf('Backtrace matrix is incorrect'));
    failedTests = failedTests + 1;
end
if abs(optCost2 - optCost1) > 1e-10
    disp(sprintf('Optimal cost incorrect: %f',optCost2));
    failedTests = failedTests + 1;
end
if optOffset1 ~= optOffset2
    disp(sprintf('Optimal offset incorrect: %d',optOffset2));
    failedTests = failedTests + 1;
end
if ~isequal(P1,P2)
    disp(sprintf('Backtrace path incorrect'));
    failedTests = failedTests + 1;
end
if failedTests == 0
    disp(sprintf('Subseq DTW implementation passes non subseq test'));
end

