%% Test subsequence DTW implementation

% Set up cost matrix
addpath ../DTW_AudioLabs;
clear all;
close all;

% Create cost matrix test cases
C = {};

% Test #1
C{1} = randi(100,[9,20]);
C{1}(1,5) = 0.1; % define min cost path (total cost = .9)
C{1}(2,7) = 0.1;
C{1}(3,9) = 0.1;
C{1}(4,10) = 0.1;
C{1}(5,11) = 0.1;
C{1}(7,12) = 0.1;
C{1}(9,13) = 0.1;

% Test #2
load('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat');
V = f_CENS';
Q = V(:,112:180);
load('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat');
R = f_CENS';
C{2} = 1-Q'*R;

for testCaseIdx=1:length(C),

    fprintf('==> Testing with subsequence\n');

    % Run AudioLabs subsequence DTW code as reference
    parameter.dn = int32([1 1 2]); % allowable steps
    parameter.dm = int32([1 2 1]);
    parameter.dw = [1.0 1.0 2.0];
    parameter.SubSequence = true;
    [S1,B1] = TH_DTW_C_to_DE(C{testCaseIdx},parameter);
    [optCost1,optOffset1] = min(S1(end,:));
    parameter.EndIndex = int32(optOffset1); % must specify endpoint
    P1 = TH_DTW_E_to_Warpingpath(B1,parameter);

    S1
    parameter.EndIndex
    optCost1
    P1 - 1


    % Run MIRlab subseq DTW code
    stepsQ = int32([1 1 2]);
    stepsR = int32([1 2 1]);
    weights = [1.0 1.0 2.0];
    [S2,optCost2,B2,optOffset2] = DTW_New_costMatrix(C{testCaseIdx},stepsQ,stepsR,weights, true);
    P2 = DTW_backtrace(B2,stepsQ,stepsR,optOffset2);

    % Compare results to reference implementation
    failedTests = 0;
    if max(max(abs(S1-S2))) ~= 0
        disp(sprintf('Cumulative cost matrix is incorrect'));
        failedTests = failedTests + 1;
    end
    if max(max(abs(B2 - int8(B1)))) ~= 0
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
        disp(sprintf('Subseq DTW implementation passes test suite #%d, subsequence',testCaseIdx));
    end

    fprintf('==> Testing without subsequence\n');

    % NOW RUN WITH NON-SUBSEQUENCE
    % Run AudioLabs DTW code as reference
    parameter.dn = int32([1 1 2 1]); % allowable steps
    parameter.dm = int32([1 2 1 3]);
    parameter.dw = [1.0 1.0 2.0 3.0];
    parameter.SubSequence = false;
    [S1,B1] = TH_DTW_C_to_DE(C{testCaseIdx},parameter);
    % pick the best cost as the corner
    optCost1 = S1(end,end);
    optOffset1 = length(S1(end,:));
    parameter.EndIndex = int32(optOffset1); % must specify endpoint
    P1 = TH_DTW_E_to_Warpingpath(B1,parameter);


    % Run MIRlab DTW code
    stepsQ = int32([1 1 2 1]);
    stepsR = int32([1 2 1 3]);
    weights = [1.0 1.0 2.0 3.0];
    [S2,optCost2,B2,optOffset2] = DTW_New_costMatrix(C{testCaseIdx},stepsQ,stepsR,weights, false);
    P2 = DTW_backtrace(B2,stepsQ,stepsR,optOffset2);

    % Compare results to reference implementation
    failedTests = 0;
    if max(max(abs(S1-S2))) ~= 0
        disp(sprintf('Cumulative cost matrix is incorrect'));
        failedTests = failedTests + 1;
    end
    if max(max(abs(B2 - int8(B1)))) ~= 0
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
        disp(sprintf('Subseq DTW implementation passes test suite #%d, non-subseqeuence',testCaseIdx));
    end


end
