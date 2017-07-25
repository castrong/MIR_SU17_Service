addpath('../DTW_MIRLab')


% load the data
load('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat');
V = f_CENS';
Q = V(:,112:180);
load('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat');
R = f_CENS';
C{2} = 1-Q'*R;

% Run MIRlab subseq DTW code
stepsQ = int32([1 1 2 1]);
stepsR = int32([1 2 1 3]);
weights = [1.0 1.0 2.0];
subsequence = true;

numIterations = 100;
firstStepTimes = [];
secondStepTimes = [];
totalTimes = [];

% the two steps - time these
for i=1:numIterations
	startTic = tic;
	[S2,optCost2,B2,optOffset2] = DTW_New_costMatrix(C{2},stepsQ,stepsR,weights, subsequence);
	firstStepTime = toc(startTic);
	secondStepStartTic = tic;
	P2 = DTW_backtrace(B2,stepsQ,stepsR,optOffset2);
	secondStepTime = toc(secondStepStartTic);

	firstStepTimes = [firstStepTimes firstStepTime];
	secondStepTimes = [secondStepTimes secondStepTime];
	totalTimes = [totalTimes firstStepTime+secondStepTime];
end
fprintf('Ran subsequence test for %g iterations\n', numIterations);
fprintf('First Step: %f\n', mean(firstStepTimes));
fprintf('Second Step: %f\n', mean(secondStepTimes));
fprintf('Total Times: %f\n', mean(totalTimes));


% Run MIRlab subseq DTW code
stepsQ = int32([1 1 2 1]);
stepsR = int32([1 2 1 3]);
weights = [1.0 1.0 2.0];
subsequence = false;

numIterations = 100;
firstStepTimes = [];
secondStepTimes = [];
totalTimes = [];

% the two steps - time these
for i=1:numIterations
	startTic = tic;
	[S2,optCost2,B2,optOffset2] = DTW_New_costMatrix(C{2},stepsQ,stepsR,weights, subsequence);
	firstStepTime = toc(startTic);
	secondStepStartTic = tic;
	P2 = DTW_backtrace(B2,stepsQ,stepsR,optOffset2);
	secondStepTime = toc(secondStepStartTic);

	firstStepTimes = [firstStepTimes firstStepTime];
	secondStepTimes = [secondStepTimes secondStepTime];
	totalTimes = [totalTimes firstStepTime+secondStepTime];
end
fprintf('Ran non subsequence test for %g iterations\n', numIterations);
fprintf('First Step: %f\n', mean(firstStepTimes));
fprintf('Second Step: %f\n', mean(secondStepTimes));
fprintf('Total Times: %f\n', mean(totalTimes));