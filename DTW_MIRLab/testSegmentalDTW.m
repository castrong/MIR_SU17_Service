%% test segmental DTW backtrace

addpath ../DTW_AudioLabs;
clear all;
close all;

% set for each test
C = {};
refEndPts = {};
qLens = {};

% Test #1
C{1} = randi(100,[9,1000]);
refEndPts{1} = [50,100,160,250,300,600,800,825,990];
queryPathCost = .1;
queryLength = 50;
totalCost = queryPathCost * length(refEndPts);
for i=1:length(refEndPts{1})
    C{1}(i,refEndPts{1}(i)) = queryPathCost;
end
qLens{1} = int32(repmat(queryLength,size(refEndPts{1})));

% Test #2: taken from real data
load 'testdata_segmentalDTW.mat'; % loads pathCostMat and queryLengths
C{2} = pathCostMat;
refEndPts{2} = [580,898,1726,2399,2753];
qLens{2} = int32(queryLengths);

for i=1:length(C)
    
    [predEndPts,~,~] = segmentalDTW_backtrace(C{i},qLens{i});
    failedTests = 0;
    if ~isequal(refEndPts{i},predEndPts)
        disp(sprintf('Predicted end points are incorrect'));
        failedTests = failedTests + 1;
    end
    if failedTests == 0
        disp(sprintf('Segmental DTW implementation passes test #%d',i));
    end

end

