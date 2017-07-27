clear all;
close all hidden;

load('test_Shostakovich_JazzSuite2_6_Waltz2_Chailly_CENS_41_10.mat');
V = f_CENS';
load('test_Shostakovich_JazzSuite2_6_Waltz2_Yablonsky_CENS_41_10.mat');
W = f_CENS';

M = size(V,2);
N = size(W,2);

% cost matrix
C = 1-V'*W;

% allowed step sizes
% Sigma_1
parameter.dn = int32([1 0 1]);
parameter.dm = int32([0 1 1]);
% Sigma_2: (no horizontal or vertical steps)
% parameter.dn = int32([1 2 1]);
% parameter.dm = int32([2 1 1]);

% step weights
parameter.dw = [1 1 1];
% parameter.dw = [1.5 1.5 2];

showFigure = 1;



%% Run TH Matlab version
parameter.SubSequence = false;
c = cputime;
[D_TH_MAT,E_TH_MAT] = TH_DTW_C_to_DE(C,parameter);
fprintf('Running Time TH Matlab: %.5f\n',cputime-c);
if showFigure == 1
    figure(3);imagesc(D_TH_MAT); axis xy; colorbar;
    Delta_TH_Matlab = D_TH_MAT(end,:)/M;
    figure(4); imagesc(E_TH_MAT); axis xy; colorbar;
end

%% Run TH C/C++ version
parameter.SubSequence = false;
c = cputime;
[D_TH_CPP,E_TH_CPP] = TH_C_DTW_C_to_DE(C,parameter);
fprintf('Running Time TH C/C++: %.5f\n',cputime-c);
if showFigure == 1
    figure(5);imagesc(D_TH_CPP); axis xy; colorbar;
    Delta_TH_DLL = D_TH_CPP(end,:)/M;
    figure(6); imagesc(E_TH_CPP); axis xy; colorbar;
end



WP1 = TH_DTW_E_to_Warpingpath(E_TH_CPP,parameter);
WP2 = TH_DTW_E_to_Warpingpath(E_TH_CPP,parameter);

visDTW(D_TH_CPP,WP1);








