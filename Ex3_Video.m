clear;
clc;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for video
R = [1 20 20 20 1]; rng(100)
sigma1 = 0;  sigma2 = 0.5;
tol1 = 1e-5; 
maxiters1 = 50;
J1 = 0.1 * prod(10*ones(1,length(R)-1));

path = fullfile(pwd, 'data_sets','skate-110734.mp4');
v1 = VideoReader(path);
frame = read(v1,[1 30]);
X = double(frame);

% run--sigma=0
fprintf('\n------Running with sigma1 = 0-------\n')
[G1, err1, rel1, T1] = tt_ts(X, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

% run--sigma=0.5
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G4, err4, rel4, T4] = tt_ts(X, R, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5, err5, rel5, T5] = tt_random(X, R, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% numerical experiment3 plot Figure 11
tiledlayout(2,3);
nexttile(1);imshow(frame(:,:,:,1));title('Skate(original)');
nexttile(2); vecG2 = tt_contraction(G2); I2=uint8(reshape(vecG2,size(frame)));imshow(I2(:,:,:,1)); title('TT- Random with sketch size= 1000, \sigma= 0, iteration: 50');
nexttile(3); vecG1 = tt_contraction(G1); I3=uint8(reshape(vecG1,size(frame)));imshow(I3(:,:,:,1)); title('TT- TS with sketch size= 1000, \sigma= 0, iteration: 50');

nexttile(4);imshow(frame(:,:,:,1));title('Skate(original)');
nexttile(5); vecG5 = tt_contraction(G5); I5=uint8(reshape(vecG5,size(frame)));imshow(I5(:,:,:,1)); title('TT- Random with sketch size= 1000, \sigma= 0.5, iteration: 50');
nexttile(6); vecG4 = tt_contraction(G4); I4=uint8(reshape(vecG4,size(frame)));imshow(I4(:,:,:,1)); title('TT- TS with sketch size= 1000, \sigma= 0.5, iteration: 50');