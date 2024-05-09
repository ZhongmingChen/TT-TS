clear;
clc;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for RGB image
R1 = [1 50 3 1]; R2 = [1 100 3 1];
rng(100)
sigma1 = 0; 
tol1 = 1e-5; 
maxiters1 = 50;
J1 = 0.1 * prod(10*ones(1,length(R1)-1));
J2 = 0.2 * prod(10*ones(1,length(R1)-1));
J3 = 0.3 * prod(10*ones(1,length(R1)-1));

path = fullfile(pwd, 'data_sets','pompoms_ms','pompoms_ms','pompoms_RGB.bmp');
I = imread(path);
X = double(I);

fprintf('\n------Running with sigma1 = 0-------\n')
%% sketch size = 100
[G1, err1, rel1, T1] = tt_ts(X, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% sketch size = 200
[G3, err3, rel3, T3] = tt_ts(X, R1, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G4, err4, rel4, T4] = tt_random(X, R1, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% sketch size = 300
[G5, err5, rel5, T5] = tt_ts(X, R1, J3, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G6, err6, rel6, T6] = tt_random(X, R1, J3, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% numerical experiment3 plot Figure 7
tiledlayout(2,4);
nexttile(1);imshow(I);title('Original image with noise = 0');
nexttile(2); vecG2 = tt_contraction(G2); I2=uint8(reshape(vecG2,size(I)));imshow(I2); title('TT- Random with sketch size= 100, \sigma= 0');
nexttile(3); vecG4 = tt_contraction(G4); I4=uint8(reshape(vecG4,size(I)));imshow(I4); title('TT- Random with sketch size= 200, \sigma= 0');
nexttile(4); vecG6 = tt_contraction(G6); I6=uint8(reshape(vecG6,size(I)));imshow(I6); title('TT- Random with sketch size= 300, \sigma= 0');
nexttile(5);imshow(I);title('Original image with noise = 0');
nexttile(6); vecG1 = tt_contraction(G1); I1=uint8(reshape(vecG1,size(I)));imshow(I1); title('TT- TS with sketch size= 100, \sigma= 0');
nexttile(7); vecG3 = tt_contraction(G3); I3=uint8(reshape(vecG3,size(I)));imshow(I3); title('TT- TS with sketch size= 200, \sigma= 0');
nexttile(8); vecG5 = tt_contraction(G5); I5=uint8(reshape(vecG5,size(I)));imshow(I5); title('TT- TS with sketch size= 300, \sigma= 0');