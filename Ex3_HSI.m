clear;
clc;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for hyperspectral image 
R1 = [1 20 20 1]; rng(50)
sigma1 = 0;  sigma2 = 0.5;
tol1 = 1e-5; 
maxiters1 = 50;
J1 = 0.5 * prod(10*ones(1,length(R1)-1));
J2 = 1 * prod(10*ones(1,length(R1)-1));

path = fullfile(pwd, 'data_sets','Indian_pines.mat');
load(path);
X = double(indian_pines);

%% Sketch size = 500
% run--sigma=0
fprintf('\n------Running with sigma1 = 0-------\n')
[G1, err1, rel1, T1] = tt_ts(X, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G3, err3, rel3, T3] = tt_als(X, R1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

% run--sigma=0.5
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G4, err4, rel4, T4] = tt_ts(X, R1, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5, err5, rel5, T5] = tt_random(X, R1, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G6, err6, rel6, T6] = tt_als(X, R1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% Sketch size = 1000
% run--sigma=0
fprintf('\n------Running with sigma1 = 0-------\n')
[G1_1, err1_1, rel1_1, T1_1] = tt_ts(X, R1, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2_1, err2_1, rel2_1, T2_1] = tt_random(X, R1, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

% run--sigma=0.5
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G4_1, err4_1, rel4_1, T4_1] = tt_ts(X, R1, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5_1, err5_1, rel5_1, T5_1] = tt_random(X, R1, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% numerical experiment3 plot Figure 9
%%  sigma = 0
figure; %
subplot(121); [width,height,band] = size(X); N = width * height;data = reshape(X,N,band);plot(data(1,:),'m','LineWidth',1);hold on;
vecG1 = tt_contraction(G1);X1=reshape(vecG1,size(X));[width,height,band] = size(X1); N = width * height;data = reshape(X1,N,band);plot(data(1,:),'b','LineWidth',1);hold on;
vecG2 = tt_contraction(G2);X2=reshape(vecG2,size(X));[width,height,band] = size(X2); N = width * height;data = reshape(X2,N,band);plot(data(1,:),'r','LineWidth',1);hold on;
vecG3 = tt_contraction(G3);X3=reshape(vecG3,size(X));[width,height,band] = size(X3); N = width * height;data = reshape(X3,N,band);plot(data(1,:),'g','LineWidth',1);hold on;
legend('Original spectral curve','TT-TS','TT-Random','TT-ALS'); xlim([0 220]);
title('Hyperspectral image with sketch size =500, r=20, \sigma = 0, iteration: 50'); ylim([0 7000]);

subplot(122); [width,height,band] = size(X); N = width * height;data = reshape(X,N,band);plot(data(1,:),'m','LineWidth',1);hold on;
vecG1_1 = tt_contraction(G1_1);X1=reshape(vecG1_1,size(X));[width,height,band] = size(X1); N = width * height;data = reshape(X1,N,band);plot(data(1,:),'b','LineWidth',1);hold on;
vecG2_1 = tt_contraction(G2_1);X2=reshape(vecG2_1,size(X));[width,height,band] = size(X2); N = width * height;data = reshape(X2,N,band);plot(data(1,:),'r','LineWidth',1);hold on;
vecG3 = tt_contraction(G3);X3=reshape(vecG3,size(X));[width,height,band] = size(X3); N = width * height;data = reshape(X3,N,band);plot(data(1,:),'g','LineWidth',1);hold on;
legend('Original spectral curve','TT-TS','TT-Random','TT-ALS'); xlim([0 220]);
title('Hyperspectral image with sketch size =1000, r=20, \sigma = 0, iteration: 50'); ylim([0 6000]);

%% numerical experiment3 plot Figure 10
%% sigma = 0.5
figure; % 
subplot(121); [width,height,band] = size(X); N = width * height;data = reshape(X,N,band);plot(data(1,:),'m','LineWidth',1);hold on;
vecG4 = tt_contraction(G4);X4=reshape(vecG4,size(X));[width,height,band] = size(X4); N = width * height;data = reshape(X4,N,band);plot(data(1,:),'b','LineWidth',1);hold on;
vecG5 = tt_contraction(G5);X5=reshape(vecG5,size(X));[width,height,band] = size(X5); N = width * height;data = reshape(X5,N,band);plot(data(1,:),'r','LineWidth',1);hold on;
vecG6 = tt_contraction(G6);X6=reshape(vecG6,size(X));[width,height,band] = size(X6); N = width * height;data = reshape(X6,N,band);plot(data(1,:),'g','LineWidth',1);hold on;
legend('Original spectral curve','TT-TS','TT-Random','TT-ALS'); xlim([0 220]);
title('Hyperspectral image with sketch size =500, r=20, \sigma = 0.5, iteration: 50');ylim([0 6000]);

subplot(122); [width,height,band] = size(X); N = width * height;data = reshape(X,N,band);plot(data(1,:),'m','LineWidth',1);hold on;
vecG4_1 = tt_contraction(G4_1);X1=reshape(vecG4_1,size(X));[width,height,band] = size(X1); N = width * height;data = reshape(X1,N,band);plot(data(1,:),'b','LineWidth',1);hold on;
vecG5_1 = tt_contraction(G5_1);X2=reshape(vecG5_1,size(X));[width,height,band] = size(X2); N = width * height;data = reshape(X2,N,band);plot(data(1,:),'r','LineWidth',1);hold on;
vecG6 = tt_contraction(G6);X6=reshape(vecG6,size(X));[width,height,band] = size(X6); N = width * height;data = reshape(X6,N,band);plot(data(1,:),'g','LineWidth',1);hold on;
legend('Original spectral curve','TT-TS','TT-Random','TT-ALS'); xlim([0 220]);
title('Hyperspectral image with sketch size =1000, r=20, \sigma = 0.5, iteration: 50'); ylim([0 6000]);