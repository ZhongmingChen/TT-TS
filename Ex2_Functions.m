clear;
clc;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for one-dimensional function
rng(50)
R1 = [1 5*ones(1,5) 1]; R2 = [1 20*ones(1,5) 1];
sigma1 = 0; sigma2 = 0.5;
tol1 = 1e-10; 
maxiters1 = 100;
J1 = 0.000045 * prod(10*ones(1,length(R1)-1));
J2 = 0.001 * prod(10*ones(1,length(R2)-1));

% generate tensor
x1 = linspace(-5,5,10^6);
y1 = sinc(x1);
X1 = reshape(y1,10*ones(1,6));

x2 = linspace(0.01,1,10^6);
y2 = sin(4./x2).*cos(x2.^2);
X2 = reshape(y2,10*ones(1,6));

%% run X1 with sigma = 0
% Run X1 with sketch size = 45
fprintf('\n------Running with X1-------\n')
[G1, err1, rel1, T1] = tt_ts(X1, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X1, R1, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

% Run X2 with  sketch size = 1000
fprintf('\n------Running with X2-------\n')
[G3, err3, rel3, T3] = tt_ts(X2, R2, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G4, err4, rel4, T4] = tt_random(X2, R2, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% run X2 with sigma = 0.5
% Run X2 with sketch size = 45
fprintf('\n------Running with sigma1 = 0-------\n')
[G1_1, err1_1, rel1_1, T1_1] = tt_ts(X1, R1, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2_1, err2_1, rel2_1, T2_1] = tt_random(X1, R1, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G, err, rel, T] = tt_als(X1, R1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

% Run X2 with  sketch size = 1000
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G3_1, err3_1, rel3_1, T3_1] = tt_ts(X2, R2, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G4_1, err4_1, rel4_1, T4_1] = tt_random(X2, R2, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G_1, err_1, rel_1, T_1] = tt_als(X2, R2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% numerical experiment2 plot Figure 4
figure; % sigma = 0
subplot(231);plot(x1,y1);hold on; title('y = sinc(x)')
subplot(232);vecG2 = tt_contraction(G2);plot(x1,vecG2,'-','LineWidth',1);hold on;legend('TT-Random'); title('Sketch size =45, r= 5, \sigma= 0, iteration: 100')
subplot(233);vecG1 = tt_contraction(G1);plot(x1,vecG1,'LineWidth',1);hold on;legend('TT-TS');title('Sketch size =45, r= 5, \sigma= 0, iteration: 100');ylim([-0.4 1])
subplot(234);plot(x2,y2);hold on; title('y = sin(4./x).*cos(x^2)')
subplot(235);vecG4 = tt_contraction(G4);plot(x2,vecG4,'LineWidth',1);hold on;legend('TT-Random');title('Sketch size =1000, r= 20, \sigma= 0, iteration: 100')
subplot(236);vecG3 = tt_contraction(G3);plot(x2,vecG3,'LineWidth',1);hold on;legend('TT-TS');title('Sketch size =1000, r= 20, \sigma= 0, iteration: 100');ylim([-1 1])

%% numerical experiment2 plot Figure 5
figure; % sigma = 0.5
subplot(231);plot(x1,y1);hold on; title('y = sinc(x)')
subplot(232);vecG2 = tt_contraction(G2_1);plot(x1,vecG2,'-','LineWidth',1);hold on;legend('TT-Random'); title('Sketch size =45, r= 5, \sigma= 0.5, iteration: 100')
subplot(233);vecG1 = tt_contraction(G1_1);plot(x1,vecG1,'LineWidth',1);hold on;legend('TT-TS');title('Sketch size =45, r= 5, \sigma= 0.5, iteration: 100');ylim([-0.4 1])
subplot(234);plot(x2,y2);hold on; title('y = sin(4./x).*cos(x^2)')
subplot(235);vecG4 = tt_contraction(G4_1);plot(x2,vecG4,'LineWidth',1);hold on;legend('TT-Random');title('Sketch size =1000, r= 20, \sigma= 0.5, iteration: 100')
subplot(236);vecG3 = tt_contraction(G3_1);plot(x2,vecG3,'LineWidth',1);hold on;legend('TT-TS');title('Sketch size =1000, r= 20, \sigma= 0.5, iteration: 100');ylim([-1 1])

%% numerical experiment2 plot Figure 6
figure; % time
subplot(121);plot([0 cumsum(T1_1)], [0 log10(rel1_1)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T2_1)], [0 log10(rel2_1)],'r','LineWidth',1.5); plot([0 cumsum(T)], [0 log10(rel)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.025])
 ylabel('Relative error(log10)'); ylim([-3.5 0.2])
 title('y=sinc(x), sketch size = 45, r=5, \sigma = 0.5')
 legend('TT-TS','TT-Random','TT-ALS')
  
 subplot(122);plot([0 cumsum(T3_1)], [0 log10(rel3_1)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T4_1)], [0 log10(rel4_1)],'r','LineWidth',1.5); plot([0 cumsum(T_1)], [0 log10(rel_1)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 1.6])
 ylabel('Relative error(log10)'); ylim([-10 0.2])
 title('y = sin(4./x).*cos(x^2), sketch size = 1000, r=20, \sigma = 0.5')
 legend('TT-TS','TT-Random','TT-ALS')