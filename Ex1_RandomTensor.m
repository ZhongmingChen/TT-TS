clear;
clc;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for randomly generated tensor
rng(60)
R = [1 5*ones(1,5) 1];
sigma1 = 0; sigma2 = 0.5;
tol1 = 1e-16; 
maxiters1 = 300; maxiters2 = 500;
J1 = 0.00015 * prod(10*ones(1,length(R)-1));
J2 = 0.0002 * prod(10*ones(1,length(R)-1));

% generate dense tensor
fprintf('Generating dense tensor X1...\n ');
[X1, cores1] = generate_low_rank_ttensor([10 10 10 10 10 10], [1 5*ones(1,5) 1], 0);
fprintf('Done!\n');
fprintf('Generating dense tensor X2 with noise = 0.01...\n ');
[X2, cores2] = generate_low_rank_ttensor([10 10 10 10 10 10], [1 5*ones(1,5) 1], 0.01);
fprintf('Done!\n');
fprintf('Generating dense tensor X3 with noise = 0.1...\n ');
[X3, cores3] = generate_low_rank_ttensor([10 10 10 10 10 10], [1 5*ones(1,5) 1], 0.1);
fprintf('Done!\n');

%% sketch size = 150
% run--sigma=0
fprintf('\n------Running with sigma1 = 0-------\n')
[G1, err1, rel1, T1] = tt_ts(X1, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X1, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G3, err3, rel3, T3] = tt_als(X1, R, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

fprintf('\n------Running with noise = 0.01-------\n')
[G1_1, err1_1, rel1_1, T1_1] = tt_ts(X2, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G2_1, err2_1, rel2_1, T2_1] = tt_random(X2, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G3_1, err3_1, rel3_1, T3_1] = tt_als(X2, R, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);

fprintf('\n------Running with noise = 0.1-------\n')
[G1_2, err1_2, rel1_2, T1_2] = tt_ts(X3, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G2_2, err2_2, rel2_2, T2_2] = tt_random(X3, R, J1, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G3_2, err3_2, rel3_2, T3_2] = tt_als(X3, R, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);


%% sketch size = 200
% run--sigma=0
fprintf('\n------Running with sigma1 = 0-------\n')
[G21, err21, rel21, T21] = tt_ts(X1, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G22, err22, rel22, T22] = tt_random(X1, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

fprintf('\n------Running with noise = 0.01-------\n')
[G21_1, err21_1, rel21_1, T21_1] = tt_ts(X2, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G22_1, err22_1, rel22_1, T22_1] = tt_random(X2, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);

fprintf('\n------Running with noise = 0.1-------\n')
[G21_2, err21_2, rel21_2, T21_2] = tt_ts(X3, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);
[G22_2, err22_2, rel22_2, T22_2] = tt_random(X3, R, J2, sigma1,'tol', tol1, 'maxiters', maxiters2, 'verbose', true);

%% numerical experiment1 plot Figure 2
iter1 = [0:maxiters1]; iter2 = [0:maxiters2]; 

%% sketch size = 150
figure;
 subplot(231);p1=plot(iter1,[0 log10(rel1)],'b->','MarkerSize',8,'LineWidth',1);p1.MarkerSize = 8; p1.MarkerIndices = 1:15:length(rel1);hold on;p1=plot(iter1,[0 log10(rel2)],'r-*','MarkerSize',8,'LineWidth',1);p1.MarkerSize = 8; p1.MarkerIndices = 1:25:length(rel2); plot(iter1,[log10(min(rel3))*ones(1,length(iter1))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 150, \sigma = 0, noise = 0')
  legend('TT-TS','TT-Random','TT-ALS')
 
 subplot(232);p2=plot(iter2,[0 log10(rel1_1)],'b->','MarkerSize',8,'LineWidth',1);p2.MarkerSize = 8; p2.MarkerIndices = 1:15:length(rel1_1);hold on;p2=plot(iter2,[0 log10(rel2_1)],'r-*','MarkerSize',8,'LineWidth',1);p2.MarkerSize = 8; p2.MarkerIndices = 1:25:length(rel2_1); plot(iter2,[log10(min(rel3_1))*ones(1,length(iter2))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 150, \sigma = 0, noise = 0.01')
  legend('TT-TS','TT-Random','TT-ALS')
 
 subplot(233);p3=plot(iter2,[0 log10(rel1_2)],'b->','MarkerSize',8,'LineWidth',1);p3.MarkerSize = 8; p3.MarkerIndices = 1:15:length(rel1_2);hold on;p3=plot(iter2,[0 log10(rel2_2)],'r-*','MarkerSize',8,'LineWidth',1);p3.MarkerSize = 8; p3.MarkerIndices = 1:25:length(rel2_2); plot(iter2,[log10(min(rel3_2))*ones(1,length(iter2))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 150, \sigma = 0, noise = 0.1')
 legend('TT-TS','TT-Random','TT-ALS')
 
 %% sketch size = 200
  subplot(234);p1=plot(iter1,[0 log10(rel21)],'b->','MarkerSize',8,'LineWidth',1);p1.MarkerSize = 8; p1.MarkerIndices = 1:15:length(rel21);hold on;p1=plot(iter1,[0 log10(rel22)],'r-*','MarkerSize',8,'LineWidth',1);p1.MarkerSize = 8; p1.MarkerIndices = 1:25:length(rel22); plot(iter1,[log10(min(rel3))*ones(1,length(iter1))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 200, \sigma = 0, noise = 0')
  legend('TT-TS','TT-Random','TT-ALS')
 
 subplot(235);p2=plot(iter2,[0 log10(rel21_1)],'b->','MarkerSize',8,'LineWidth',1);p2.MarkerSize = 8; p2.MarkerIndices = 1:15:length(rel21_1);hold on;p2=plot(iter2,[0 log10(rel22_1)],'r-*','MarkerSize',8,'LineWidth',1);p2.MarkerSize = 8; p2.MarkerIndices = 1:25:length(rel22_1); plot(iter2,[log10(min(rel3_1))*ones(1,length(iter2))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 200, \sigma = 0, noise = 0.01')
  legend('TT-TS','TT-Random','TT-ALS')
 
 subplot(236);p3=plot(iter2,[0 log10(rel21_2)],'b->','MarkerSize',8,'LineWidth',1);p3.MarkerSize = 8; p3.MarkerIndices = 1:15:length(rel21_2);hold on;p3=plot(iter2,[0 log10(rel22_2)],'r-*','MarkerSize',8,'LineWidth',1);p3.MarkerSize = 8; p3.MarkerIndices = 1:25:length(rel22_2); plot(iter2,[log10(min(rel3_2))*ones(1,length(iter2))],'g-','LineWidth',2);
 xlabel('Iteration steps');
 ylabel('Relative error(log10)');
 title('Sketch size = 200, \sigma = 0, noise = 0.1')
 legend('TT-TS','TT-Random','TT-ALS')

%% numerical experiment1 plot Figure 3
% Time v.s. error
maxiters1 = 100;

%% run--sigma=0.5

%% sketch size = 150 
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G1, err1, rel1, T1] = tt_ts(X1, R, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2, err2, rel2, T2] = tt_random(X1, R, J1, sigma2,'tol', tol1, 'maxiters', 1000, 'verbose', true);
[G3, err3, rel3, T3] = tt_als(X1, R, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

fprintf('\n------Running with noise = 0.01-------\n')
[G1_1, err1_1, rel1_1, T1_1] = tt_ts(X2, R, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2_1, err2_1, rel2_1, T2_1] = tt_random(X2, R, J1, sigma2,'tol', tol1, 'maxiters', 500, 'verbose', true);
[G3_1, err3_1, rel3_1, T3_1] = tt_als(X2, R, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

fprintf('\n------Running with noise = 0.1-------\n')
[G1_2, err1_2, rel1_2, T1_2] = tt_ts(X3, R, J1, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G2_2, err2_2, rel2_2, T2_2] = tt_random(X3, R, J1, sigma2,'tol', tol1, 'maxiters', 500, 'verbose', true);
[G3_2, err3_2, rel3_2, T3_2] = tt_als(X3, R, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);

%% sketch size = 200
fprintf('\n------Running with sigma1 = 0.5-------\n')
[G4, err4, rel4, T4] = tt_ts(X1, R, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5, err5, rel5, T5] = tt_random(X1, R, J2, sigma2,'tol', tol1, 'maxiters', 1000, 'verbose', true);

fprintf('\n------Running with noise = 0.01-------\n')
[G4_1, err4_1, rel4_1, T4_1] = tt_ts(X2, R, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5_1, err5_1, rel5_1, T5_1] = tt_random(X2, R, J2, sigma2,'tol', tol1, 'maxiters', 500, 'verbose', true);

fprintf('\n------Running with noise = 0.1-------\n')
[G4_2, err4_2, rel4_2, T4_2] = tt_ts(X3, R, J2, sigma2,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
[G5_2, err5_2, rel5_2, T5_2] = tt_random(X3, R, J2, sigma2,'tol', tol1, 'maxiters', 500, 'verbose', true);

 %% plot sketch size = 150
 figure;
subplot(231);plot([0 cumsum(T1)], [0 log10(rel1)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T2)], [0 log10(rel2)],'r','LineWidth',1.5);plot([0 cumsum(T3)], [0 log10(rel3)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.5])
 ylabel('Relative error(log10)'); ylim([-12 0])
 title('Sketch size = 150, \sigma = 0.5, noise = 0')
  legend('TT-TS','TT-Random','TT-ALS')

subplot(232);plot([0 cumsum(T1_1)], [0 log10(rel1_1)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T2_1)], [0 log10(rel2_1)],'r','LineWidth',1.5);plot([0 cumsum(T3_1)], [0 log10(rel3_1)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.15])
 ylabel('Relative error(log10)'); ylim([-3 0])
 title('Sketch size = 150, \sigma = 0.5, noise = 0.01')
  legend('TT-TS','TT-Random','TT-ALS')

subplot(233);plot([0 cumsum(T1_2)], [0 log10(rel1_2)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T2_2)], [0 log10(rel2_2)],'r','LineWidth',1.5);plot([0 cumsum(T3_2)], [1 log10(rel3_2)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.1])
 ylabel('Relative error(log10)'); ylim([-2 0])
 title('Sketch size = 150, \sigma = 0.5, noise = 0.1')
  legend('TT-TS','TT-Random','TT-ALS')
  
%% plot sketch size = 200
subplot(234);plot([0 cumsum(T4)], [0 log10(rel4)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T5)], [0 log10(rel5)],'r','LineWidth',1.5);plot([0 cumsum(T3)], [0 log10(rel3)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.5])
 ylabel('Relative error(log10)'); ylim([-12 0])
 title('Sketch size = 200, \sigma = 0.5, noise = 0')
  legend('TT-TS','TT-Random','TT-ALS')

subplot(235);plot([0 cumsum(T4_1)], [0 log10(rel4_1)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T5_1)], [0 log10(rel5_1)],'r','LineWidth',1.5);plot([0 cumsum(T3_1)], [0 log10(rel3_1)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.15])
 ylabel('Relative error(log10)'); ylim([-3 0])
 title('Sketch size = 200, \sigma = 0.5, noise = 0.01')
  legend('TT-TS','TT-Random','TT-ALS')

subplot(236);plot([0 cumsum(T4_2)], [0 log10(rel4_2)],'b','LineWidth',1.5);hold on;plot([0 cumsum(T5_2)], [0 log10(rel5_2)],'r','LineWidth',1.5);plot([0 cumsum(T3_2)], [0 log10(rel3_2)],'g','LineWidth',1.5);
 xlabel('Time(s)'); xlim([0 0.1])
 ylabel('Relative error(log10)'); ylim([-2 0])
 title('Sketch size = 200, \sigma = 0.5, noise = 0.1')
  legend('TT-TS','TT-Random','TT-ALS')