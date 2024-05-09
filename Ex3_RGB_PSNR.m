clc;
clear;
addpath(genpath('tensor_toolbox-master1'));
addpath(genpath('tucker-tensorsketch-master'));

%% numerical experiment for RGB image
path = fullfile(pwd, 'data_sets','pompoms_ms','pompoms_ms','pompoms_RGB.bmp');
I = imread(path);
X = double(I);
R1 = [1 50 3 1]; R2 = [1 100 3 1];
rng(220)
sigma1 = 0; 
tol1 = 1e-5; 
maxiters1 = 50;
psnr1=zeros(1,5); psnr2=zeros(1,5);  psnr4=zeros(1,5);  psnr5=zeros(1,5);
%% TT-ranks = (1,50,3,1)
for i =250:250:1000
      [G1, err1, rel1, T1] = tt_ts(X, R1, i, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
      vecG1 = tt_contraction(G1);J=uint8(reshape(vecG1,size(I)));
      psnr2(i/250+1)=PSNR(I,J);
end
for i =250:250:1000
      [G2, err2, rel2, T2] = tt_random(X, R1, i, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
      vecG2 = tt_contraction(G2);J=uint8(reshape(vecG2,size(I)));
      psnr1(i/250+1)=PSNR(I,J);
end
[G3, err3, rel3, T3] = tt_als(X, R1, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG3 = tt_contraction(G3);J=uint8(reshape(vecG3,size(I))); psnr3=PSNR(I,J)*ones(1,5); 
[G1_1, err1_1, rel1_1, T1_1] = tt_ts(X, R1, 100, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG1 = tt_contraction(G1_1);J=uint8(reshape(vecG1,size(I))); psnr2(1)=PSNR(I,J);
[G2_1, err2_1, rel2_1, T2_1] = tt_random(X, R1, 100, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG2 = tt_contraction(G2_1);J=uint8(reshape(vecG2,size(I))); psnr1(1)=PSNR(I,J);

%% numerical experiment3 plot Figure 8
figure;
%bar1
subplot(121)
%data
Y1=[psnr1(1),psnr2(1),psnr3(1);psnr1(2),psnr2(2),psnr3(2);psnr1(3),psnr2(3),psnr3(3);psnr1(4),psnr2(4),psnr3(4);psnr1(5),psnr2(5),psnr3(5)];
X1=1:5;
h1=bar(X1,Y1,1);
set(gca,'XTickLabel',{'100','250','500','750','1000'},'FontSize',10,'FontName','Times New Roman');
set(h1(1),'FaceColor',[162,214,249]/255)
set(h1(2),'FaceColor',[30,150,252]/255)
set(h1(3),'FaceColor',[252,243,0]/255)
ylim([20,40]); 
ylabel('PSNR(dB)');
xlabel('Sketch size');
legend('TT-Random','TT-TS','TT-ALS');
title('r = (1,50,3,1)');
set(gca,'FontSize',11);

%% TT-ranks = (1,100,3,1)
clear Gs_hat;
for i =250:250:1000
      [G4, err4, rel4, T4] = tt_ts(X, R2, i, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
      vecG4 = tt_contraction(G4);J=uint8(reshape(vecG4,size(I)));
      psnr5(i/250+1)=PSNR(I,J);
end
for i =250:250:1000
      [G5, err5, rel5, T5] = tt_random(X, R2, i, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true);
      vecG5 = tt_contraction(G5);J=uint8(reshape(vecG5,size(I)));
      psnr4(i/250+1)=PSNR(I,J);
end
[G6, err6, rel6, T6] = tt_als(X, R2, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG6 = tt_contraction(G6);J=uint8(reshape(vecG6,size(I))); psnr6=PSNR(I,J)*ones(1,5); 
[G4_1, err4_1, rel4_1, T4_1] = tt_ts(X, R2, 100, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG4 = tt_contraction(G4_1);J=uint8(reshape(vecG4,size(I))); psnr5(1)=PSNR(I,J);
[G5_1, err5_1, rel5_1, T5_1] = tt_random(X, R2, 100, sigma1,'tol', tol1, 'maxiters', maxiters1, 'verbose', true); vecG5 = tt_contraction(G5_1);J=uint8(reshape(vecG5,size(I))); psnr4(1)=PSNR(I,J);

%bar2 
subplot(122)
%data
Y2=[psnr4(1),psnr5(1),psnr6(1);psnr4(2),psnr5(2),psnr6(2);psnr4(3),psnr5(3),psnr6(3);psnr4(4),psnr5(4),psnr6(4);psnr4(5),psnr5(5),psnr6(5)];
X2=1:5;
h2=bar(X2,Y2,1);
set(gca,'XTickLabel',{'100','250','500','750','1000'},'FontSize',10,'FontName','Times New Roman');
set(h2(1),'FaceColor',[162,214,249]/255)
set(h2(2),'FaceColor',[30,150,252]/255)
set(h2(3),'FaceColor',[252,243,0]/255)
ylim([20,45]);      
ylabel('PSNR(dB)');
xlabel('Sketch size');
legend('TT-Random','TT-TS','TT-ALS');
title('r = (1,100,3,1)');
set(gca,'FontSize',11);