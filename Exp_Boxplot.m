% Experiment: Boxplot of Evaluation Metrics
% Figure 10
% Ruogu Fang 12/17/2014 Advanced Multimedia Laboratory

close all; clear; clc;

%% PSNR
PSNR = [2.51	14.19	16.03	14.19	14.33	20.95	20.95
5.02	16.31	16.58	16.31	16.33	23.9	23.9
1.5	14.93	13.73	14.93	15.18	17.94	17.94
3.08	12.71	14.15	12.71	12.71	19.15	19.15
4.57	11.98	12.6	11.98	12	18.34	18.34
3.99	16.37	14.62	16.37	16.39	21.6	21.6
8.17	19.05	19.58	19.05	19.44	22.03	22.03
6.35	16.99	16.64	16.99	17	24.53	24.53
6.79	18.85	19.16	18.85	19.17	22.66	22.66
6.08	18.06	15.87	18.06	18.1	22.44	22.44
9.34	18.58	19.34	18.58	18.93	23.18	23.18
3.19	17.42	13.42	17.42	17.46	21.86	21.86];
labels= {'sSVD','bSVD','Tikhonov','TIPS+bSVD','SPD','TTV','TIPS+TTV'};

h1 = figure;
set(h1, 'DefaultTextFontSize', 20);
set(gca,'FontSize',20);
bp = boxplot(PSNR,'labels',labels,'notch', 'on');
ylabel('PSNR','FontSize',20);
for i = 1:size(bp,2), set(bp(:,i),'linewidth',2); end

[h_psnr,p_psnr]=ttest2(PSNR(:,end-1),PSNR(:,end-2));

figpath = '../../../Figures/pdf';
w = 12; h = 10;

set(h1, 'papersize', [w h]);
set(h1, 'paperposition', [0 0 w h]);
print(h1,fullfile(figpath,'PSNR.pdf'),'-dpdf');

%% Lin's CCC
CCC = [0.074	0.306	0.348	0.306	0.312	0.664	0.664
0.209	0.583	0.537	0.583	0.62	0.605	0.605
0.113	0.503	0.391	0.503	0.517	0.676	0.676
0.116	0.407	0.397	0.407	0.414	0.609	0.609
0.163	0.444	0.432	0.444	0.465	0.529	0.529
0.048	0.283	0.218	0.283	0.284	0.493	0.493
0.194	0.539	0.504	0.539	0.559	0.644	0.644
0.224	0.377	0.324	0.377	0.387	0.403	0.403
0.159	0.561	0.502	0.561	0.579	0.712	0.712
0.169	0.415	0.39	0.415	0.428	0.265	0.265
0.226	0.531	0.485	0.531	0.55	0.727	0.727
0.073	0.366	0.195	0.366	0.368	0.459	0.459];
labels= {'sSVD','bSVD','Tikhonov','TIPS+bSVD','SPD','TTV','TIPS+TTV'};

h2 = figure;
set(h2, 'DefaultTextFontSize', 20);
set(gca,'FontSize',20);
bp = boxplot(CCC,'labels',labels,'notch', 'on');
ylabel('CCC','FontSize',20);
for i = 1:size(bp,2), set(bp(:,i),'linewidth',2); end

[h_ccc,p_ccc]=ttest2(CCC(:,end-1),CCC(:,end-2));


set(h2, 'papersize', [w h]);
set(h2, 'paperposition', [0 0 w h]);
print(h2,fullfile(figpath,'CCC.pdf'),'-dpdf');
