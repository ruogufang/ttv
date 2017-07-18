function Main_SPD
%Main_SPD is the main function for Sparse Perfusion Deconvolution algorithm
%
%   Ruogu Fang 6/26/2012
%   Advanced Multimedia Processing (AMP) Lab
%   Department of Electrical and Computer Engineering
%   Cornell University
%
%  Main_SPD reads a PCT data, adds random white noise and denoises it using
%  online dictionary learning and Lasso. The input and output PSNR are
%  compared, and the trained dictionary is displayed.
%
% References:
%
% 1. Ruogu Fang, Tsuhan Chen, Pina Sanelli. Towards Robust Deconvolution of
% Low-Dose Perfusion CT: Sparse Perfusion Deconvolution Using Online
% Dictionary Learning. Medical Image Analysis, Volumn 17, Issue 4, Pages
% 417-428, 2013.
%
% 2. Ruogu Fang, Tsuhan Chen, Pina Sanelli. Sparsity-Based Deconvolution of
% Low-Dose Perfusion CT Using Learned Dictionaries. MICCAI'12, The 15th
% Annual International Conference on Medical Image Computing and Computer
% Assisted Intervention, 2012. Lecture Notes in Computer Science Volume
% 7510, 2012, pp 272-280.
%
% Please cite the above papers if you use code in this SPD package.


clc; clear; close all;
warning off;

% Set paths
addpath(genpath('Utilities')); % Add the Utitlity path. Need to compile SPAMS library first
addpath('Data'); % Add data path

%% =================================
% Parameters for you to change
mA = 15; % Low-dose X-ray tube current level
isTrain = 1; % flag to indicate retrain a dictionary or use existing dictionary. isTrain = 0 - use existing dictionary; isTrain=1, train a dictionary from training samples
lambda = 0.06;   % Truncation parameter for cTSVD [Scalar 0<lambda<1]
m = 2; % factor of extension in block-circulant version
s = 6; % scaling factor for standard devaiation sigma in denoising

x0 = 1; y0 = 1; w = 512; h = 512; % Whole brain region
% x0 = 200; y0 = 200; w = 80; h = 80; % RACA
% x0 = 155; y0 = 170; w = 60; h = 170; % RMCA
% ==================================

%% Load CT Perfusion data used in the MedIA 2013 Paper

load ctp;
% Data format:
%   V: CTP data [T x X x Y] AIFx, AIFy, VOFx, VOFy: aif and vof
%   coordinations PRE: Pre-enhancement cutoff (first frame included in
%   calc's) POST: Post-enhancement cutoff (last frame included in calc's)

% noisy level
mA0 = 190; % Original X-ray tube current level
sigma = pct_mA2sigma(mA,mA0); % Compute the standard deviation of the Gaussian noise to simulate low-dose


% Compute Brain Mask
B = squeeze(mean(V(1:10,:,:),1));
Mask = pct_brainMask(B,0,120);

% Add Gaussian noise to simulate low-dose
Vn = pct_noise(permute(V,[2 3 1]),[],sigma,'g',Mask);
Vn = permute(Vn,[3 1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AIF and VOF search for unnoisy data
% (Uncomment the code below to redefine AIF and VOF)

% r = 5;           % Search radius for AIF and VOF
% [AIF0,AIFx,AIFy] = pct_aifsearch(C0, r);
% [VOF0,VOFx,VOFy] = pct_aifsearch(C0, r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Whole Brain CBF
[C0,AIF0] = pct_preprocess(V, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noiseless data
[C,AIF1] = pct_preprocess(Vn, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noisy data

% Crop data to Region of Interest
C0 = C0(:,y0:y0+h-1,x0:x0+w-1);
C = C(:,y0:y0+h-1,x0:x0+w-1);
Mask = Mask(y0:y0+h-1,x0:x0+w-1);

% Find the minimum bounding box
bb = pct_minBoundingBox(Mask);
C0 = C0(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
C = C(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Mask = Mask(bb(1,1):bb(1,2),bb(2,1):bb(2,3));

%% Method 0: circulant Truncated Singular Value Decomposition (cTSVD)
% Unnoisy
tid_TSVD = tic;
[R0, CBF0_TSVD, CBV0_TSVD, MTT0_TSVD,TTP0_TSVD] = pct_map(C0, AIF0, Mask, lambda, m);
t_TSVD = toc(tid_TSVD);

% Noisy
[R, CBF_TSVD, CBV_TSVD, MTT_TSVD,TTP_TSVD] = pct_map(C, AIF1, Mask, lambda, m);

%% Method 2: Sparse Perfusion Deconvolution (SPD)

% set parameters for SPD methods %
im = CBF0_TSVD;
imnoise = CBF_TSVD;
params.x = imnoise;
params.numThreads=2; % number of threads
params.K = 256; % dictionary size
params.iter=20;  % let us see what happens after 1000 iterations.
params.sigma = sigma/s; % estimated noise standard deviation
params.maxval = 255;
params.trainnum = 10000;
params.C = C;
params.AIF = AIF1;
params.R = R;
params.verbose = false;
params.rho = 1.05;
params.m = m;
params.beta = 3e4; % weight of prior
params.blocksize = 8;
overlap = 7; % patch overlap
params.stepsize = params.blocksize - overlap;


if isTrain
    % Training data of 2 patients load CBF_001.mat
    load('train001.mat');
    imtrain = CBF(100:300,100:450);
    load('train002.mat');
    imtrain = [imtrain;CBF(100:300,100:450)];
    params.xdic = double(imtrain);
end

%% Method 2.1:KSVD SPD
% Use KSVD for dictionary learning and OMP for sparse coding

params.lambda = 1; % weight of the input signal in KSVD package
params.J = []; % initialize energy function
params.round = 2;

if ~isTrain
    %Load pre-trained dictionary
    load D_KSVD_8x8.mat
    params.D = dict_ksvd;
end

% K-SVD deconvolution
tid_ksvd = tic;
disp('K-SVD Deconvolution');
[imout_ksvd, params] = ksvdpct(params);
J_ksvd = params.J;
t_ksvd = toc(tid_ksvd)

%% Method 2.2: Online Learning SPD
% Use online dictionary learning and Lasso for sparse coding

% Parameters
params.lambda = 2; % weight of the sparsity term in SPAMS package
params.mixture = 1; % weight of input signal
params.J = []; % initialize energy function
params.round = 5;

if ~isTrain
    %Load pre-trained dictionary
    load D_ODL_8x8.mat
    params.D = dict_online;
end

% Online deconvolution
tid_online = tic;
disp('Performing SPD with online learning ...');
[imout_online, params] = spd(params);
J_online = params.J;
t_online = toc(tid_online)


% Filter non brain tissue of the maps
dict = params.dict;
im = pct_mask(im,Mask);
imnoise = pct_mask(imnoise,Mask);
imout_ksvd = pct_mask(imout_ksvd, Mask);
imout_online = pct_mask(imout_online, Mask);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Printing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

dictimg = showdict(dict,[1 1]*params.blocksize,round(sqrt(params.K)),round(sqrt(params.K)),'lines','highcontrast');
figure; imshow(imresize(dictimg,2,'nearest'));
title('Trained dictionary');

% Cost function convergence
% figure; plot(J_ksvd); title('Cost function for KSVD');
% figure; plot(J_online); title('Cost function for ODL');


% Show maps
c = 50; % color range for ctshow function
figure;ctshow(im,Mask,[0 c]);title('High-dose');
figure;ctshow(imnoise,Mask,[0 c]); title('Low-dose');
figure;ctshow(imout_ksvd,Mask,[0 c]); title('KSVD');
figure;ctshow(imout_online,Mask,[0 c]); title('Online SPD');

% MSE to reference standard (cTSVD @ 190mA)
imnoise_diff = imnoise-im;
imout_ksvd_diff = imout_ksvd - im;
imout_online_diff = imout_online - im;

imnoise_mse = mean(reshape(imnoise_diff(Mask).^2,[],1))
imout_ksvd_mse = mean(reshape(imout_ksvd_diff(Mask).^2,[],1))
imout_online_mse = mean(reshape(imout_online_diff(Mask).^2,[],1))

% PSNR
imnoise_psnr = 20 * log10(max(im(Mask))/sqrt(imnoise_mse))
imout_ksvd_psnr = 20 * log10(max(im(Mask))/sqrt(imout_ksvd_mse))
imout_online_psnr = 20 * log10(max(im(Mask))/sqrt(imout_online_mse))


