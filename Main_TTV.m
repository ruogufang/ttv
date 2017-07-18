% TTV-CTP CT perfusion deconvolution using total variation regularization on
% the residue functions with non-iso regularization parameter lambda
% 
% Figure 7 & 8
%
% Ruogu Fang 9/21/2013 Advanced Multimedia Laboratory

close all; clear; clc;

% Set paths
addpath(genpath('toolbox')); % Include toolbox package

%%%%%%%%%% Setting parameters %%%%%%%%%%%%%%%%%%
mA = 15; % tube current-exposure time product
rt = 1; % downsampling rate in time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tid_all = tic;

% Parameters
m = 2; % block-circulant Ca extended for m times
lambda = 0.1; % cTSVD truncation parameter for singular values
rho = 1.05; %Average brain tissue density
Tn = 30; % truncated time in TTV

% Load data
imgname = 'CTP.mat';
load(imgname);
% Data format:
%   V: CTP data [T x X x Y] AIFx, AIFy, VOFx, VOFy: aif and vof
%   coordinations PRE: Pre-enhancement cutoff (first frame included in
%   calc's) POST: Post-enhancement cutoff (last frame included in calc's)

load acf;

% Region of Interest
x0 = 200; y0 = 200; w = 80; h = 80; % RACA
% x0 = 150; y0 = 150; w = 200; h = 200; % RMCA
% x0 = 1; y0 = 1; w = 512; h = 512; % Whole image

% noisy level
mA0 = 190;
sigma = pct_mA2sigma(mA,mA0);

% Remove negative values
V(V<0) = 0;

% Find Brain Mask
B = squeeze(mean(V(1:10,:,:),1));
Mask = pct_brainMask(B,0,120);

% Add correlated Gaussian (spectral) noise to simulate low-dose
Vn = pct_noise(permute(V,[2 3 1]),acf,sigma,'s',Mask);
Vn = permute(Vn,[3 1 2]);

% PCT preprocess
% Whole Brain CBF
[C,AIF0] = pct_preprocess(V, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noiseless data
[Cn,AIF] = pct_preprocess(Vn, PRE, POST, AIFx,AIFy,VOFx,VOFy); % noisy data

% Crop data to Region of Interest
if exist('x0','var')
    C = C(:,y0:y0+h-1,x0:x0+w-1);
    Cn = Cn(:,y0:y0+h-1,x0:x0+w-1);
    Mask = Mask(y0:y0+h-1,x0:x0+w-1);
end

% Find the minimum bounding box
bb = pct_minBoundingBox(Mask);
C = C(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Cn = Cn(:,bb(1,1):bb(1,2),bb(2,1):bb(2,3));
Mask = Mask(bb(1,1):bb(1,2),bb(2,1):bb(2,3));

clear V Vn VOF PRE POST B acf

% Method 0: block-circulant SVD (bSVD) no noise as reference
RIF0_bSVD = pct_bsvd(C,AIF0,1,lambda,m,Mask);
CBF0_bSVD = pct_cbf(RIF0_bSVD,rho);
CBV0_bSVD = pct_cbv(RIF0_bSVD,rho);
MTT0_bSVD = pct_mtt(RIF0_bSVD);
% Mask with vasuclar elimination
MaskV = pct_velim(CBV0_bSVD,Mask);
% CBF0_bSVD(~Mask)=0; CBV0_bSVD(~Mask) = 0; MTT0_bSVD(~Mask)=0;

% Method 1: Standard Truncated Singular Value Decomposition (sSVD)
tic;
RIF_sSVD = pct_ssvd(Cn,AIF,rt,lambda,Mask);
CBF_sSVD = pct_cbf(RIF_sSVD,rho);
CBV_sSVD = pct_cbv(RIF_sSVD,rho);
MTT_sSVD = pct_mtt(RIF_sSVD);
t_sSVD = toc;

% Method 2: block-circulant SVD (bSVD)
tic;
RIF_bSVD = pct_bsvd(Cn,AIF,rt,lambda,m,Mask);
CBF_bSVD = pct_cbf(RIF_bSVD,rho);
CBV_bSVD = pct_cbv(RIF_bSVD,rho);
MTT_bSVD = pct_mtt(RIF_bSVD);
t_bSVD = toc;

% Method 3: Tikhonov
tic;
RIF_tikh = pct_tikh(Cn,AIF,rt,lambda,1,Mask);
CBF_tikh = pct_cbf(RIF_tikh,rho);
CBV_tikh = pct_cbv(RIF_tikh,rho);
MTT_tikh = pct_mtt(RIF_tikh);
t_tikh = toc;


%% Method 4: Total Variation Regularization (TTV)
% parameters
[T,h,w]=size(Cn);
% input.reg = 1e-4; % when rt=1
input.reg = [1e-4 1] * 1e-4;
input.maxitr=500;
input.l=-inf; input.u=inf;
input.no = 5; % max number of iterations allowed
[Ca,Cc] = pct_circ(AIF,Cn,[],m); % block-circulant version
Tc = T*m;
input.A=Ca;
input.b=reshape(Cc,Tc,[]);
input.tt = Tn;
input.n1=h; input.n2=w; input.nt = Tc;


% TTV
tic;
out = TTV_FCSA(input);
t_TTV = toc;
RIF_TTV = reshape(out.y,[Tc,h,w]);
RIF_TTV = RIF_TTV(1:T,:,:); % keep only first T time points
CBF_TTV = pct_cbf(RIF_TTV,rho); CBF_TTV(~Mask) = 0;
CBV_TTV = pct_cbv(RIF_TTV,rho); CBV_TTV(~Mask) = 0;
MTT_TTV = pct_mtt(RIF_TTV); MTT_TTV(~Mask) = 0;
fprintf(1,'TTV: Iter_time = %.2fsec, Iteration Num = %d\n',out.xtime(end), input.no);


%% Plot results

CBF = [CBF0_bSVD CBF_sSVD CBF_bSVD CBF_tikh CBF_TTV];
% CBV = [CBV0_bSVD CBV_sSVD CBV_bSVD CBV_tikh CBV_TTV];
% MTT = [MTT0_bSVD MTT_sSVD MTT_bSVD MTT_tikh MTT_TTV];

figure; ctshow(CBF,[],[0 80]); colorbar;

% figure;
% set(gca,'FontSize',20);
% subplot(311); ctshow(CBF,[],[0 80]); colorbar;
% subplot(312); ctshow(CBV,[],[0 15]); colorbar;
% subplot(313); ctshow(MTT,[],[0 20]); colorbar;

figure; plot(out.funv); title('TTV funv');

% RMSE and Lin's Concordance
rmse_CBF_sSVD = pct_rmse(CBF_sSVD(MaskV),CBF0_bSVD(MaskV));
rmse_CBF_bSVD = pct_rmse(CBF_bSVD(MaskV),CBF0_bSVD(MaskV));
rmse_CBF_tikh = pct_rmse(CBF_tikh(MaskV),CBF0_bSVD(MaskV));
rmse_CBF_TTV = pct_rmse(CBF_TTV(MaskV),CBF0_bSVD(MaskV));

rmse_CBV_sSVD = pct_rmse(CBV_sSVD(MaskV),CBV0_bSVD(MaskV));
rmse_CBV_bSVD = pct_rmse(CBV_bSVD(MaskV),CBV0_bSVD(MaskV));
rmse_CBV_tikh = pct_rmse(CBV_tikh(MaskV),CBV0_bSVD(MaskV));
rmse_CBV_TTV = pct_rmse(CBV_TTV(MaskV),CBV0_bSVD(MaskV));

rmse_MTT_sSVD = pct_rmse(MTT_sSVD(MaskV),MTT0_bSVD(MaskV));
rmse_MTT_bSVD = pct_rmse(MTT_bSVD(MaskV),MTT0_bSVD(MaskV));
rmse_MTT_tikh = pct_rmse(MTT_tikh(MaskV),MTT0_bSVD(MaskV));
rmse_MTT_TTV = pct_rmse(MTT_TTV(MaskV),MTT0_bSVD(MaskV));

[lin_CBF_sSVD,ci_CBF_sSVD] = pct_lincon(CBF_sSVD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_bSVD,ci_CBF_bSVD] = pct_lincon(CBF_bSVD(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_tikh,ci_CBF_tikh] = pct_lincon(CBF_tikh(MaskV),CBF0_bSVD(MaskV));
[lin_CBF_TTV,ci_CBF_TTV] = pct_lincon(CBF_TTV(MaskV),CBF0_bSVD(MaskV));

[lin_CBV_sSVD,ci_CBV_sSVD] = pct_lincon(CBV_sSVD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_bSVD,ci_CBV_bSVD] = pct_lincon(CBV_bSVD(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_tikh,ci_CBV_tikh] = pct_lincon(CBV_tikh(MaskV),CBV0_bSVD(MaskV));
[lin_CBV_TTV,ci_CBV_TTV] = pct_lincon(CBV_TTV(MaskV),CBV0_bSVD(MaskV));

[lin_MTT_sSVD,ci_MTT_sSVD] = pct_lincon(MTT_sSVD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_bSVD,ci_MTT_bSVD] = pct_lincon(MTT_bSVD(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_tikh,ci_MTT_tikh] = pct_lincon(MTT_tikh(MaskV),MTT0_bSVD(MaskV));
[lin_MTT_TTV,ci_MTT_TTV] = pct_lincon(MTT_TTV(MaskV),MTT0_bSVD(MaskV));

P_sSVD = polyfit(CBF_sSVD,CBF0_bSVD,1);
P_bSVD = polyfit(CBF_bSVD,CBF0_bSVD,1);
P_tikh = polyfit(CBF_tikh,CBF0_bSVD,1);
P_TTV = polyfit(CBF_TTV,CBF0_bSVD,1);

t_all = toc(tid_all);
