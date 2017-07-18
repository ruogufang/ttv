function params = spd_sd(params)
%SPD_SD estimates the blood flow map using a steepest descent algorithm
%
% INPUT:
%       PARAMS - structure of parameters (for detailed field definition,
%       see spd.m)
%
% OUTPUT:
%       PARAMS - structure of parameters with fields f, J updated
%
% SPD_SD function is called by SPD_MAP
%
%   Ruogu Fang 6/26/2012
%   Advanced Multimedia Processing (AMP) Lab
%   Department of Electrical and Computer Engineering
%   Cornell University

if ~isfield(params,'dt')
    dt = 1;
else
    dt=params.dt;
end

tid = tic;

% parse input arguments %
f = reshape(params.f,[],1); % vecorize f
prior = reshape(params.prior,[],1);
beta = params.beta;
rho = params.rho;
m = params.m;
C = params.C; % C is the tissue concentraiton curves of N voxels [TxHxW]
AIF = params.AIF; % AIF is the arterial input function [Tx1]
[T,H,W] = size(C);
R = params.R ./ repmat(max(params.R),[T 1 1]); % normalize R !!! zero divider !!!
R(isnan(R)) = 0;
R = reshape(R,T,[]);

% Block-circulant
L = m * T;
AIF = cat(1,AIF,zeros(L-T,1));
Ca = pct_bcaif(AIF,1);
C = cat(1,C,zeros(L-T,H,W));
R = cat(1,R,zeros(L-T,H*W));
K = Ca * R; % convolution of Ca with residue R 
C = C/dt; % correct \delta t for sampling rate

% Solve minimize J = beta*||f-prior(f)||_2^2+||C-K*diag(f)||_2^2

% Step 1: Data preprocessing
% Sparse version of b. Case: p=2
b = sparse([reshape(C,[],1);beta.*prior]); % Vectorize Y
% Sparse verion of A
idx_y = (1:(L+1)*(H*W))';
idx_x = [reshape(repmat(1:H*W,L,1),[],1);(1:H*W)'];
s = [K(:);beta.*ones(H*W,1)];
A = sparse(idx_y,idx_x,s);


% Step 2: Steepest descent optimization
e = sparse(A'*(b-A*f));
step=(e'*e)/((A*e)'*(A*e));

% Convert from mL/g/s to mL/100g/min
scaling_factor = 60 * 100 / rho;
step = step/scaling_factor;
f_new = f + step*e;
params.f = reshape(f_new,size(params.f));

% params.f = params.f * scaling_factor; % convert parameters from mL/g/s to mL/100g/min
J_new = norm(b-A*f_new,2)^2+beta*norm(f_new-f,2)^2;
params.J = [params.J J_new];

t=toc(tid);
fprintf('Optimization time:%.2f\t Cost function:%.2e\n',t,J_new);

end
