function [params] = spd_prior(params)
%SPD_PRIOR estimates the prior blood flow map from the learned dictionary
%
% INPUT:
%       PARAMS - structure of parameters (for detailed field definition,
%       see spd.m)
%
% OUTPUT:
%       PARAMS - structure of parameters with fields prior, nz updated
%
% SPD_PRIOR function is called by SPD_MAP
%
%   Ruogu Fang 6/26/2012
%   Advanced Multimedia Processing (AMP) Lab
%   Department of Electrical and Computer Engineering
%   Cornell University
%
% Code derived from ksvddenoise.m of
%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs


% parse input arguments %
f = params.f;
D = params.dict;
blocksize = params.blocksize;

% blocksize %
if (numel(blocksize)==1)
    blocksize = ones(1,2)*blocksize;
end

% stepsize %
if (isfield(params,'stepsize'))
    stepsize = params.stepsize;
    if (numel(stepsize)==1)
        stepsize = ones(1,2)*stepsize;
    end
else
    stepsize = ones(1,2);
end
if (any(stepsize<1))
    error('Invalid step size.');
end


% lambda %
if (~isfield(params,'lambda'))
    params.lambda = params.maxval/(10*params.sigma);
end

% lambda %
if (~isfield(params,'dt'))
    params.dt = 1;
end


% compute G %
G = D'*D;


% verify dictionary normalization %
atomnorms = diag(G);
if (any(abs(atomnorms-1) > 1e-2))
    error('Dictionary columns must be normalized to unit length');
end


% denoise the signal %

nz = 0;  % count non-zeros in block representations
k = 0; % counter

% the denoised signal
y = zeros(size(f));

blocknum = prod(floor((size(f)-blocksize)./stepsize) + 1);

for j = 1:stepsize(2):size(y,2)-blocksize(2)+1
    k = k+1;
    % the current batch of blocks
    blocks = im2colstep(f(:,j:j+blocksize(2)-1),blocksize,stepsize);
    
    % remove DC
    [blocks, dc] = remove_dc(blocks,'columns');
    
    gamma = zeros(size(D,2),size(blocks,2));
    % denoise the blocks using Lasso
    if any(blocks(:))
        gamma = mexLasso(blocks,D,params);
    else
        fprintf('Empty block\n');
    end
    nz = nz + nnz(gamma);
    cleanblocks = add_dc(D*gamma, dc, 'columns');
    
    cleanim = col2imstep(cleanblocks,[size(y,1) blocksize(2)],blocksize,stepsize);
    y(:,j:j+blocksize(2)-1) = y(:,j:j+blocksize(2)-1) + cleanim;
    
    %     printf('%d/%d blocks processed.\n',size(blocks,2)*k, blocknum);
    
end

params.nz = nz/blocknum;  % average number of non-zeros

% average the denoised and noisy signals
cnt = countcover(size(f),blocksize,stepsize);
y = (y+params.mixture*f)./(cnt + params.mixture);
% y(y<0)=0;
params.prior = y;

end

