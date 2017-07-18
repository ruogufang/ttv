function [cbf,cbv,mtt,ttp,params] = spd_map(params)
%SPD_MAP optimizes 2-D signals of blood flow map using Maximum A Posterior.
%Steepest descent optimiation is used.
%
% Optimization of convolution model with dictionary induced prior using
% steepest descent
%
% INPUT:
%       PARAMS - structure of parameters (for detailed field definition,
%       see spd.m)
%
% OUTPUT:
%       MAPS.
%       CBF      - Estimated CBF map [H x W]
%       MTT      - Estimated MTT map [H x W]
%       CBV      - Estimated CBV map [H x W]
%       TTP      - Estimated TTP map [H x W] (optional)
%       PARAMS - structure of parameters with fields f,J,prior,nz updated
%
%
% SPD_MAP calls two functions: SPD_PRIOR and SPD_SD for prior estimation
% and steepest descent optimiztion
%
%
%   Ruogu Fang 4/10/2014
%   Advanced Multimedia Processing (AMP) Lab
%   Department of Electrical and Computer Engineering
%   Cornell University

cbf_prior = params.cbf0; % copy the noisy cbf map
mtt_prior = params.mtt0; % copy of the noisy mtt map
cbv_prior = params.cbv0; % copy of the noisy cbv map
if isfield(params,'ttp0')
    ttp_prior = params.ttp0; % copy of the noisy cbv map
end

% Iteration of two steps
for i = 1 : params.round
    % CBF
    %     params.lambda = 10; % weight of the sparsity term in SPAMS package
    params.f = cbf_prior;
    params = spd_prior(params);
    cbf_prior = params.prior;
    params = spd_sd(params);
    % MTT
    %     params.lambda = 1; % weight of the sparsity term in SPAMS package
    params.f = mtt_prior;
    params = spd_prior(params);
    mtt_prior = params.prior;
    % CBV
    params.f = cbv_prior;
    params = spd_prior(params);
    cbv_prior = params.prior;
    % TTP
    if isfield(params,'ttp0')
        params.f = ttp_prior;
        params = spd_prior(params);
        ttp_prior = params.prior;
    end
    
end

% Output of restored perfusion map
params = spd_prior(params);
cbf = cbf_prior;
mtt = mtt_prior;
cbv = cbv_prior;
if isfield(params,'ttp0')
    ttp = ttp_prior;
end

end