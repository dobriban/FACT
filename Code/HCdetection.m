function [H, stat] = HCdetection(p, alpha, pvalcut)
% The function HCdetection realizes the Higher Criticism Detection method 
% to detect whether there is signal in high dimensional data, when the 
% signal is weak and sparse. The procedure is given by Donoho and Jin (2004)
% 
% Function: [H, stat] = HCdetection(p, alpha, pvalcut)
% 
% Inputs:
%   p        n-by-1 array of p-values from data
%   alpha    0 < alpha < 1, indicates the smallest alpha*n p-values will be
%             used to calculate the HC statistic
%            default 1/2 
%   pvalcut  0 < pvalcut < 1, indicates those small p-values (smaller than
%             pvalcut) will be taken away to avoid heavy tails of test
%             statistic
%            default 1/n
% 
% Outputs
%   H        "0" or "1" scalar, indicates whether H_0 is accepted ("0") or
%            rejected ("1").
%   stats    3-by-1 struct, including the test statistic
%      stats.HCT        Higher Criticism test statistic
%      stats.numselect  Number of variables at which HCT is achieved
%      stats.HC         HC score for each feature
% 
% 
% Example:
%  n = 1e+6;
%  eps = 1/1000;
%  mu = sqrt(2*0.15*log(n));
% 
%  data = randn(n, 1);
%  p = 2*(1 - normcdf(abs(data)));
%  [H01, stat01] = HCdetection(p, 1, 0);
%  [H02, stat02] = HCdetection(p, 1/2, 1/n);
% 
%  rng(0,'twister');
%  r = randi([1 1e+6],1,1000);
%  data(r) = randn(1000, 1) + mu;
%  p = 2*(1 - normcdf(abs(data)));
%  [H11, stat11] = HCdetection(p, 1, 0);
%  [H12, stat12] = HCdetection(p, 1/2, 1/n);
% 
% Reference: 
% Donoho and Jin (2004) Higher criticism for detecting sparse heterogeneous 
% mixtures




% Error checking

if (nargin<1 || isempty(p))
  error 'Please input the p-value array'
end

if (max(p) > 1 || min(p) < 0)
  error 'The p-value array has unexpected values'
end

n = length(p);
if (nargin<2||isempty(alpha))
    alpha = 1/2;
end
if (alpha > 1 || alpha < 0)
  error 'alpha must be between 0 and 1'
end

if (nargin<3||isempty(pvalcut))
    pvalcut = 1/n;
end
if (pvalcut > 1 || pvalcut < 0)
  error 'pvalcut must be between 0 and 1'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %Find the HC threshold with input p-value
kk    = (1:n)'/(1 + n);
psort = sort(p);
HC  = sqrt(n)*(kk - psort)./sqrt(psort - psort.^2);
HC(psort < pvalcut)  = -inf;
HC(round(alpha*n)+1:n) = -inf;
HCT = max(HC);
numselect   = find(HC == HCT);

stat.HCT = HCT; stat.numselect = numselect; stat.HC = HC;
%T = 1;
T = 1;
H = (HCT/sqrt(2*log(log(n)))/(1 + 1/sqrt(log(log(n)))) > T)*1;
end