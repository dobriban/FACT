%% Detect Sparse Heterogeneous Mixtures with Higher Criticism Statistic
% To detect sparse and weak signals in mixtures, Donoho and Jin (2004)
% proposed that Higher Criticism statistic can be used, and has optimal
% behavior. The function HCdetection is to realize this procedure. It is
% assumed that the size alpha_n slowly converges to 0 for this test.


%% An Example with Simulated Data
% Generate the data with large sample size n, and the same distribution
n = 1e+6;
data = randn(n, 1);

%%
% Calculate two-sided p-value for the data
p = 2*(1 - normcdf(abs(data)));

%%
% Use the p-value to do test, the first one is classicial Higher Criticism
% test, and the second one is a refinement
[H01, stat01] = HCdetection(p, 1, 0);
[H02, stat02] = HCdetection(p, 1/2, 1/n);
H01
H02
%% 
% The result shows that the null hypothesis is accepted.

%%
% In the alternative case, 
% Generate the data such that with 1/1000 probability, the sample comes
% from normal distribution with a small mean.
eps = 1/1000;
mu = sqrt(2*0.15*log(n));
r = binornd(1, eps, n, 1);
data(r == 1) = data(r == 1) + mu;

%%
% Calculate two-sided p-value for the data
p = 2*(1 - normcdf(abs(data)));

%%
% Use the p-value to do test, the first one is classicial Higher Criticism
% test, and the second one is a refinement
[H11, stat11] = HCdetection(p, 1, 0);
[H12, stat12] = HCdetection(p, 1/2, 1/n);
H11
H12
%%
% For both statistics, the null hypothesis is rejected. 