%Experiment with closed testing
cd('C:\Dropbox\Projects\Closure\Experiments\Experiment 3 - Mixture')
addpath('C:\Dropbox\Projects\Closure\Code')
beep off
%% Set parameters
n_mc = 1e2;
n = 100;
sigma2_global = 2*n;
%sparsity = [1,10,100];
sparsity = [10,30,50,70,90,100];
l = length(sparsity);
effect_size_arr = [0,1,2];
%effect_size_arr = 0;
m = length(effect_size_arr);

for j=1:m
    effect_size = effect_size_arr(j);
    mean_global = effect_size*sigma2_global^(1/2);
    %%
    rng(2);
    mean_n_sig_f = zeros(l,1);
    mean_n_sig_s = zeros(l,1);
    mean_n_sig_w = zeros(l,1);
    fwer_f = zeros(l,1);
    fwer_s = zeros(l,1);
    fwer_w = zeros(l,1);
    for k=1:l
        s = sparsity(k);
        mu = (mean_global^2/s)^(1/2);
        mu_v = zeros(n,1);
        mu_v(1:s) = mu;
        ti = tic;
        n_sig_f = zeros(n_mc,1);
        n_sig_s = zeros(n_mc,1);
        n_sig_w = zeros(n_mc,1);
        for i=1:n_mc
            X = mu_v+randn(n,1);
            P = 1-normcdf(X);
            
            n_sig_f(i) = fct(P,'fisher');
            n_sig_s(i) = fct(P,'simes');
            n_sig_w(i) = fct(P,'mix-simes-hc',s,n);
        end
        toc(ti);
        mean_n_sig_f(k) = mean(n_sig_f(1:s));
        mean_n_sig_s(k) = mean(n_sig_s(1:s));
        mean_n_sig_w(k) = mean(n_sig_w(1:s));
        if (j==1)
            fwer_f(k) = mean(n_sig_f>0);
            fwer_s(k) = mean(n_sig_s>0);
            fwer_w(k) = mean(n_sig_w>0);
        end
    end
    %%
    rng(2);
    savefigs =1; a = {'-','--','-.',':'};
    figure, hold on
    if j>1 %power
        h1 = plot(sparsity,mean_n_sig_f,'linewidth',3,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        h2 = plot(sparsity,mean_n_sig_s,'linewidth',3,'color',rand(1,3));
        set(h2,'LineStyle',a{2});
        h3 = plot(sparsity,mean_n_sig_w,'linewidth',3,'color',rand(1,3));
        set(h3,'LineStyle',a{3});
        ylabel('# True Discoveries')
    else %fwer
        h1 = plot(sparsity,fwer_f,'linewidth',3,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        h2 = plot(sparsity,fwer_s,'linewidth',3,'color',rand(1,3));
        set(h2,'LineStyle',a{2});
        h3 = plot(sparsity,fwer_w,'linewidth',3,'color',rand(1,3));
        set(h3,'LineStyle',a{3});
        ylabel('FWER')
    end
    
    
    xlabel('Sparsity')
    xticks(sparsity)
    set(gca,'fontsize',20)
    xlim([min(sparsity), max(sparsity)]);
    
    legend([h1,h2,h3],{'Fisher','Simes','Simes-HC'},'location','Best')
    
    if savefigs==1
        filename = sprintf( './FCT-mix-n=%d-effect-size=%.2f.png',n,effect_size);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        %close(gcf)
    end
end
