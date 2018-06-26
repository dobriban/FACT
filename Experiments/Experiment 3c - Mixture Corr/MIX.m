%Experiment with closed testing
cd('C:\Dropbox\Projects\Closure\Experiments\Experiment 3c - Mixture Corr')
addpath('C:\Dropbox\Projects\Closure\Code')
beep off
%% Set parameters
n_mc = 1e3;
n = 100;
sigma2_global = 2*n;

sparsity = [0,50,90];
l = length(sparsity);

effect_size_arr = 1;
m = length(effect_size_arr);

for corr_type=1:2
%corr_type=2;
switch corr_type
    case 1
        corr = [-0.99/(n-1),linspace(0,0.9,10)];
    case 2
        corr = linspace(-0.9,0.9,10);
end
p = length(corr);
%%

for j=1:m
    effect_size = effect_size_arr(j);
    mean_global = effect_size*sigma2_global^(1/2);
    %%
    rng(2);
    mean_n_sig_hc = zeros(l,p);
    fwer = zeros(l,p);
    for k=1:l
        s = sparsity(k);
        mu = 0;
        if s>0
            mu = (mean_global^2/s)^(1/2);
        end
        mu_v = zeros(n,1);
        mu_v(1:s) = mu;
        ti = tic;
        for q=1:p
            co = corr(q);
            n_sig_hc = zeros(n_mc,1);
            
            switch corr_type
                case 1
                    Sigma = (1-co)*eye(n)+co*ones(n);
                case 2
                    a = co.^(0:n-1);
                    Sigma = toeplitz(a);
            end
            Sigma12 = Sigma^(1/2);
            
            for i=1:n_mc
                X = mu_v+Sigma12*randn(n,1);
                P = 1-normcdf(X);
                
                [~,~,which_rej] = fct(P,'mix-simes-hc',s,n);
                n_sig_hc(i) = sum(which_rej(s+1:n));
            end
            toc(ti);
            mean_n_sig_hc(k,q) = mean(n_sig_hc);
            fwer(k,q) = mean(n_sig_hc>0);
        end
    end
    %%
    rng(2);
    savefigs =1; a = {'-','--','-.',':'};
    figure, hold on
    h1 = plot(corr,fwer,'linewidth',3);
    set(h1,'LineStyle',a{1});
    xlabel('Correlation')
    ylabel('FWER')
    set(gca,'fontsize',20)
    xlim([min(corr), max(corr)]);
    xticks(linspace(min(corr), max(corr),5));
    
    legendCell = cellstr(num2str(sparsity', 's=%-d'));
    legend(h1,legendCell,'location','Best')
    
    if savefigs==1
        filename = sprintf( './FCT-corr-n=%d-effect-size=%.2f-corr-type=%d.png',n,effect_size,corr_type);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        %close(gcf)
    end
end
end