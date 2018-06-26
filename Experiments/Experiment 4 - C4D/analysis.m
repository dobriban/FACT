cd('C:\Dropbox\Projects\Closure\Experiments\Experiment 4 - C4D');
addpath '../../Code'
%% Analysis
data = 'C4D-CAD.C4D.Consortium.Nat.Genet.2011';
load(['../../Data/Processed/' data '.mat'],'P','SNP','chr','pos');
P = P(chr==9);
SNP = SNP(chr==9);
pos = pos(chr==9);
%
pos_CDKNA = 21967752;
b = pos_CDKNA+10e5*[-1,1];
ind = (pos >= b(1))&(pos <= b(2));
P = P(ind);
SNP = SNP(ind);
pos = pos(ind);
J = length(P);

[n_sig_s,~,p_thresh_s] = fct(P,'simes');
s = 0.1*J;
[n_sig_m,~,p_thresh_m] = fct(P,'mix-simes-hc',s,J);

%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
Z = -norminv(P);
figure, hold on
h1 = plot(pos,Z,'linewidth',1,'color',rand(1,3));
set(h1,'LineStyle',a{1});

pa= -norminv(p_thresh_s); %your point goes here
y=[pa,pa];
x=get(gca,'Xlim');
h2 = plot(x,y,'linewidth',1.5,'color',rand(1,3));
set(h2,'LineStyle',a{2});
pa= -norminv(p_thresh_m); %your point goes here
y=[pa,pa];
x=get(gca,'Xlim');
h3 = plot(x,y,'linewidth',1.5,'color',rand(1,3));
set(h3,'LineStyle',a{3});

xlim([min(pos), max(pos)]);
ylim([min(Z), max(Z)]);

xlabel('SNP Position')
ylabel('Z-score')
set(gca,'fontsize',20)
legend([h2,h3],{'Simes','Simes-HC'},'location','Best')

if savefigs==1
    filename = sprintf( './C4D.png');
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end
