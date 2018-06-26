function [n_sig,ind_rej,which_rej,p_thresh] = fct(P,type,varargin)
% Fast Closed Testing method
%type:
% 'bonf' - Bonferroni
% 'fisher' - Fisher
% 'stouffer' - Stouffer
% 'wilkinson' - Wilkinson
% Etc...

[Q,ind]  = sort(P);
m = length(Q);
for k=1:m+1
    for j=k:m
        index = [k,(j+1):m];
        T = global_test(Q(index),type,varargin);
        if (T == 0)
            break
        end
    end
    if (T == 0)
        break
    end
end

ind_rej = ind(1:k-1);
which_rej = zeros(m,1);
which_rej(ind_rej)=1;
n_sig = k-1;
if k>2
    p_thresh = Q(k-1);
else
    p_thresh=0;
end



