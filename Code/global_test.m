function R = global_test(Q,type,varargin)
% Global Testing method
%type:
% 'bonf' - Bonferroni
% 'fisher' - Fisher
% 'stouffer' - Stouffer
% 'wilkinson' - Wilkinson
% etc...
%Q is assumed sorted

alpha = 0.05;
m = length(Q);
switch type
    case 'bonf'
        T = Q(1);
        c = alpha/m;
        
    case 'fisher'
        T = 2*sum(log(Q));
        c = -chi2inv(1-alpha,2*m);
        
    case 'stouffer'
        T = sum(norminv(Q));
        c = m^(1/2)*norminv(alpha);
        
    case 'wilkinson'
        d = alpha;
        T = -sum(Q <= d);
        c = -binoinv(1-alpha,m,d);
        
    case 'simes'
        A = (1:m)';
        T = min(Q./A);
        c = alpha/m;
        
    case 'HC'
        th = 0; %1/n
        al = 1; %1/2
        T = -HCdetection(Q, al, th);
        c = -1;
        
    case 'mix-simes-fisher'
        args = varargin{1};
        s = args{1};
        n = args{2};
        th = n-s+1;
        if (m <= th) %the number of tests is small, so can only use one significant p-value
            A = (1:m)';
            T = min(Q./A);
            c = alpha/m;
        else  %the number of tests is large, so can use many significant p-values
            T = 2*sum(log(Q));
            c = -chi2inv(1-alpha,2*m);
        end
        
    case 'mix-simes-hc'
        args = varargin{1};
        s = args{1};
        n = args{2};
        th = n-s+1;
        if (m <= th) %the number of tests is small, so can only use one significant p-value
            A = (1:m)';
            T = min(Q./A);
            c = alpha/m;
        else  %the number of tests is large, so can use many significant p-values
            th = 0; %1/n
            al = 1; %1/2
            T = -HCdetection(Q, al, th);
            c = -1;
        end
end

if T <=c
    R = 1;
else
    R = 0;
end