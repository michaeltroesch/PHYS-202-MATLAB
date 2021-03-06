function X2=chisqr(x,f,error)
N=length(x);
w=1./(error.^2);
del=sum(w)*sum(w.*(x.^2))-(sum(w.*x))^2;
a=(sum(w.*(x.^2))*sum(w.*f)-sum(w.*x)*sum(w.*x.*f))/del;
b=(sum(w)*sum(w.*x.*f)-sum(w.*x)*sum(w.*f))/del;
chisqr_dof=(1/(N-2))*sum(w.*(f - a - b*x).^2);
sig_a=sqrt((sum(w.*(x.^2)))/del);
sig_b=sqrt(sum(w)/del);
X2=chisqr_dof;