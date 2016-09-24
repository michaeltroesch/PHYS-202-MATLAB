function [i,h]=m_rule(f,a,b,n)
h=(b-a)/n;
x=linspace(a,b,n+1);
xeval=linspace(a+h/2,b-h/2,n);
fn=feval(f,xeval);
i=0;
for j=1:n
    i=i+fn(j)*(x(j+1) - x(j));
end
end