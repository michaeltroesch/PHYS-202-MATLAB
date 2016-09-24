function [t,y]=rk4(dydt,tn,dt,y0)
t=(0:dt:tn)';
n=length(t);
y=y0*ones(n,1);
dt2=dt/2; dt3=dt/3; dt6=dt/6;
for j=2:n
    k1=feval(dydt,t(j-1),y(j-1));
    k2=feval(dydt,t(j-1)+dt2,y(j-1)+dt2*k1);
    k3=feval(dydt,t(j-1)+dt2,y(j-1)+dt2*k2);
    k4=feval(dydt,t(j-1)+dt,y(j-1)+dt*k3);
    y(j)=y(j-1)+dt6*(k1+k4)+dt3*(k2+k3);
end
end