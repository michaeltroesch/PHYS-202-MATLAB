function [t,y]=heun_method(dydt,tinit,yinit,tfinal,n)
dt=(tfinal-tinit)/n;
t=[tinit zeros(1,n)];
y=[yinit zeros(1,n)];
for j=1:n
    t(j+1)=t(j)+dt;
    ynew=y(j)+dt*dydt(t(j),y(j));
    y(j+1)=y(j)+(dt/2)*(dydt(t(j),y(j))+dydt(t(j),ynew));
end
end
