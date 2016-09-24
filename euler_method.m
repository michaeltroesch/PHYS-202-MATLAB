function [t,y]=euler_method(dydt,tinit,yinit,tfinal,n)
dt=(tfinal-tinit)/n;
t=[tinit zeros(1,n)];
y=[yinit zeros(1,n)];
for i=1:n
    t(i+1)=t(i)+dt;
    y(i+1)=y(i)+dt*dydt(t(i),y(i));
end
end