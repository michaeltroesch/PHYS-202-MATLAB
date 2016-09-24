function [x,y]=enxplot()
x=0:15;
y=[];
k=1;
while length(y) < length(x)
   y(k) = enxsum(k-1);
   k=k+1;
end
plot(x,y,'b--o')
title('The geometric series of e','FontSize',20)
ylabel('f(x) (no units)','FontSize',16)
xlabel('n (no units)','FontSize',16)
end