fid=fopen('expsin.dat','w');   %Open file: expsin.dat in cd for write
x=linspace(0,4*pi,400);        %Define linear space with 400 values from 0 to 4pi
f_x=3*exp((-x.^2)/10).*sin(x); %our function of x (Force in Newtons)            
for i=1:length(x)
fprintf(fid,'%0.5f        %0.5f\n',x(i),f_x(i));
end
fclose(fid);
load('expsin.dat');
plot(expsin(:,1),expsin(:,2),'k');
xlabel('x (meters)','FontSize',16);
ylabel('F (Newtons)','FontSize',16);
title('F(x) = 3*exp((-x^2)/10)*sin(x)','FontSize',18);

