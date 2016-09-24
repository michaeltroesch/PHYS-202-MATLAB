%% Problem 1.
% 1a.  Charging the capacitor in an RC circuit, voltage across
% the capacitor is measured with time t.
t =[0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
v = [0, 5.85, 10.3, 11.9, 14.5, 15, 16, 16.2, 16.3, 16.5, 16.55];
tq=0:0.01:30;
v_nearest=interp1(t,v,tq,'nearest');
plot(t,v,'b.',tq,v_nearest,'c-','MarkerSize',14);
title('Comparison of Interpolation Techniques','FontSize',16)
ylabel('Capacitor Potential (V)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
v_linear=interp1(t,v,tq,'linear');
plot(tq,v_linear,'m-','MarkerSize',14)
v_pchip=interp1(t,v,tq,'pchip');
plot(tq,v_pchip,'y-','MarkerSize',14)
v_cubicspline=interp1(t,v,tq,'spline');
plot(tq,v_cubicspline,'k-','MarkerSize',14)
legend({'Data','Nearest','Linear','PCHIP','Cubic Spline'},'Location','NorthWest');
str1=strcat('Nearest Neighbor=',num2str(v_nearest(751)));
str2=strcat('Linear=',num2str(v_linear(751)));
str3=strcat('Piecewise Hermite Cubic=',num2str(v_pchip(751)));
str4=strcat('Cubic Spline=',num2str(v_cubicspline(751)));
text(12,11,'Values for t=7.5 s by differing methods:')
text(12,10,str1)
text(12,9,str2)
text(12,8,str3)
text(12,7,str4)
hold off
%%
% 1b. Same as part (a), but with different data.
t =[0 .1 .499 .5 .6 1.0 1.4 1.5 1.899 1.9 2.0];
v = [0 .06 .17 .19 .21 .26 .29 .29 .30 .31 .31];
tq=0:0.01:2;
v_nearest=interp1(t,v,tq,'nearest');
plot(t,v,'b.',tq,v_nearest,'b-','MarkerSize',14);
title('Comparison of Interpolation Techniques','FontSize',16)
ylabel('Capacitor Potential (V)','FontSize',14)
xlabel('Time (s)','FontSize',14)
hold on
v_linear=interp1(t,v,tq,'linear');
plot(tq,v_linear,'r-','MarkerSize',14)
v_pchip=interp1(t,v,tq,'pchip');
plot(tq,v_pchip,'c-','MarkerSize',14)
v_cubicspline=interp1(t,v,tq,'spline');
plot(tq,v_cubicspline,'k-','MarkerSize',14)
legend({'Data','Nearest','Linear','PCHIP','Cubic Spline'},'Location','SouthEast');
hold off
%% Commentary on Graph (1b)
% Looking at the graph above, we can see that the most off-base
% interpolation comes from the cubic spline method--we had data points that
% were in some cases very close together and matching the second
% derivatives in those cases sent our interpolated line out on wild
% tangents!  I might consider using a linear method on this data for that
% reason.  PCHIP has less issues than Cubic Spline, but it still shows some
% errors I'm somewhat uncomfortable with.  Linear also has the advantage of
% computing speed.  I could trust PCHIP and Linear in this case however.
% Nearest Neighbor always has some big problems, and Cubic Spline doesn't
% seem to handle this sort of unevenly spaced data very well at all.
%% Problem 2.
C=10:39;
R=[3.239,3.118,3.004,2.897,2.795,2.7,2.61,2.526,2.446,2.371,2.3,2.233...
2.169,2.11,2.053,2,1.95,1.902,1.857,1.815,1.774,1.736,1.7,1.666,1.634...
1.603,1.574,1.547,1.521,1.496];
Rq=3.239:-0.001:1.496;
C_nearest=interp1(R,C,Rq,'nearest');
C_linear=interp1(R,C,Rq,'linear');
C_pchip=interp1(R,C,Rq,'pchip');
C_cubicspline=interp1(R,C,Rq,'spline');
plot(R,C,'b.',Rq,C_nearest,'b-','MarkerSize',14);
title('Interpolation of Thermistor Data','FontSize',16)
ylabel('Temperature (C)','FontSize',14)
xlabel('Resistance (Megaohms)','FontSize',14)
hold on
plot(Rq,C_linear,'r-','MarkerSize',14)
plot(Rq,C_pchip,'c-','MarkerSize',14)
plot(Rq,C_cubicspline,'k-','MarkerSize',14)
legend({'Data','Nearest','Linear','PCHIP','Cubic Spline'},'Location','SouthWest');
hold off
%%
%prompt='Please input the resistance in megaohms!';
%x=input(prompt);
%y=find(x == Rq);
%if isempty(y)
%    error('Outside of interpolated data range!')
%end
%str1=strcat('Nearest Neighbor=',num2str(C_nearest(y)));
%str2=strcat('Linear=',num2str(C_linear(y)));
%str3=strcat('Piecewise Hermite Cubic=',num2str(C_pchip(y)));
%str4=strcat('Cubic Spline=',num2str(C_cubicspline(y)));
%disp(str1)
%disp(str2)
%disp(str3)
%disp(str4)
%%
Nearest=[34;27;14];
Linear=[33.5938;26.7708;14.3368];
PCHIP=[33.5885;26.7656;14.3296];
Spline=[33.5891;26.7662;14.3297];
r1='1.647 Megaohms';
r2='1.913 Megaohms';
r3='2.763 Megaohms';
Resistances={r1,r2,r3};
T=table(Nearest,Linear,PCHIP,Spline,'RowNames',Resistances);
disp(T)
%% Problem 3.
% 3a.  Importing Lineardata.dat
filename='Lineardata.dat';
headerLinesIn=2;
delimiterIn=' ';
A=importdata(filename,delimiterIn,headerLinesIn);
t=transpose(A.data(:,1));
y=transpose(A.data(:,2));
plot(t,y,'b*')
xlabel('t (seconds)','FontSize',14)
ylabel('y (meters)','FontSize',14)
title('Lineardata.dat','FontSize',16)
%%
%function yq=m_linear(x,y,xq)
%    %x,y: x and y data points to be interpolated from
%    %xq: points to evaluate via interpolation
%    %yq: returned value for line @ points given by xq  
%    d=diff(y)./diff(x);
%    n=length(x);
%    k=ones(size(xq));
%    for i=1:n-1
%        if x(i+1) > x(i)
%            k(x(i) <= xq) = i;
%        else
%            k(x(i) >= xq) = i;
%        end
%    end
%    s=xq-x(k);
%    yq=y(k) + s.*d(k);
%end
%%
tq=0:0.001:20;
yq=m_linear(t,y,tq);

% m_linear is outlined above.
% It is commented out here, because I am publishing this script with
% MATLAB's built-in publishing function.

plot(t,y,'b*',tq,yq,'k-')
xlabel('t (seconds)','FontSize',14)
ylabel('y (meters)','FontSize',14)
title('Manual Linear Interpolation','FontSize',16)
%%
%prompt='Please input the value of t for which to evaluate y(t)!';
%x=input(prompt);
%y=find(x==tq);
%if isempty(y)
%    error('Outside of interpolated data range!')
%end
%fprintf('\nTime value to evaluate via interpolation: %6.4f\n interpolated y data: %6.4f',tq(y),yq(y))
%%
t_i=[2.73;24.512;11.785;18.993];
y_i=[yq(2731);NaN;yq(11786);yq(18994)];

%These values were obtained by using the prompt script above.  It is
%commented out so I could use MATLAB's publishing function without errors.

T=table(t_i,y_i);
disp(T)
%% Problem 4.
% 4a. f(x) = x * ln(x), beginning with an comparison of numerical methods
% of estimating the value of the function's first derivative.
x=1;
for n=1:0.5:5
    dx=10^-n;
    dfdx_exact=log(x) + 1;
    dfdx_forward((n*2)-1)=(((x+dx) * log(x+dx)) - (x * log(x)))/dx;
    dfdx_back((n*2)-1)=((x * log(x)) - ((x-dx) * log(x-dx)))/dx;
    dfdx_center((n*2)-1)=(((x+dx) * log(x+dx)) - ((x-dx) * log(x-dx)))/(2*dx);
end

for n=1:0.5:5
    dx((n*2)-1)=10^-n;
end

Step_Size=transpose(dx);
Forward=transpose(dfdx_forward);
Backward=transpose(dfdx_back);
Center=transpose(dfdx_center);
T=table(Step_Size,Forward,Backward,Center);
disp(T)
%% Comments on the table (4a)
% The table above shows us that the difference between numerical
% approximation and the exact value is less than 10^-5 at different values
% of n for different methods.  For the forward method: n=4 (10^-4), backward
% method: n=5 (10^-5), center method: n=2.5 (10^-2.5)
%%
% 4b. A similar problem to part (a), but instead of first derivative's
% numerical approximation, we will be evaluating and comparing instead for
% methods of approximating second and third derivatives of the same
% function; f(x) = x * ln(x)
x=1;
for n=0.5:0.5:4
    dx=10^-n;
    d2_fx(n*2)=((x+dx)*log(x+dx)-(2*(x*log(x)))+((x-dx)*log(x-dx)))/(dx^2);
    d3_fx(n*2)=(1/2)*(((x+dx+dx)*log(x+dx+dx))-(2*((x+dx)*log(x+dx)))+(2*((x-dx)*log(x-dx)))-((x-dx-dx)*log(x-dx-dx)))/dx^3;
end


for n=0.5:0.5:4
    dx(n*2)=10^-n;
end

Step_Size=transpose(dx);
Second_Derivative=transpose(d2_fx);
Third_Derivative=transpose(d3_fx);
T=table(Step_Size,Second_Derivative,Third_Derivative);
disp(T)
%% Comments on the table (4b)
% Our new table of values shows when our numerical approximations are within
% 10^-5 of the actual exact value... for the 2nd derivative, it happens when
% n=2.  For the third derivative, it's when n=2.5.  Nice!
%% Problem 5.
% 5a.  Plotting the displacement of a spring over time.
clear all;
t=0:0.0001:20;
p=2.*exp(-t./6).*cos(3*t);
plot(t,p)
title('Position of a spring over time','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Position (m)','FontSize',14)
%%
% 5b. Calculating the exact and numerical approximations of velocity and
% acceleration.  Plotting the velocity.  A step size of 10^-3 achieves a
% residual between the approximation and exact value less than 10^-4, but I
% used 10^-4 step size, because that was required for that level of
% precision on the numerically approximated graph of acceleration.

dpdt_exact= -(cos(3*t).*exp(-t/6))./3 - 6*sin(3*t).*exp(-t/6);
d2pdt_exact=2.*sin(3.*t).*exp(-t./6) - (323.*cos(3.*t).*exp(-t./6))/18;

%Calculating the exact value of the derivative.
%Functions were determined with MATLAB symbolic toolbox.

dt=10^-4;
N=length(t);
dpdt=ones(size(dpdt_exact));
dpdt(2:N-1)=(p(3:N)-p(1:N-2))/(2*dt);
dpdt(1)=2*dpdt(2)-dpdt(3);
dpdt(N)=2*dpdt(N-1)-dpdt(N-2);
d2pdt(2:N-1)=(p(3:N)-2*p(2:N-1)+p(1:N-2))/(dt)^2;
d2pdt(1)=2*d2pdt(2)-d2pdt(3);
d2pdt(N)=2*d2pdt(N-1)-d2pdt(N-2);
residual=dpdt-dpdt_exact;
residual2=d2pdt-d2pdt_exact;
plot(t,dpdt_exact,'b*-')
hold on
plot(t,dpdt,'ro-','MarkerSize',2)
title('Exact & Numerical Approx. graphs of velocity','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Velocity (m/s)','FontSize',14)
legend({'Exact','Numerical'},'location','SouthEast','FontSize',16)
hold off
%%
plot(t,residual,'r*','MarkerSize',1)
title('Residual, approx. and exact velocity','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Numerical approx minus exact','FontSize',14)
%%
% 5c. Plotting the acceleration of the mass-spring system.  As mentioned
% before, the step size to achieve the desired accuracy of a residual no
% greater than 10^-4 at any point on the graph was 0.0001 (10^-4).
plot(t,d2pdt_exact,'b*-')
hold on
plot(t,d2pdt,'ro-','MarkerSize',2)
title('Exact & Numerical Approx. graphs of acceleration','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Acceleration (m/s^2)','FontSize',14)
legend({'Exact','Numerical'},'location','SouthEast','FontSize',16)
hold off
%%
plot(t,residual2,'r*','MarkerSize',1)
title('Residual, approx. and exact accel.','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Numerical approx minus exact','FontSize',14)
%% Problem 6.
% 6a.  A plot of the displacement above the earth's surface as a function
% of time.  It takes 122.13~ seconds to hit the Earth (there is some
% rounding error here, on the order of magnitude 10^-5.
clear all;
t=0:0.01:122.13;
h=10000;
k=0.1;
g=9.81;
v_o=100;
z=h-(g.*t)./k + ((k*v_o + g)./k^2).*(1-exp(-k.*t));
plot(t,z)
title('Displacement vs. time, v\o= 100m/s','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Displacement (m)','FontSize',14)
%%
% 6b. I do all the calculations for part (c) within this code.  This plots
% the velocity as both a numerical approximation and as an exact derivative
% function, and then we will plot the residual between them.  A step size
% of 0.01 gave us the required precision of a residual no greater than
% 10^-5 at any point over the interval.

dt=0.01;
dzdt_exact= (1981*exp(-t./10))./10 - 981/10;
d2z_exact= -(1981.*exp(-t./10))./100;

N=length(t);
dzdt=ones(size(dzdt_exact));
dzdt(2:N-1)=(z(3:N)-z(1:N-2))/(2*dt);
dzdt(1)=(3*dzdt(2))-(3*dzdt(3))+dzdt(4);
dzdt(N)=(3*dzdt(N-1))-(3*dzdt(N-2))+dzdt(N-3);
d2z(2:N-1)=(z(3:N)-2*z(2:N-1)+z(1:N-2))/(dt)^2;
d2z(1)=(3*d2z(2))-(3*d2z(3))+d2z(4);
d2z(N)=(3*d2z(N-1))-(3*d2z(N-2))+d2z(N-3);
residual=dzdt-dzdt_exact;
residual2=d2z-d2z_exact;
plot(t,dzdt_exact,'b*-')
hold on
plot(t,dzdt,'ro-','MarkerSize',2)
title('Exact & Numerical Approx. graphs of velocity','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Velocity (m/s)','FontSize',14)
legend({'Exact','Numerical'},'location','SouthEast','FontSize',16)
hold off
%used MATLAB symbolic toolbox to calculate derivative functions.
%%
plot(t,residual,'r*','MarkerSize',1)
title('Residual, approx. and exact accel.','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Numerical approx. minus exact','FontSize',14)
%%
% 6c. Plotting our previously calculated data for the acceleration (second
% derivative) of the particle falling in Earth's atmosphere.  The step
% size, as before, is 0.01 to achieve a residual no greater than 10^-5 at
% any point on the time interval.
plot(t,d2z_exact,'b*-')
hold on
plot(t,d2z,'ro-','MarkerSize',2)
title('Exact & Numerical Approx. graphs of acceleration','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Acceleration (m/s^2)','FontSize',14)
legend({'Exact','Numerical'},'location','SouthEast','FontSize',16)
hold off
%%
plot(t,residual2,'r*','MarkerSize',1)
title('Residual, approx. and exact accel.','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Numerical approx. minus exact','FontSize',14)
%%
% 6d. What happens when the initial velocity changes?  Lets use some plots
% to investigate.  The time it takes to hit the ground is obviously
% shorter, on the order of 95.627 seconds, with a slight error somewhere on
% the order of 10^-5.  The terminal velocity, also unsurprisingly, doesn't
% change!  A plot of the velocity (via numerical approximation technique)
% will be provided to show that the terminal velocity is the same as the
% previous data.

clear all;
t=0:0.001:95.627;
h=10000;
k=0.1;
g=9.81;
v_o=-160;
z=h-(g.*t)./k + ((k*v_o + g)./k^2).*(1-exp(-k.*t));
plot(t,z)
title('Displacement vs. time, v\o= -160m/s','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Displacement (m)','FontSize',14)
%%
% Here's the code that approximates the velocity for our new starting
% velocity.  It shows the velocity approaching the same value as
% previously, which makes sense with the shape of our graph and considering
% the reality that the function is trying to model. (6d)
N=length(t);
dt=0.001;
dzdt(2:N-1)=(z(3:N)-z(1:N-2))/(2*dt);
dzdt(1)=(3*dzdt(2))-(3*dzdt(3))+dzdt(4);
dzdt(N)=(3*dzdt(N-1))-(3*dzdt(N-2))+dzdt(N-3);
plot(t,dzdt,'bo-','MarkerSize',1)
title('Numerical approx. of velocity, v\o= -160m/s','FontSize',16)
xlabel('Time (s)','FontSize',14)
ylabel('Velocity (m/s)','FontSize',14)
