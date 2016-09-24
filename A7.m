%% Assignment #7
%% Problem 1.
% Integrating the complex-square of the wavefunction to evaluate the
% probability it is between L/2 and L/3 in an infinite well of length L.
% In parts b and c, we consider the case of a particle as a simple harmonic
% oscillator.
%% 1a.

%m_rule is a user-defined function that is as follows-

%function [i,h]=m_rule(f,a,b,n)
%   %Parameters
%   %f: Function to be evaluted by midpoint rule
%   %a: Beginning point of integration
%   %b: End point of integration
%   %n: Number of bins/steps
%   %Returns
%   %i: Evaluated integral (midpoint rule)
%   %h: Step size, optional.
%
%h=(b-a)/n;
%x=linspace(a,b,n+1);
%xeval=linspace(a+h/2,b-h/2,n);
%fn=feval(f,xeval);
%i=0;
%for j=1:n
%    i=i+fn(j)*(x(j+1) - x(j));
%end
%end

syms f(x) g(x);
f(x)=2*sin(pi*x)^2; g(x)=2*sin(2*pi*x)^2;
a=1/3; b=1/2;
fint_0=int(f(x),x,a,b); exact_0=eval(fint_0);
fint_1=int(g(x),x,a,b); exact_1=eval(fint_1);
fun_0=@(x) 2.*sin(pi*x).^2; fun_1=@(x) 2.*sin(2*pi*x).^2;
h_0=12; h_1=151;
est_0=m_rule(fun_0,a,b,h_0);
est_1=m_rule(fun_1,a,b,h_1);

%% Number of bins and percentage error
% The bins required is 12 for the ground state case, and 151 for the first
% excited state.  These are the minimum number of bins for an error that is
% less than 0.1%.
%% Plotting the complex-square wavefunction
xp=linspace(0,1,1111); fx_0=2.*sin(pi*xp).^2; fx_1=2.*sin(2*pi*xp).^2;
plot(xp,fx_0,'r')
hold on
plot(xp,fx_1,'b')
title('Complex-square of Wavefunction, Infinite Well, n=1,2','FontSize',14)
xlabel('x (position)','FontSize',12)
ylabel('|\psi(x)|^2','FontSize',12)
ylim([0 2]);
line([1/2 1/2],get(gca,'ylim'));
line([1/3 1/3],get(gca,'ylim'));
legend({'n=1','n=2'},'Location','NorthEastOutside')
hold off

%% Explanation of the graph for the infinite well
% Looking at the graph, notice that while n=1 has a maximum at L/2, n=2 has a
% minimum at that point, with absolutely no chance of the particle being found
% there, and likewise, there is an increasing probability from L/3 to L/2
% in the n=1 case and a rapidly decreasing probability in the n=2 case.
% The area under consideration is delineated with vertical lines.
% This matches up with our results from integration--it's about 2.7 times
% more likely that we'll find the n=1, ground state particle in that region
% as compared to the first excited state, n=2.

%% 1b.
% Now considering the particle in the simple harmonic oscillator.
fx=@(x) (exp(-x.^2/2)).^2; gx=@(x) (x.*exp(-x.^2/2)).^2;
fx_un=m_rule(fx,-1000,1000,100000); gx_un=m_rule(gx,-1000,1000,100000);
A_0=sqrt(1/sqrt(pi)); A_1=sqrt(2/sqrt(pi));
syms f(x) g(x);
f(x)=exp(-x^2/2)^2; g(x)=(x*exp(-x^2/2))^2;
sym_0=int(f(x),x,-1000,1000);
sym_1=int(g(x),x,-1000,1000);
fx_sym=eval(sym_0); gx_sym=eval(sym_1);

%Integrating numerically with my own m_rule function.
%Checking the answer with symbolic MATLAB tools.
%The same normlization constant for each!
%We can trust our numerical method now.

%% 1c.
% Calculating the probability of finding the particle in the harmonic
% oscillator between x=-0.2 and x=0.2 for the ground state and then
% the first excited state.
fx_n=@(x) (A_0.*exp(-x.^2/2)).^2;
gx_n=@(x) (A_1.*x.*exp(-x.^2/2)).^2;
a=-0.2; b=0.2; N=15;
exactP_0=eval(int(A_0^2*f(x),x,a,b));
exactP_1=eval(int(A_1^2*g(x),x,a,b));
P_0=m_rule(fx_n,a,b,N);
P_1=m_rule(gx_n,a,b,N);
xq=linspace(-1,1,1000);
fx_p=(A_0.*exp(-xq.^2/2)).^2; gx_p=(A_1.*xq.*exp(-xq.^2/2)).^2;

plot(xq,fx_p,'r')
hold on
plot(xq,gx_p,'b')

title('Complex-square of Wavefunction, S.H.O., n=0,1','FontSize',14)
xlim([-1 1]);
ylim([0 0.57]);
x1=-0.2; x2=0.2;
line([x1 x1],get(gca,'ylim'));
line([x2 x2],get(gca,'ylim'));
%Drawing lines on the graph to show the area we are interested in.

xlabel('x (position)','FontSize',12)
ylabel('|\psi(x)|^2','FontSize',12)
legend({'n=0','n=1'},'Location','NorthEastOutside')
hold off

%% Explanation of the simple harmonic oscillator graph
% Our calculations above yielded the result that while there is a 22.27%
% chance to find it in the region x=-0.2 to x=0.2 for the ground state
% case, there is only a .59% chance for the first excited state to be found
% in this region.  Looking at the graphs does demonstrate this difference.
% The area underneath the curve for the ground state is clearly many times
% larger than the same range for the excited state. It has a peak % at x=0, while the first excited state has 0 probability of being found at
% the point x=0.  The area under consideration is delineated with vertical lines. 10 bins was enough for the ground state function to
% evaluate numerically with an error < 0.1%.  It took 15 bins to increase
% the accuracy for the n=1, first excited state to the same level--error of
% less than 0.1%.
%% Problem 2.
% A charged cylinder in the x-y plane, producing an electric field at a
% certain observation point.  The distance in terms of y doesn't change,
% and the total field may be considered as a summation of point charges
% over dx.  There is some pen and paper work attached for this problem to
% show the logic of certain functions that are used.
%% 2a.
% For part a, a uniform charge density of 1x10^-6 C/m

k=8.987551787997912e9; lambda=1e-6;
theta_min=.46365; theta_max=pi/4; yp=2;
syms a(x) b(y);
a(x) = ((k*lambda)/2)*sin(x);
b(y) = ((k*lambda)/2)*cos(y);
Ex_exact=eval(int(a(x),x,theta_min,theta_max));
Ey_exact=eval(int(b(y),y,theta_min,theta_max));
% Determination of Ex and Ey vectors by symbolic integration.
ax=@(x) ((k*lambda)/yp)*sin(x);
by=@(y) ((k*lambda)/yp)*cos(y);
steps=300;
Ex_numerical=m_rule(ax,theta_min,theta_max,steps);
%At 200 steps, the error of the midpoint rule was about 10^-4.
%Increasing to 300 got rid of that 10^-4 error.
%Our numerical estimate of Ex is now accurate to at least 4 decimal points.

Ey_numerical=m_rule(by,theta_min,theta_max,steps);
%300 steps also keeps our Ey, numerically evaluated, to the same accuracy.

fprintf('Exact E-field, x-direction: %5.5f V/m\n', Ex_exact);
fprintf('Numerical approximation, x-direction: %5.5f V/m\n', Ex_numerical);
fprintf('Exact E-field, y-direction: %5.5f V/m\n', Ey_exact);
fprintf('Numerical approximation, y-direction: %5.5f V/m\n', Ey_numerical);
%% Further explanation for 2a.
% While 200 steps kept us almost within 10^-4, with only an error of 1 in
% the 4th decimal spot, 300 steps eliminated even that error for Ex and Ey,
% keeping our numerical method to a very, very accurate level.  Our
% symbolic and numerical results thus agree, giving us a total electric
% field at point P(5,3) as 841.7711 x-hat + 1168.8921 y-hat in unit
% vector notation.  The units of the electric field are Volts/meter.
%% 2b.
% Charge density is no longer constant throughout the charged rod.
% Instead, it is a function of x, lambda(x)= (2e-6*x)/b^2, where b=3.
% Because we're integrating over theta and not x, we'll have to use some
% cleverness to get the proper values of x.  The relation of x to a given
% value of theta is as follows: x=5-(2/tan(theta));  With this relation we
% can now express the charge at a point on the rod in terms of theta, not
% just x.

k=8.987551787997912e9; lambda_0=2e-6;
theta_min=.46365; theta_max=pi/4; yp=2;
b=3;
syms c(m) d(n);
c(m) = ((k*(lambda_0/b^2)*(5-(2/tan(m))))/2)*sin(m);
d(n) = ((k*(lambda_0/b^2)*(5-(2/tan(n))))/2)*cos(n);
Ex_exact2=eval(int(c(m),m,theta_min,theta_max));
Ey_exact2=eval(int(d(n),n,theta_min,theta_max));
% Determination of Ex and Ey vectors by symbolic integration.
steps=800;
cm=@(m) ((k*(lambda_0/b^2).*(5-(2./tan(m))))./yp).*sin(m);
dn=@(n) ((k*(lambda_0/b^2).*(5-(2./tan(n))))./yp).*cos(n);
Ex_numerical2=m_rule(cm,theta_min,theta_max,steps);
Ey_numerical2=m_rule(dn,theta_min,theta_max,steps);
%Determination of Ex and Ey vectors by numerical analysis.

fprintf('Exact E-field, x-direction: %5.5f V/m\n', Ex_exact2);
fprintf('Numerical approximation, x-direction: %5.5f V/m\n', Ex_numerical2);
fprintf('Exact E-field, y-direction: %5.5f V/m\n', Ey_exact2);
fprintf('Numerical approximation, y-direction: %5.5f V/m\n', Ey_numerical2);
%% Further explanation for 2b.
% Despite having to change the terms due to integrating with respect to theta
% rather than x, the same methods worked just fine as long as we reframed the
% charge distribution as a function of theta rather x as well. In terms of unit
% vector notation, the field at P can be expressed as 416.2381 x-hat + 548.8207 y-hat.
% This is a weaker field than from part a, a quick consideration of the
% strength of the field makes this an assuring result that makes some
% sense.  Interesting to note is the number of steps we had to take to get a similar
% error, less in magnitude for each vector than 10^-4... 300 steps was no
% longer enough, but rather, we had to increase the number to 800 steps to
% achieve a similar level of accuracy!  This makes sense considering that
% non-uniform charge density complicates things and requires a higher
% amount of samples to achieve the same level of absolute error.
%% Problem 3.
% Plotting a particle with a known velocity function, and using numerical
% and symbolic MATLAB tools to calculate its position and acceleration.
%% 3a.
dt=0.005;
v_o=2.3; v_a=10; omega=2; gamma=0.2; tq=0:dt:10;
v_func=v_o+v_a.*cos(omega.*tq).*exp(-gamma.*tq);
plot(tq,v_func,'b*','MarkerSize',1)
title('Velocity of a particle with time','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Velocity (m/s)','FontSize',12)
%% 3b.
N=length(tq);
dvdt(2:N-1)=(v_func(3:N)-v_func(1:N-2))/(2*dt);
%Center finite difference method of numerical approximation.

dvdt(1)=2*dvdt(2)-dvdt(3);
dvdt(N)=2*dvdt(N-1)-dvdt(N-2);
%Manual linear extrapolation at the boundaries of our interval.

syms v(t)
v(t)=v_o+v_a*cos(omega*t)*exp(-gamma*t);
v_M=@(t) v_o+v_a.*cos(omega*t).*exp(-gamma*t);
sym_accel=diff(v(t),t);
accel_func= -2.*cos(2.*tq).*exp(-tq./5) -20.*sin(2.*tq).*exp(-tq./5);
plot(tq,accel_func,'b*','MarkerSize',12)
hold on
plot(tq,dvdt,'r*','MarkerSize',4)
title('Numerical and Exact calculation of acceleration','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Acceleration (m/s^2)','FontSize',12)
legend({'Exact function','Numerical approximation'},'Location','SouthEast')
hold off

%%
residual=abs(accel_func-dvdt);
plot(tq,residual,'b*')
title('Absolute residual of acceleration','FontSize',14)
ylabel('Difference of numerical and exact (m/s^2)','FontSize',12)
xlabel('Time (s)','FontSize',12)

%% Commentary on the residual and step-size
% This plot of the residual shows that 0.05 is a small enough step-size.
% Our goal was to achieve an error smaller than 10^-3 for the entire
% interval.  Even with bigger error on the very ends due to the linear
% extrapolation, this step size achieves the desired level of accuracy.
%% 3c.
a=0; N=length(tq); bins=200;
pos_exact=(23.*tq)/10 - (50.*cos(2.*tq).*exp(-tq./5))./101 + (500.*sin(2.*tq).*exp(-tq./5))./101 + 1.495;
pos_approx=ones(1,N);
for k=1:N
 pos_approx(k)=m_rule(v_M,a,a+((k-1)*dt),bins)+1;
end
plot(tq,pos_exact,'b*','MarkerSize',12)
hold on
plot(tq,pos_approx,'r*','MarkerSize',4)
title('Numerical and Exact calculation of position','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Displacement (m)','FontSize',12)
legend({'Exact function','Numerical approximation'},'Location','SouthEast')
hold off
%Solving the initial value problem; x(0)=1m

%%
residual=abs(pos_approx-pos_exact);
plot(tq,residual,'b*')
title('Absolute residual of position','FontSize',14)
ylabel('Difference of numerical and exact (m)','FontSize',12)
xlabel('Time (s)','FontSize',12)
%Plotting the residual for 3c.

%% Examining the results
% Our approximation gives us a really good answer with 200 bins.  The error
% is always on the order of 10^-4 in terms of absolute magnitude, that's
% less than a 0.01% error across the entire interval!
%% Problem 4.
% A disc of uniform charge density denoted by sigma.  There is pen & paper
% work attached to the back of this assignment that shows what the proper
% integral for evaluating the electric field at a point located above the
% center of the disc, a distance z away, can be evaluated as, as an
% integral summing over small changes in radius, dr.
%% 4a.
% Comparison of MATLAB's symbolic calculation with the closed form solution
% of the integral obtained by pencil & paper method.

sigma=1e-6; R=0.1; z=1; R_min=0;
e_0=8.85418782e-12;
syms E(r);
E(r)=(sigma*z*r)/(2*e_0*(z^2+r^2)^(3/2));
sym_E=eval(int(E(r),r,0,R));
fprintf('MATLAB symbolic evaluation, z-direction: %5.5f V/m\n',sym_E)

%MATLAB's calculation agrees with our pen&paper work.

%% 4b.
% Now, we turn our attempt to the Monte Carlo technique.  A sampling of one
% million random numbers in the range of our integral seems to give us a
% rather consistent answer, with error less than 0.01% reliably.  Note that
% the method is *random* however, and some runs will produce an error
% greater than this.  Our exact solution, in unit vector notation, is
% 280.2521 z-hat.

MCN=1e6;
er1=@(r) (sigma*z.*r)./(2*e_0.*(z^2+r.^2).^(3/2));
rq1=R_min+(R-R_min)*rand(1,MCN);
field_MC=(R-R_min)/MCN*sum(er1(rq1));
fprintf('Electric field, z-direction, for first case: %5.5f V/m',field_MC)
%% 4c.
% As with part b, we will be using the Monte Carlo technique, but the
% charge density is no longer uniform, but rather, dependent upon the point
% on our disc as defined by radius r, squared.  Luckily, dealing with this
% change will be fairly simple--we'll just need to introduce an additional
% r^2 term in to our function, and run the same sort of Monte Carlo
% technique as before.  We'll keep the number of random sampling on our
% interval at one million--the method has a great many advantages, but a
% quick convergence is not one of them!  We get a much weaker field than in
% part b--this makes sense though, because r is a small number, and r^2 is
% of course then, even smaller!  In unit vector notation, the field at
% point P, a distance z above the center of the disc, is 1.3973 z-hat.

MCN2=1e6;
er2=@(r) (sigma*z.*r.^3)./(2*e_0.*(z^2+r.^2).^(3/2));
rq2=R_min+(R-R_min)*rand(1,MCN2);
field_MC2=(R-R_min)/MCN2*sum(er2(rq2));
fprintf('Electric field, z-direction, for second case: %5.5f V/m',field_MC2)
%% Problem 5.
a=2; b=0; c=pi; d=0; e=2*pi; f=0; N=1000000;
%Our limits for integration, and number of random #s for Monte Carlo
%numerical integration
rq=b+(a-b)*rand(1,N);
tq=d+(c-d)*rand(1,N);
dV_MCI=((a-b)*(c-d)*(e-f))/N*sum(rq.^2.*sin(tq));
%Monte Carlo approximation of known analytical solution.
dV_exact=(4/3)*pi*a^3;
error_pct=abs(100-(dV_MCI/dV_exact * 100));
fprintf('Monte Carlo result: %5.5f m^3\n',dV_MCI);
fprintf('Analytical result: %5.5f m^3\n',dV_exact);
fprintf('Error: %5.5f percent\n\n',error_pct);
%To get the answer to within 0.1% of the analytical result consistently,
%around 1 million random numbers are needed by the Monte Carlo method.
%% 5b.
a=100; b=0; c=pi; d=0; e=2*pi; f=0; N=1000000;
rq=b+(a-b)*rand(1,N);
tq=d+(c-d)*rand(1,N);
I_MC=((a-b)*(c-d)*(e-f))/N*sum(exp(-rq).^2.*rq.*sin(tq));

syms f(r) g(t);
f(r)=r*exp(-r).^2;
g(t)=sin(t);
r_0=int(f(r),r,b,a);
t_0=int(g(t),t,d,c);
r_exact=eval(r_0); t_exact=eval(t_0); p_exact=2*pi;
I_exact=r_exact*t_exact*p_exact;
fprintf('Monte Carlo result: %5.5f\n',I_MC);
fprintf('Symbolic result: %5.5f\n\n', I_exact);

%%

a=100; b=0; c=pi; d=0; e=2*pi; f=0; N=1000000;
rq=b+(a-b)*rand(1,N);
tq=d+(c-d)*rand(1,N);
I_MC2=((a-b)*(c-d)*(e-f))/N*sum(rq.^2.*sin(tq).*(exp(-rq./2).*cos(tq).*rq).^2);
syms f(r) g(t);
f(r)= r^4*exp(-r/2)^2;
g(t)=sin(t)*cos(t)^2;
r_02=int(f(r),r,b,a);
t_02=int(g(t),t,d,c);
p_exact2=2*pi;
r_exact2=eval(r_02); t_exact2=eval(t_02);
I_exact2=r_exact2*t_exact2*p_exact2;
fprintf('Monte Carlo result: %5.5f\n',I_MC2);
fprintf('Symbolic result: %5.5f\n\n', I_exact2);

A_100=1/sqrt(pi); A_210=1/sqrt(32*pi);
fprintf('Normalization constant for ground state: %5.5f\n',A_100)
fprintf('Normalization constant for excited state: %5.5f\n',A_210)
%
%% Normalization constants for hydrogen atom
% We've determined both normalization constants now!  For the ground state,
% it's 1/sqrt(pi).  For the 2,1,0 state, it's 1/sqrt(32*pi).
%% 5c.
for i=1:6
a=i*2; b=(i-1)*2; c=pi; d=0; e=2*pi; f=0; N=1000000;
rq=b+(a-b)*rand(1,N);
tq=d+(c-d)*rand(1,N);
I_MC=((a-b)*(c-d)*(e-f))/N*sum((1/pi).*exp(-rq).^2.*rq.*sin(tq));
syms f(r) g(t);
f(r)=(1/pi)*r*exp(-r).^2;
g(t)=sin(t);
r_0=int(f(r),r,b,a);
t_0=int(g(t),t,d,c);
r_exact=eval(r_0); t_exact=eval(t_0); p_exact=2*pi;
I_exact=r_exact*t_exact*p_exact;
PROB_100(i)=I_MC*100;
PROB_100_exact(i)=I_exact*100;
end
%%
for i=1:6
a=i*2; b=(i-1)*2; c=pi; d=0; e=2*pi; f=0; N=1000000;
rq=b+(a-b)*rand(1,N);
tq=d+(c-d)*rand(1,N);
I_MC2=((a-b)*(c-d)*(e-f))/N*sum((1/(32*pi)).*rq.^2.*sin(tq).*(exp(-rq./2).*cos(tq).*rq).^2);
syms f(r) g(t);
f(r)=(1/(32*pi))*r^4*exp(-r/2)^2;
g(t)=sin(t)*cos(t)^2;
p_exact2=2*pi;
r_02=int(f(r),r,b,a);
t_02=int(g(t),t,d,c);
r_exact2=eval(r_02); t_exact2=eval(t_02);
I_exact2=r_exact2*t_exact2*p_exact2;
PROB_210(i)=I_MC2*100;
PROB_210_exact(i)=I_exact2*100;
end
Reg={'0 to 2', '2 to 4', '4 to 6', '6 to 8', '8 to 10', '10 to 12'};
Region=transpose(Reg);
nlm100=transpose(PROB_100);
nlm100exact=transpose(PROB_100_exact);
nlm210=transpose(PROB_210);
nlm210exact=transpose(PROB_210_exact);
T=table(Region,nlm100,nlm100exact,nlm210,nlm210exact);
T.Properties.VariableUnits = {'' '%' '%' '%' '%'};
disp(T)
%% Commentary on the table of probabilities
% Each row consists of a region of space for the possible radius, 0 to 2, 2
% to 4, 4 to 6, 6 to 8, 8 to 10, and 10 to 12.  The units for the probability is percentage.  The probabilities were
% converted to percentages by multiplying the probability mass function by
% 100.  For instance, an electron in the ground state will be found in the
% region of r between 0 and 2 around 90.9% of the time; while for the
% excited 2,1,0 state it only has a ~5.2% chance of being found in the same
% radial region.
%%
c=pi; d=0; e=2*pi; f=0; N=1000000;
for i=1:6
    
r_min=(i-1)*2;
r_max=i*2;
dr=0.001;
r_i=r_min:dr:r_max;
psi_100=r_i.^2.*A_100^2.*exp(-r_i).^2;
psi_210=(1/2)*r_i.^2.*(A_210.*r_i.*exp(-r_i/2)).^2;
plot(r_i,psi_100,'r*','MarkerSize',1)
title('Radial probability functions, n,l,m=1,0,0 & 2,1,0','FontSize',14)
xlabel('r (no units)','FontSize',12)
ylabel('r^2|R(r)|^2','FontSize',12)
hold on
plot(r_i,psi_210,'b*','MarkerSize',1)
end
legend({'1,0,0','2,1,0'},'Location','NorthEast');
line([2 2],get(gca,'ylim'));
line([4 4],get(gca,'ylim'));
line([6 6],get(gca,'ylim'));
line([8 8],get(gca,'ylim'));
line([10 10],get(gca,'ylim'));
%% Commentary on the graph of radial probability
% One thing to notice is where the functions have maximums--at 1 'r' for the
% ground state, and 4 'r' for the 2,1,0 state--this is good news for us,
% because it matches up with the idea of the Bohr radius (We didn't include
% units here, but it still matches up with the proper number of Bohr radii
% for these electron states.  The amount of area underneath each curve for
% each section of 0 to 2, 2 to 4, etc. r values seems to be a match for the
% table of probabilities.