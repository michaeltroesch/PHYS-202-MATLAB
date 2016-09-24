%% Assignment #8
%% Problem 1.
% The exact analytical solution, and a comparison to its solution via the
% finite difference method of numerical approximation--for the second-order
% linear differential equation, y'' + y' - 2y = 2t.
%% 1a. and 1b.

clear all;
dt=0.004; t=0:dt:4;
y_analytical=exp(t)-(exp(-2*t)/2)-(1/2)-t;
%Closed form solution worked out in paperwork attached to the back.

plot(t,y_analytical,'b*','MarkerSize',6)
hold on
title('d^2y/dt^2 + dy/dt - 2y = 2t, y(0)=0, dy/dt(0)=1','FontSize',16)
xlabel('t (seconds)','FontSize',14)
ylabel('y (meters)','FontSize',14)

N=length(t);
y(1)=0;
v(1)=1;
for j=1:N-1
    y(j+1)=y(j)+(v(j)*dt);
    v(j+1)=v(j)+((2*t(j))+(2*y(j))-v(j))*dt;
end

plot(t,y,'r*','MarkerSize',1)
legend({'Analytical','Finite Difference Method'},'FontSize',12,'Location','NorthEast')
hold off

%% Step size and the numerical method used
% My initial step size was dt=0.01, which did indeed satisfy not having an
% absolute error no greater than 0.5m over the entire interval.  A quick
% check of other options showed that while dt=0.02 also satisfied our
% requirement for accuracy, dt=0.03 was just barely too large near the end
% of the interval.  However, this was for Runge-Kutta's 2nd order method.
% The Euler's method that we switched to requires a smaller step-size.  To
% reduce the absolute error at any point to less than 0.5m for Euler's
% method, required a step-size dt=0.004, much smaller!
%% Problem 2.
% Modeling a projectile with resistive force and at different initial
% velocities and starting angles.  Determining the starting angle that
% produces the greatest range for the projectile.
%% 2a.

clear all;
dt=0.001;
t=0:dt:17;
N=length(t);

b=0;
m=1;

g=9.8;
h=0;

theta_0=pi/4;
v(1)=45;

vx(1)=v(1)*cos(theta_0);
vy(1)=v(1)*sin(theta_0);

sx_analytical=vx(1).*t;
sy_analytical=h + vy(1).*t - (1/2)*g*t.^2;

ax(1)=(-b*sqrt(vx(1)^2+vy(1)^2)*vx(1))*(1/m);
ay(1)=(-b*sqrt(vx(1)^2+vy(1)^2)*vy(1))*(1/m) - g;

sx(1)=0;
sy(1)=h;

R_marker=0;

for j=1:N-1
    vx(j+1)=vx(j) + ax(j)*dt;
    vy(j+1)=vy(j) + ay(j)*dt;
    ax(j+1)=(-b*sqrt(vx(j+1)^2+vy(j+1)^2)*vx(j+1))*(1/m);
    ay(j+1)=(-b*sqrt(vx(j+1)^2+vy(j+1)^2)*vy(j+1))*(1/m) - g;
    sx(j+1)=sx(j) + vx(j)*dt;
    sy(j+1)=sy(j) + vy(j)*dt;
    
    if(sy(j+1) <= 0 && R_marker == 0)
        R_marker=j+1;
    end
end

R=sx(R_marker);
R_analytical=sx_analytical(R_marker);
difference=abs(sy_analytical-sy);


plot(sx_analytical,sy_analytical,'b*','MarkerSize',6)
hold on
plot(sx,sy,'r*','MarkerSize',1)
text(50,12,'Range, analytical = 206.6696m')
text(50,10,'Range, numerical = 206.6696m')
ym=max(sy)+4.2;
ylim([0 ym])
title('Trajectory of a projectile in x-y','FontSize',14)
xlabel('x (meters)','FontSize',14)
ylabel('y (meters)','FontSize',14)
hold off

%% The step-size dt required for a certain accuracy
% Our goal of accuracy against the analytical solution was a maximum value
% less than 0.15m across the entire interval.  A step-size of 0.002 was not
% sufficient to achieve this accuracy with the chosen numerical method,
% Euler's finite difference method.  At dt=0.001, we get the desired
% accuracy across the entire interval.
%% The analytical solution of theta max
% For the case with no air resistance, we can look just at the theta which
% maximizes velocity in x and y simultaneously.  Aside from the constant
% starting velocity which we decide, x and y are functions of cosine and
% sine.  The first derivative of cos + sin (multiplied by a constant, in
% this case the initial velocity) is cos - sin.  The second derivative
% then, is -cos - sin.  Given that we only need to consider the range of 0
% to pi/2, we quickly find that there is a local maxima at pi/4 which is
% our theta max.
%% 2b.

clear all;
dt=0.001;
t=0:dt:17;
N=length(t);

b=0;
m=1;

g=9.8;
h=2;

percents=[.98 1 1.2];
for k=1:length(percents)
theta_0=percents(k)*(pi/4);
v(1)=45;

vx(1)=v(1)*cos(theta_0);
vy(1)=v(1)*sin(theta_0);

sx_analytical=vx(1).*t;
sy_analytical=h + vy(1).*t - (1/2)*g*t.^2;

ax(1)=(-b*sqrt(vx(1)^2+vy(1)^2)*vx(1))*(1/m);
ay(1)=(-b*sqrt(vx(1)^2+vy(1)^2)*vy(1))*(1/m) - g;

sx(1)=0;
sy(1)=h;

R_marker=0;

for j=1:N-1
    vx(j+1)=vx(j) + ax(j)*dt;
    vy(j+1)=vy(j) + ay(j)*dt;
    ax(j+1)=(-b*sqrt(vx(j+1)^2+vy(j+1)^2)*vx(j+1))*(1/m);
    ay(j+1)=(-b*sqrt(vx(j+1)^2+vy(j+1)^2)*vy(j+1))*(1/m) - g;
    sx(j+1)=sx(j) + vx(j)*dt;
    sy(j+1)=sy(j) + vy(j)*dt;
    
    if(sy(j+1) <= 0 && R_marker == 0)
        R_marker=j+1;
    end
end

R(k)=sx(R_marker);
R_analytical(k)=sx_analytical(R_marker);
difference=abs(sy_analytical-sy);


plot(sx_analytical,sy_analytical,'*','MarkerSize',3)
hold on
plot(sx,sy,'*','MarkerSize',1)
end

str1=strcat('Range 98%, numerical =', num2str(R(1)),'m');
str1a=strcat('Range 98%, numerical =', num2str(R_analytical(1)),'m');
str2=strcat('Range 100%, numerical =', num2str(R(2)),'m');
str2a=strcat('Range 100%, analytical =', num2str(R_analytical(2)),'m');
str3=strcat('Range 102%, numerical =', num2str(R(3)),'m');
str3a=strcat('Range 102%, analytical =', num2str(R_analytical(3)),'m');

text(50,24,str1)
text(50,21,str1a)
text(50,18,str2)
text(50,15,str2a)
text(50,12,str3)
text(50,9,str3a)
title('Trajectory of a projectile in x-y','FontSize',14)
xlabel('x (meters)','FontSize',14)
ylabel('y (meters)','FontSize',14)
ym=max(sy)+4.2;
ylim([0 ym])
hold off

pdiff1=100-100*(R(1)/R(2));
pdiff2=100-100*(R(3)/R(2));

%% The calculated percentage differences
% Subtracting 2% from theta max gives us a mere .0212 percentage difference
% in the range of our projectile for these initial conditions.  Adding two
% percent to theta max results in a 5.1118 percentage difference in the
% range of the projectile.  The ranges are denoted on the plot.
%% 2c.

clear all;
dt=0.001;
t=0:dt:17;
N=length(t);
db=0.02;

m=1;

g=9.8;
h=2;

theta_0=pi/4;
v(1)=40;

vx(1)=v(1)*cos(theta_0);
vy(1)=v(1)*sin(theta_0);

sx_analytical=vx(1).*t;
sy_analytical=h + vy(1).*t - (1/2)*g*t.^2;

num_trajectories=5;

for k=0:num_trajectories
    b(k+1)=k*db;
    ax(1)=(-b(k+1)*sqrt(vx(1)^2+vy(1)^2)*vx(1))*(1/m);
    ay(1)=(-b(k+1)*sqrt(vx(1)^2+vy(1)^2)*vy(1))*(1/m) - g;

    sx(1)=0;
    sy(1)=h;

    R_marker=0;

    for j=1:N-1
        vx(j+1)=vx(j) + ax(j)*dt;
        vy(j+1)=vy(j) + ay(j)*dt;
        ax(j+1)=(-b(k+1)*sqrt(vx(j+1)^2+vy(j+1)^2)*vx(j+1))*(1/m);
        ay(j+1)=(-b(k+1)*sqrt(vx(j+1)^2+vy(j+1)^2)*vy(j+1))*(1/m) - g;
        sx(j+1)=sx(j) + vx(j)*dt;
        sy(j+1)=sy(j) + vy(j)*dt;
    
        if(sy(j+1) <= 0 && R_marker == 0)
            R_marker=j+1;
        end
    end

    R(k+1)=sx(R_marker);
    R_analytical=sx_analytical(R_marker);
    difference=abs(sy_analytical-sy);

    plot(sx,sy)
    hold on
end

str1=strcat('b=0.00',' N m^2/s^2 ',' R=',num2str(R(1)),'m');
str2=strcat('b=',num2str(b(2)),' N m^2/s^2 ',' R=',num2str(R(2)),'m');
str3=strcat('b=',num2str(b(3)),' N m^2/s^2 ',' R=',num2str(R(3)),'m');
str4=strcat('b=',num2str(b(4)),' N m^2/s^2 ',' R=',num2str(R(4)),'m');
str5=strcat('b=',num2str(b(5)),' N m^2/s^2 ',' R=',num2str(R(5)),'m');
text(55,25,str1)
text(55,23,str2)
text(55,21,str3)
text(55,19,str4)
text(55,17,str5)
ylim([0 50])
title('Fifty projectiles in x-y with varied air resistance','FontSize',14)
xlabel('x (meters)','FontSize',14)
ylabel('y (meters)','FontSize',14)
hold off

%% Determination of the max range for b = 0.05 N m^2 s^-2
% Altering the script in part 2b slightly, checking different percentages
% of pi/4 to find a maximum for this value of b, a maximum was found at 68
% percent of pi/4, which would give us a theta max for this drag
% coefficient of 0.17 * pi radians.
%% Problem 3.
% A non-linear simple harmonic oscillator.  New from the consideration in a
% previous worksheet is the addition of a damping and driving force.  Three
% new constants are introduced that can be tweaked to influence the
% physical scenario--gamma, the damping force factor.  A--which influences
% the amplitude of the driving force, and B (capital Omega) which influences the
% period of the driving force.  We'll also plot kinetic, potential, and
% total energy over the time interval (all these energy values will be on a
% single plot).
%% 3a.
% A similar problem to the one we solved in class--the damping and driving
% terms are included, but the constants that govern them are all set to
% zero.  Graphs for initial angles of pi/4, pi/2, and .95*pi.

clear all;
t=0:0.01:200;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0; B=0; gamma=0;
starting_angles=[pi/4, pi/2, .95*pi];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    KE(1)=(1/2)*m*L^2*omega(1)^2;
    PE(1)=m*g*((1-cos(theta(1)))*L);
    TE(1)=KE(1)+PE(1);
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
        KE(j+1)=(1/2)*m*L^2*omega(j+1)^2;
        PE(j+1)=m*g*((1-cos(theta(j+1)))*L);
        TE(j+1)=KE(j+1)+PE(j+1);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on

    figure(3)
    plot(t,omega)
    hold on
    
    figure(k+3)
    plot(t,PE,'b')
    hold on
    plot(t,KE,'k')
    plot(t,TE,'r')
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Nonlinear SHO','FontSize',14)
legend({'\theta_o=\pi/4 rad','\theta_o=\pi/2 rad','\theta_o=.95\pi rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Nonlinear SHO','FontSize',14)
legend({'\theta_o=\pi/4 rad','\theta_o=\pi/2 rad','\theta_o=.95\pi rad'},'Location','NorthEast')
hold off

figure(3)
xlabel('time (s)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('\omega vs. time for Nonlinear SHO','FontSize',14)
legend({'\theta_o=\pi/4 rad','\theta_o=\pi/2 rad','\theta_o=.95\pi rad'},'Location','NorthEast')
hold off

figure(4)
xlabel('time (s)','FontSize',14)
ylabel('Energy (Joules)','FontSize',14)
title('Energy vs. time for Nonlinear SHO','FontSize',14)
legend({'Potential','Kinetic','Total'},'Location','NorthEast')
hold off

figure(5)
xlabel('time (s)','FontSize',14)
ylabel('Energy (Joules)','FontSize',14)
title('Energy vs. time for Nonlinear SHO','FontSize',14)
legend({'Potential','Kinetic','Total'},'Location','NorthEast')
hold off

figure(6)
xlabel('time (s)','FontSize',14)
ylabel('Energy (Joules)','FontSize',14)
title('Energy vs. time for Nonlinear SHO','FontSize',14)
legend({'Potential','Kinetic','Total'},'Location','NorthEast')
hold off

%% The shape of the phase space in part a.
% For the higher starting angle, .95*pi, we get a sort of teardrop,
% eyeball, or almond shape--so it's worthwhile to consider what we are
% looking at when we see the phase space graph.  The x-axis is theta, while
% the y-axis is omega.  For the lower angles, pi/4 and pi/2, the graph is
% very much a circle.  The critical point where the graph begins to shape
% is after pi/2--rather than a smooth transition between theta and omega,
% we're looking at a ball that takes much longer to reach zero angular
% velocity near the end of its arc.  That's the reason for the shape.  Also
% worth considering is that it being on the outside of the other phase
% spaces makes sense--the maximums for omega and theta both become greater
% for greater initial values of theta. (Although theta=pi will be an
% unstable equilibrium corresponding to the ball being placed directly
% above the pivot vertically.)
%% Energy conservation in part a.
% Yes, energy is conserved.  There are many ways to reach the logical
% conclusion, outside of looking at graphs of Energy over time for the
% oscillator in this case.  My personal argument it this--the pendulum does
% absolutely no work on its environment.  No outside forces are applied.
% It is not surprising then, that total energy is constant over the
% interval.
%% 3b.

t=0:0.01:200;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0; B=0; gamma=0.06;
starting_angles=[.95*pi];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    KE(1)=(1/2)*m*L^2*omega(1)^2;
    PE(1)=m*g*((1-cos(theta(1)))*L);
    TE(1)=KE(1)+PE(1);
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
        KE(j+1)=(1/2)*m*L^2*omega(j+1)^2;
        PE(j+1)=m*g*((1-cos(theta(j+1)))*L);
        TE(j+1)=KE(j+1)+PE(j+1);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on

    figure(3)
    plot(t,omega)
    hold on
    
    figure(k+3)
    plot(t,PE,'b')
    hold on
    plot(t,KE,'k')
    plot(t,TE,'r')
end

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=5;
    KE(1)=(1/2)*m*L^2*omega(1)^2;
    PE(1)=m*g*((1-cos(theta(1)))*L);
    TE(1)=KE(1)+PE(1);
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
        KE(j+1)=(1/2)*m*L^2*omega(j+1)^2;
        PE(j+1)=m*g*((1-cos(theta(j+1)))*L);
        TE(j+1)=KE(j+1)+PE(j+1);
    end

    figure(5)
    plot(theta,omega)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Damped Nonlinear SHO \gamma=0.06','FontSize',14)
legend({'\theta_o=.95\pi'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Damped Nonlinear SHO \gamma=0.06','FontSize',14)
legend({'\theta_o=.95\pi rad'},'Location','NorthEast')
hold off

figure(3)
xlabel('time (s)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('\omega vs. time for Damped Nonlinear SHO \gamma=0.06','FontSize',14)
legend({'\theta_o=.95\pi rad'},'Location','NorthEast')
hold off

figure(4)
xlabel('time (s)','FontSize',14)
ylabel('Energy (Joules)','FontSize',14)
title('Energy vs. time for Damped Nonlinear SHO \gamma=0.06','FontSize',14)
legend({'Potential','Kinetic','Total'},'Location','NorthEast')
hold off

%% Energy conservation in part b.
% A quick look at the graphs could tell you, in this situation, energy is
% not conserved--but why?  Consider the damping force that is at work here.
% As the pendulum moves, it does work on the environment as whatever
% medium it is travelling through pushes back against it.  The energy
% declines over time, approaching 0 as t goes to infinity.
%% 3b. (Initial angular velocity added)

figure(5)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Damped Nonlinear SHO \gamma=0.06','FontSize',14)
legend({'\theta_o=.95\pi rad, \omega_o=5 rad/s'},'Location','NorthEast')
hold off

%% Phase space w/ initial angular velocity of 5 radians/second
% So we have our interesting case of the damped nonlinear oscillator now
% having an initial angular velocity.  The motion that occurs in this case
% can be described by looking at the phase space graph.  The mass makes
% many full rotations around the pivot, but as time goes on, the damping
% force reduces the total energy, so it eventually falls in to the same
% pattern as the oscillator did at an initial angular velocity of zero.
%% 3c. (Driving turned off.)

t=0:0.01:200;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0; B=1.1; gamma=0.00;
starting_angles=[0.1];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    KE(1)=(1/2)*m*L^2*omega(1)^2;
    PE(1)=m*g*((1-cos(theta(1)))*L);
    TE(1)=KE(1)+PE(1);
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

%% 3c. (Driving turned on.)

t=0:0.01:200;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0.101; B=1.1; gamma=0.00;
starting_angles=[0.1];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Driven Nonlinear SHO A=0.101 \Omega=1.1','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Driven Nonlinear SHO A=0.101 \Omega=1.1','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

%% Explanation of the differences with the driving force, part c.
% The phase space graph in this case can be difficult to interpret.  But
% generally, it shows that there is a range of amplitudes for our motion
% over the interval.  The graph of angular position shows us how that looks
% in a bit more of a clear manner--the movement starts out small, grows
% larger, and then goes down again.  This is because of the nature of the
% driving force--sometimes, it is in phase with the motion of our pendulum,
% pushing it to greater heights.  However, as it falls out of phase, it
% begins lowering the range of its motion.
%% 3d. (Driving turned off.)

t=0:0.01:400;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0; B=1.1; gamma=0.00;
starting_angles=[0.01];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    KE(1)=(1/2)*m*L^2*omega(1)^2;
    PE(1)=m*g*((1-cos(theta(1)))*L);
    TE(1)=KE(1)+PE(1);
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, \theta vs. \omega, Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.1 rad'},'Location','NorthEast')
hold off

%% 3d. (Driving turned on.)
t=0:0.01:400;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0.02; B=1.1; gamma=0.00;
starting_angles=[0.01];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Driven Nonlinear SHO A=0.02, \Omega=1.1','FontSize',14)
legend({'\theta_o=0.01 rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Driven Nonlinear SHO A=0.02, \Omega=1.1','FontSize',14)
legend({'\theta_o=0.01 rad'},'Location','NorthEast')
hold off

%% Explanation of the differences with the driving force, part d.
% When a force is resonant to the oscillator, it is corresponding to the
% peaks and troughs appropriately, growing the amplitude of its motion.
% That is what occurs here.  With a starting angle of 0.01 radians, that
% remains the maximum for the non-driven case.  The introduction of the
% resonant force however, drives this all the way up to 0.2 radians before
% falling back down again--a twenty-fold increase.
%% 3e.

t=0:0.01:400;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0.1; B=2/3; gamma=0.06;
starting_angles=[0.95*pi];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Driven-Damped Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.95\pi rad'},'Location','NorthEast')
hold off

%% Identifying the transient and steady-state parts of the graph
% The transient part is the graph up until around 200~250 seconds--after
% that point, it falls in to a steady motion.  Before that point, it's
% maximum and minimums are varying, trending downward towards the
% steady-state function.  If the initial angle is decreased, it will reach
% the steady-state part at a faster rate as compared to starting with theta
% equal to 0.95*pi radians.
%% 3f.

t=0:0.01:400;
N=length(t);
dt=(t(end)-t(1))/N;
g=9.8; L=9.8; m=9.8;

A=0.7; B=0.66; gamma=0.2;
starting_angles=[0.6 0.6001];

for k=1:length(starting_angles)
    theta(1)=starting_angles(k);
    omega(1)=0;
    
    for j=1:N-1
        theta_m=theta(j)+(omega(j)*(dt/2));
        omega_m=omega(j) - (((g/L * sin(theta(j))) + (gamma*omega(j)) - (A*cos(B*t(j))))*(dt/2));
        theta(j+1)=theta(j) + omega_m*dt;
        omega(j+1)=omega(j) - (((g/L * sin(theta_m)) + (gamma*omega_m) - (A*cos(B*t(j))))*dt);
    end

    figure(1)
    plot(theta,omega)
    hold on

    figure(2)
    plot(t,theta)
    hold on
end

figure(1)
xlabel('\theta (radians)','FontSize',14)
ylabel('\omega (radians/sec)','FontSize',14)
title('Phase space, Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.6 rad','\theta_o=0.6001 rad'},'Location','NorthEast')
hold off

figure(2)
xlabel('time (s)','FontSize',14)
ylabel('\theta (radians)','FontSize',14)
title('\theta vs. time for Chaotic Nonlinear SHO','FontSize',14)
legend({'\theta_o=0.6 rad','\theta_o=0.6 rad'},'Location','NorthEast')
hold off

%% Explaining the results for our Chaotic, Driven-Damped, Nonlinear Oscillator
% The tiniest of changes in our initial angle, 0.6 -> 0.6001 produces a
% comparatively huge difference in the motion our simulation predicts! The
% thing to notice is that our motion in this case is *chaotic*--*not
% random*, our simulation, after all, is capable of predicing it...
% however, even the slightest of changes will throw our predictions in
% wildly different areas.
%% Problem 4.
% The modelling of atmospheric convection by three coupled ordinary linear
% differential equations--solved by a finite difference method.
%% 4a.

clear all;
sigma=10; b=8/3; r=28;
dt=0.001;
t=0:dt:40;
x(1)=1; y(1)=1; z(1)=20;
N=length(t);

for j=1:N-1
    x_m=x(j) + (sigma*(y(j)-x(j)))*(dt/2);
    y_m=y(j) + (r*x(j) - y(j) - x(j)*z(j))*(dt/2);
    z_m=z(j) + ((x(j)*y(j) - b*z(j)))*(dt/2);
    x(j+1)=x(j) + (sigma*(y_m - x_m))*dt;
    y(j+1)=y(j) + (r*x_m - y_m - x_m*z_m)*dt;
    z(j+1)=z(j) + (x_m*y_m - b*z_m)*dt;
end

figure(1)
plot(t,x,'b')
title('x(t), finite difference method','FontSize',14)
xlabel('time (seconds)','FontSize',14)
ylabel('x(t) (unspecified units)','FontSize',14)
figure(2)
plot(t,y,'b')
title('y(t), finite difference method','FontSize',14)
xlabel('time (seconds)','FontSize',14)
ylabel('y(t) (unspecified units)','FontSize',14)
figure(3)
plot(t,z,'b')
title('z(t), finite difference method','FontSize',14)
xlabel('time (seconds)','FontSize',14)
ylabel('z(t) (unspecified units)','FontSize',14)

%% 4b.

figure(4)
plot(y,x)
title('Phase Space, x vs. y','FontSize',14)
xlabel('y(t) (unspecified units)','FontSize',14)
ylabel('x(t) (unspecified units)','FontSize',14)

figure(5)
plot(x,z)
title('Phase Space, z vs. x','FontSize',14)
xlabel('x(t) (unspecified units)','FontSize',14)
ylabel('z(t) (unspecified units)','FontSize',14)

figure(6)
plot(y,z)
title('Phase Space, z vs. y','FontSize',14)
xlabel('y(t) (unspecified units)','FontSize',14)
ylabel('z(t) (unspecified units)','FontSize',14)

%% 4c.

figure(9)
plot3(x,y,z,'b')
title('x(t) vs. y(t) vs. z(t)','FontSize',14)
xlabel('x(t)','FontSize',14)
ylabel('y(t)','FontSize',14)
zlabel('z(t)','FontSize',14)

%% 4d.

sigma2=10; b2=8/3; r2=28;
dt2=0.001;
t2=0:dt2:40;
x2(1)=1; y2(1)=1; z2(1)=20.01;
N2=length(t2);

for j=1:N2-1
    x_m2=x2(j) + (sigma2*(y2(j)-x2(j)))*(dt2/2);
    y_m2=y2(j) + (r2*x2(j) - y2(j) - x2(j)*z2(j))*(dt2/2);
    z_m2=z2(j) + ((x2(j)*y2(j) - b2*z2(j)))*(dt2/2);
    x2(j+1)=x2(j) + (sigma2*(y_m2 - x_m2))*dt2;
    y2(j+1)=y2(j) + (r2*x_m2 - y_m2 - x_m2*z_m2)*dt2;
    z2(j+1)=z2(j) + (x_m2*y_m2 - b2*z_m2)*dt2;
end

figure(7)
plot(t2,x2,'b')
title('x(t), finite difference method','FontSize',14)
xlabel('time (seconds)','FontSize',14)
ylabel('x(t) (unspecified units)','FontSize',14)

%% 4e.

dist_x=x2-x;
figure(8)
plot(t2,dist_x,'Color',[0.5 0 0.5])
title('distance in x(t) created by z\_o=20 -> 20.01','FontSize',14)
xlabel('time (seconds)','FontSize',14)
ylabel('distance (unspecified units)','FontSize',14)

%% Problem 5.
% Analysis of data of the oscillation of beams in the SLAC PEP-II B-meson
% factory.  Examination of the causes for the oscillation through a fast
% Fourier transform.
%%

delimiterIn=' ';
headerLinesIn=3;
file1='BPMy-x1.dat';
file2='BPMy-x2.dat';
D1=importdata(file1,delimiterIn,headerLinesIn);
D2=importdata(file2,delimiterIn,headerLinesIn);
t1=D1.data(:,2); t2=D2.data(:,2);
y1=D1.data(:,3); y2=D2.data(:,3);
x1=D1.data(:,4); x2=D2.data(:,4);

figure(1)
plot(t1,y1,'c')
title('BPM data from PEP-II, Water Pumps On, y-data','FontSize',14)
xlabel('Time (seconds)','FontSize',14)
ylabel('y (microns)','FontSize',14)

figure(2)
plot(t1,x1,'c')
title('BPM data from PEP-II, Water Pumps On, x-data','FontSize',14)
xlabel('Time (seconds)','FontSize',14)
ylabel('x (microns)','FontSize',14)

figure(3)
plot(t2,y2,'c')
title('BPM data from PEP-II, Water Pumps Off, y-data','FontSize',14)
xlabel('Time (seconds)','FontSize',14)
ylabel('y (microns)','FontSize',14)

figure(4)
plot(t2,x2,'c')
title('BPM data from PEP-II, Water Pumps Off, x-data','FontSize',14)
xlabel('Time (seconds)','FontSize',14)
ylabel('x (microns)','FontSize',14)

ffty1=fft(y1); fftx1=fft(x1);
ffty2=fft(y2); fftx2=fft(x2);

Py1=abs(ffty1).^2; Px1=abs(fftx1).^2;
Py2=abs(ffty2).^2; Px2=abs(fftx2).^2;

N1=length(t1); N2=length(t2);

dt1=(t1(end)-t1(1))/N1; dt2=(t2(end)-t2(1))/N2;

fy1=((1/dt1)*(1:N1)/N1); fx1=((1/dt1)*(1:N1)/N1);
fy2=((1/dt2)*(1:N2)/N2); fx2=((1/dt2)*(1:N2)/N2);

figure(5)
plot(fy1,Py1,'r')
xlim([0 2e3])
title('Power vs. frequency, y-data, Pumps On','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (\omega)','FontSize',14)

figure(6)
plot(fx1,Px1,'r')
xlim([0 2e3])
title('Power vs. frequency, x-data, Pumps On','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (\omega)','FontSize',14)

figure(7)
plot(fy2,Py2,'r')
xlim([0 2e3])
title('Power vs. frequency, y-data, Pumps Off','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (\omega)','FontSize',14)

figure(8)
plot(fx2,Px2,'r')
xlim([0 2e3])
title('Power vs. frequency, x-data, Pumps Off','FontSize',14)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('Power (\omega)','FontSize',14)

%% Were the water pumps causing low frequency oscillation?
% Yes!  Looking at our power graphs, we can see peaks at exactly 120 Hz for
% the data set we collected with the pumps turned on.  This corresponds to
% the oscillation of a device with an AC power source.  When looking at the
% data with the pumps off, the peak at 120 Hz disappears and only the
% natural frequency of the beams remain.
%% The "tunes" of the beam, natural frequencies, fx and fy
% For the x-data, the natural frequency is 843 Hz.  For the y-data, it's
% 571 Hz.  There is coupling present when the pumps are turned off, but
% with the pumps on, it's not very noticable, if it's taking place at all.
% The ratio between the strength of the natural frequency and the coupling
% frequency seen in the other plane seems to be on the order of 10 to 1.