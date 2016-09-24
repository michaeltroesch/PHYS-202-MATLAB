%% Problem 1.
% Analysis of the time-constant of an RC circuit.  Parts a, b, and c are
% described below.
t=[0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30];
V=[0, 5.85, 10.3, 11.9, 14.5, 15, 16, 16.2, 16.3, 16.5, 16.55];
sigV=[0.15, 0.25, 0.35, 0.29, 0.21, 0.25, 0.55, 0.24, 0.32, 0.35, 0.31];
%Our recorded data for the circuit
%%
% 1a. Plotting our data with error bars.
errorbar(t,V,sigV,'b*','LineStyle','none','MarkerSize',3);
title('Capacitor Potential vs. Time with Error bars','FontSize',14)
ylabel('Capacitor Potential (V)','FontSize',14)
xlabel('Time (s)','FontSize',14)
xlim([0 32]);
%%
% 1b. Our X^2/DOF value of 1.192 suggests a good fit!
N=length(t); nu=2; dof=N-nu;
w=1./(sigV.^2);
[tdata,Vdata,weights]=prepareCurveData(t,V,w);
ft=fittype('A*(1-exp(B*t))','independent','t');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Weights = weights;
opts.StartPoint=[5 0];
[fitresult,~]=fit(tdata,Vdata,ft,opts);
A=fitresult.A;
B=fitresult.B;
chisqr_dof=(1/dof)*sum(weights.*(Vdata - A*(1-exp(B*tdata))).^2);
%%
% 1c. Plotting our data along with the calculated fitline.
plot(fitresult,'k',tdata,Vdata);
hold on
legend({'Voltage data','Nonlinear Fit'},'FontSize',12','Location','SouthEast');
h=errorbar(t,V,sigV,'b ','LineStyle','none');
set(h,'MarkerSize',0.5);
title('Capacitor Potential vs. Time with Error bars','FontSize',14)
ylabel('Capacitor Potential (V)','FontSize',14)
xlabel('Time (s)','FontSize',14)
xlim([0 32]);
sig_A=0.0385;
sig_B=0.01285;
str1=strcat('\chi^2 / DOF = ',num2str(chisqr_dof));
str2=strcat('V_\circ = ',num2str(A),' +/- ',num2str(sig_A),' V');
str3=strcat('\tau = ', num2str(-1/B),' +/- ',num2str(-1/B-(-1/(B-sig_B))),' s');
text(14,10,str1);
text(14,8.8,str2);
text(14,8.2,str3);
hold off
%% Problem 2.
% Description of radioactive decay of silver.  Parts a, b, and c are
% outlined in the comments.
%%
% 2a. Importing data with importdata function. Plotting data with errorbar.
headerLinesIn=1;
delimiterIn=' ';
fn='SilverDecay.dat';
SilverDecay=importdata(fn,delimiterIn,headerLinesIn);
t=SilverDecay.data(:,1);
counts=SilverDecay.data(:,2);
error_counts=SilverDecay.data(:,3);
errorbar(t,counts,error_counts,'b*','LineStyle','none','MarkerSize',3)
xlabel('Time (s)','FontSize',14);
xlim([-1 451]);
ylabel('Counts of Silver','FontSize',14);
title('Counts of Silver vs. Time','FontSize',14);
%%
% 2b. Determining all of our fit-line and statistical analysis
N=length(t); nu=2; dof=N-nu;
w=(1./(error_counts.^2));
[tdata,countsdata,weights]=prepareCurveData(t,counts,w);
ft=fittype('A*exp(B*t)','independent','t');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Weights = weights;
opts.StartPoint=[60 0];
[fitresult,~]=fit(tdata,countsdata,ft,opts);
A=fitresult.A;
B=fitresult.B;
chisqr_dof=(1/dof)*sum(weights.*(countsdata - A*exp(B*tdata)).^2);
sig_A=3.08;
sig_B=0.000260;
%Our value of Chi-Square / DOF is close to 1 and indicates a good fit.
%%
% 2c. Plotting the fit line and statistical data.
plot(fitresult,'k',tdata,countsdata)
hold on
errorbar(t,counts,error_counts,'c ','LineStyle','none')
xlabel('Time (s)','FontSize',14);
ylabel('Counts of Silver','FontSize',14);
xlim([-1 451]);
title('Counts of Silver vs. Time','FontSize',14);
legend({'Counts data','Nonlinear Fit'},'FontSize',12','Location','NorthEast');
str1=strcat('\chi^2 / DOF = ',num2str(chisqr_dof));
str2=strcat('I_\circ = ',num2str(A),' +/- ',num2str(sig_A),' counts');
str3=strcat('\tau = ', num2str(-1/B),' +/- ',num2str(-1/B-(-1/(B-sig_B))),' s');
text(200,65,str1)
text(200,60,str2)
text(200,57,str3)
hold off
%% Problem 3.
% Reading in multiple data sets of positron beam sizes from the Stanford
% Linear Collider Damping Ring.
%%
% 3a. Using MATLAB's importdata function to read in the 8 .dat files
numfiles=8;
delimiterIn=' ';
headerLinesIn=3;

for k = 1:numfiles
  myfilename = sprintf('Data%d.dat', k);
  SLC{k} = importdata(myfilename,delimiterIn,headerLinesIn);
end

for k = 1:length(SLC)
    pixelData{k} = SLC{1,k}.data(:,:);
end

x=pixelData{k}(:,1);
%we can set the x-axis because it is constant for this dataset

%%
% 3b. This code goes through all the beam data.
% It calculates X^2/DOF.  The values are all around 1.
for k = 1:length(pixelData)
    for j = 2:11
    y=pixelData{k}(:,j);
    sig_pixel=sqrt(y);
    N=length(x); nu=4; dof=N-nu;
    w=1./sig_pixel.^2;
    [xdata,ydata,weights]=prepareCurveData(x,y,w);
    ft=fittype('A + B*exp((-1/2)*((x-C)/D).^2)');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint=[100 45 100 90];
    if k == 5
        opts.StartPoint=[-500 80 80 80];
    end
    if k == 6
        opts.StartPoint=[-2000 100 80 70];
    end
    if k == 8
        opts.StartPoint=[-1000 45 100 150];
    end
    %used trial and error to determine start points.
    [fitresult,~]=fit(xdata,ydata,ft,opts);
    A=fitresult.A;
    B=fitresult.B;
    C=fitresult.C;
    D=fitresult.D;
    VBS_array(k,j-1)=abs(D);
    chisqr_dof=(1/(dof))*sum(weights.*(ydata - (A + B*exp(-1/2*((xdata-C)/D).^2))).^2);
    I=0.25*k;
    Profile=j-1;
    title_string=strcat('SLC Beam Data, I= ', num2str(I),'mA, Profile #',num2str(Profile));
    errorbar(xdata,ydata,sig_pixel)
    hold on
    plot(fitresult,'y',xdata,ydata)
    legend({'Error Bars','Data Points','Fitted Curve'},'Location','NorthEast')
    title(title_string,'FontSize',14)
    xlabel('Pixel # (microns)','FontSize',14)
    ylabel('Pixel Intensity (P)','FontSize',14)
    xlim([0 512]);
    ylim=get(gca,'ylim');
    str1=strcat('\chi^2/DOF= ', num2str(chisqr_dof));
    text(10,ylim(2)-420,str1)
    hold off
    end
end

% 3c. Plots for 0.25 mA and 1.75 mA.
for k = 1
    for j = 4
    y=pixelData{k}(:,j);
    sig_pixel=sqrt(y);
    N=length(x); nu=4; dof=N-nu;
    w=1./sig_pixel.^2;
    [xdata,ydata,weights]=prepareCurveData(x,y,w);
    ft=fittype('A + B*exp((-1/2)*((x-C)/D).^2)');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint=[100 45 100 90]; 
    
    if k == 5
        opts.StartPoint=[-500 80 80 80];
    end
    
    if k == 6
        opts.StartPoint=[-2000 100 80 70];
    end
    
    if k == 8
        opts.StartPoint=[-1000 45 100 150];
    end
    
    %used trial and error to determine start points.
    
    [fitresult,~]=fit(xdata,ydata,ft,opts);
    A=fitresult.A;
    B=fitresult.B;
    C=fitresult.C;
    D=fitresult.D;
    VBS_array(k,j-1)=abs(D);
    chisqr_dof=(1/(dof))*sum(weights.*(ydata - (A + B*exp(-1/2*((xdata-C)/D).^2))).^2);
    I=0.25*k;
    Profile=j-1;
    title_string=strcat('SLC Beam Data, I= ', num2str(I),...
    'mA, Profile #',num2str(Profile));
    errorbar(xdata,ydata,sig_pixel)
    hold on
    plot(fitresult,'y',xdata,ydata)
    legend({'Error Bars','Data Points','Fitted Curve'},'Location','NorthEast')
    title(title_string,'FontSize',14)
    xlabel('Pixel # (microns)','FontSize',14)
    ylabel('Pixel Intensity (P)','FontSize',14)
    xlim([0 512]);
    ylim=get(gca,'ylim');
    str1=strcat('\chi^2/DOF= ', num2str(chisqr_dof));
    text(10,ylim(2)-420,str1)
    hold off
    end
end
%%
for k = 7
    for j = 7
    y=pixelData{k}(:,j);
    sig_pixel=sqrt(y);
    N=length(x); nu=4; dof=N-nu;
    w=1./sig_pixel.^2;
    [xdata,ydata,weights]=prepareCurveData(x,y,w);
    ft=fittype('A + B*exp((-1/2)*((x-C)/D).^2)');
    opts=fitoptions('Method','NonlinearLeastSquares');
    opts.StartPoint=[100 45 100 90]; 
    
    if k == 5
        opts.StartPoint=[-500 80 80 80];
    end
    
    if k == 6
        opts.StartPoint=[-2000 100 80 70];
    end
    
    if k == 8
        opts.StartPoint=[-1000 45 100 150];
    end
    
    %used trial and error to determine start points.
    
    [fitresult,~]=fit(xdata,ydata,ft,opts);
    A=fitresult.A;
    B=fitresult.B;
    C=fitresult.C;
    D=fitresult.D;
    VBS_array(k,j-1)=abs(D);
    chisqr_dof=(1/(dof))*sum(weights.*(ydata - (A + B*exp(-1/2*((xdata-C)/D).^2))).^2);
    I=0.25*k;
    Profile=j-1;
    title_string=strcat('SLC Beam Data, I= ', num2str(I),'mA, Profile #',num2str(Profile));
    errorbar(xdata,ydata,sig_pixel)
    hold on
    plot(fitresult,'y',xdata,ydata)
    legend({'Error Bars','Data Points','Fitted Curve'},'Location','NorthEast')
    title(title_string,'FontSize',14)
    xlabel('Pixel # (microns)','FontSize',14)
    ylabel('Pixel Intensity (P)','FontSize',14)
    xlim([0 512]);
    ylim=get(gca,'ylim');
    str1=strcat('\chi^2/DOF= ', num2str(chisqr_dof));
    text(10,ylim(2)-420,str1)
    hold off
    end
end
%%
% 3c. Table of mean VBS and uncertainty.
for i=1:8
    VBS_mean(i)=sum(VBS_array(i,:))/10;
end
VBS_uncertainty=sqrt(VBS_mean)./sqrt(10);
I_array=[0.25 0.5 0.75 1 1.25 1.5 1.75 2];

I=transpose(I_array);
mean_VBS=transpose(VBS_mean);
error_VBS=transpose(VBS_uncertainty);
T=table(I,mean_VBS,error_VBS);
disp(T)
%%
% 3d. Plotting mean VBS w/ uncertainty.
errorbar(I_array,VBS_mean,VBS_uncertainty,'b*','MarkerSize',3);
title('Mean VBS with Error Bars','FontSize',14);
ylabel('Vertical Beam Size (microns)','FontSize',14);
xlabel('Current, I (mA)','FontSize',14);
xlim([0.2 2.05]);
hold on
%%
% 3e. Fitting the plot with a 2nd-order polynomial.
% The value of X^2/DOF suggests a bad fit!?  Maybe a cubic function is
% needed... we'd have to do more research.
N=length(VBS_mean); nu=3; dof=N-nu;
w=1./(VBS_uncertainty).^2;
[Idata,VBSdata,weights]=prepareCurveData(I_array,VBS_mean,w);
ft=fittype('A+B*I+C*I^2','independent','I');
opts=fitoptions('Method','NonlinearLeastSquares');
opts.Weights = weights;
opts.StartPoint=[7 7 7];
[fitresult,gof]=fit(Idata,VBSdata,ft,opts);
A=fitresult.A;
B=fitresult.B;
C=fitresult.C;
chisqr_dof=(1/dof)*sum(weights.*(VBSdata - (A + B.*Idata + C.*Idata.^2)).^2);
plot(fitresult,'k')
str1=strcat('\chi^2/DOF =',num2str(chisqr_dof));
sig_A=13.5;
sig_B=31.9;
sig_C=15.2;
str2=strcat('y= ',num2str(A),'+/-',num2str(sig_A),'+',num2str(B),'+/-',...
num2str(sig_B),'*I+', num2str(C),'+/-',num2str(sig_C),'*I^2 microns');
text(.22,85,str2)
text(.22,75,str1)
legend({'VBS Data','Fitted Curve'},'Location','NorthWest')
title('Fitting VBS to 2nd-order Polynomial','FontSize',14)
xlabel('Current, I (mA)','FontSize',14)
ylabel('Mean Vertical Beam Size (microns)','FontSize',14)
clear plots;
clear variables;
%% Problem 4.
% Estimating the value of pi by a Monte Carlo method.
%%
% 4a.
% Here's how you show the relation between pi, the area of the circle, and
% the area of a square it is inscribed in, defining g as a ratio between the
% area of a circle and the area of a square:
%
% $$g= \frac{A_{circle}}{A_{square}}$$
%
% $$\frac{\pi r^2}{4 r^2}=\frac{\pi } {4}$$
%
% $$\pi = 4g$$
%
% We can use this relation to estimate pi with our virtual dartboard.
% The ratio between marks landing inside the circle vs. the total
% number of marks (those inside the square) will approximate the ratio
% of the two areas, multiplying by 4 will give us an estimate of the value
% of pi.
%%
% 4b. The plot ends up looking random, seems right!
hold off;
clear all;
n=100;
for i=1:n
    x(i)= -1 + 2*rand();
    y(i)= -1 + 2*rand();
end
plot(x,y,' d','MarkerSize',3)
title('Plot of random x/y vectors','FontSize',14)
%%
% 4c. Determining what marks lay inside the circle.
Xc=[];
Yc=[];
n=100;
for i=1:n
    x(i)= -1 + 2*rand();
    y(i)= -1 + 2*rand();
    if (x(i)^2 + y(i)^2) > 1
        Xc = [Xc, x(i)];
        Yc = [Yc, y(i)];
    end
end
plot(sin(0:0.1:2*pi+0.1),cos(0:0.1:2*pi+0.1))
title('Virtual Dartboard','FontSize',14)
hold on
plot(x,y,' d','MarkerSize',3)
hold on
plot(Xc,Yc,' *','MarkerSize',6)
hold off
%% 
% 4d. Estimating pi using the relation from 4a.
num_inside=n - length(Xc);
estimate_ratio=num_inside/n;
estimate_pi=4*estimate_ratio
%% 
pi_estimates=[];
n=10;
n_iter=[];
while n <= 1000
Xc=[];
Yc=[];
for i=1:n
    x(i)= -1 + 2*rand();
    y(i)= -1 + 2*rand();
    if (x(i)^2 + y(i)^2) > 1
        Xc = [Xc, x(i)];
        Yc = [Yc, y(i)];
    end
end
num_inside=n - length(Xc);
estimate_ratio=num_inside/n;
estimate_pi=4*estimate_ratio;
pi_estimates = [pi_estimates, estimate_pi];
n_iter=[n_iter, n];
n = n + 10;
end
semilogx(n_iter,pi_estimates,'b*','LineStyle','none')
title('Estimates of \pi from n=10 to 10^3, step-size: 10','FontSize',14);
xlabel('n (Number of Trials)','FontSize',14)
ylabel('Result of \pi estimate','FontSize',14)
hold on
line([10 1000],[pi pi],'LineWidth',1)
hold off
%% Problem 5.
% Statistical analysis of dice-rolls and dealing with the random number
% functions built in to MATLAB.
%%
% 5a. Rolling a single die 10^1 - 10^6 times.
for i=1:10
    n_1(i)=ceil(6*rand());
end

for i=1:100
    n_2(i)=ceil(6*rand());
end

for i=1:1000
    n_3(i)=ceil(6*rand());
end

for i=1:10000
    n_4(i)=ceil(6*rand());
end

for i=1:100000
    n_5(i)=ceil(6*rand());
end

for i=1:1000000
    n_6(i)=ceil(6*rand());
end

bardata=zeros(1,6);
for j=1:length(n_4)
    if n_4(j)==1
        bardata(1) = bardata(1) + 1;
    elseif n_4(j)==2
        bardata(2) = bardata(2) + 1;
    elseif n_4(j)==3
        bardata(3) = bardata(3) + 1;
    elseif n_4(j)==4
        bardata(4) = bardata(4) + 1;
    elseif n_4(j)==5
        bardata(5) = bardata(5) + 1;
    else
        bardata(6) = bardata(6) + 1;
    end
end

bar(bardata)
title('10k trials, one six-sided die','FontSize',14)
ylabel('Number of rolls','FontSize',14)
xlabel('Result of roll','FontSize',14)
%%
% 5a. Our table of means and standard deviations.
rolls=[10;10^2;10^3;10^4;10^5;10^6];
mean=[mean(n_1);mean(n_2);mean(n_3);mean(n_4);mean(n_5);mean(n_6)];
sigma=[std(n_1);std(n_2);std(n_3);std(n_4);std(n_5);std(n_6)];
T=table(rolls,mean,sigma);
disp(T)
%%
% 5b. Rolling two dice 10^1 - 10^6 times.
clear variables;
for i=1:10
    n_1(i)=ceil(6*rand())+ceil(6*rand());
end

for i=1:100
    n_2(i)=ceil(6*rand())+ceil(6*rand());
end

for i=1:1000
    n_3(i)=ceil(6*rand())+ceil(6*rand());
end

for i=1:10000
    n_4(i)=ceil(6*rand())+ceil(6*rand());
end

for i=1:100000
    n_5(i)=ceil(6*rand())+ceil(6*rand());
end

for i=1:1000000
    n_6(i)=ceil(6*rand())+ceil(6*rand());
end

bardata2=zeros(1,12);
for j=1:length(n_4)
    if n_4(j)==2
        bardata2(2) = bardata2(2) + 1;
    elseif n_4(j)==3
        bardata2(3) = bardata2(3) + 1;
    elseif n_4(j)==4
        bardata2(4) = bardata2(4) + 1;
    elseif n_4(j)==5
        bardata2(5) = bardata2(5) + 1;
    elseif n_4(j)==6
        bardata2(6) = bardata2(6) + 1;
    elseif n_4(j)==7
        bardata2(7) = bardata2(7) + 1;
    elseif n_4(j)==8
        bardata2(8) = bardata2(8) + 1;
    elseif n_4(j)==9
        bardata2(9) = bardata2(9) + 1;
    elseif n_4(j)==10
        bardata2(10) = bardata2(10) + 1;
    elseif n_4(j)==11
        bardata2(11) = bardata2(11) + 1;
    else
        bardata2(12) = bardata2(12) + 1;
    end
end

bar(bardata2);
xlim([1.5 12.5])
title('10k trials, rolling two six-sided die','FontSize',14)
ylabel('Number of rolls','FontSize',14)
xlabel('Result of roll','FontSize',14)

% 5b. This if-statement handles the logic of taking a bet.
% A 6,7, or 8 doubles up, but you get nothing otherwise.
% That's why it doesn't make sense to take the bet.
if (bardata2(6) + bardata2(7) + bardata2(8))/10000 > .5
    disp('Take that bet!')
else
    disp('No thank you, I like my money.')
end