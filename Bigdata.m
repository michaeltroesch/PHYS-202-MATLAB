filename='Bigdata.dat';
delimiterIn=' ';
headerLinesIn=1;
B=importdata(filename,delimiterIn,headerLinesIn);
B_field=B.data(:,1);
B_error=B.data(:,2);
current=B.data(:,3); %Importing data and setting to my own named vars
errorbar(current,B_field,B_error,'c*');
hold on
title('Magnetic field vs. Current','FontSize',18)
ylabel('Magnetic field (T)','FontSize',16)
xlabel('Current (A)','FontSize',16)
N=length(current);
w=1./(B_error.^2);
del=sum(w)*sum(w.*(current.^2))-(sum(w.*current))^2;
a=(sum(w.*(current.^2))*sum(w.*B_field)-sum(w.*current)*sum(w.*current.*B_field))/del;
b=(sum(w)*sum(w.*current.*B_field)-sum(w.*current)*sum(w.*B_field))/del;
fit_line=a+b*current;
chisqr_dof=(1/(N-2))*sum(w.*(B_field - a - b*current).^2);
sig_a=sqrt((sum(w.*(current.^2)))/del);
sig_b=sqrt(sum(w)/del);
str1=strcat('B=(',num2str(b),'+/-',num2str(sig_b),')*x + (',num2str(a),'+/-',num2str(sig_a),') T');
str2=strcat('\mu naught=',num2str(b*2*pi),'+/-',num2str(sig_b*2*pi),' T*m/A');
str3=strcat('\chi^2/DOF=',num2str(chisqr_dof));
plot(current,fit_line,'k');
text(-3,5e-6,str1,'FontSize',12)
text(-3,4.7e-6,str2,'FontSize',12)
text(-3,4.4e-6,str3,'FontSize',12)
legend({'B-field vs. Current','Linear Fit'}, 'FontSize',16,'Location', 'SouthEast');