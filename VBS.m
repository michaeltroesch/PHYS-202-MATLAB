filename='VBS.txt';
delimiterIn=' ';
headerLinesIn=1;
A=importdata(filename,delimiterIn,headerLinesIn);
turns=A.data(:,1);
beam_size=A.data(:,2);
beam_error=A.data(:,3); %Importing data and setting to my own named vars
N=length(turns);
vbs_mean=0;
for i = 1:N %Calculating the mean of vertical beam size
    vbs_mean=vbs_mean+(beam_size(i));
end
vbs_mean=vbs_mean*(1/N);
sig_vbs=0;
for i = 1:N %Calculating the standard deviation of our beam size measures
    sig_vbs=sig_vbs+(beam_size(i)-vbs_mean)^2;
end
sig_vbs=sig_vbs*(1/(N-1));
str1=strcat('mean:  ', num2str(vbs_mean),' microns');
str2=strcat('\sigma:  ',num2str(sig_vbs));
h=errorbar(turns,beam_size,beam_error,'*:');
h.LineStyle = 'none';
ylabel('Beam Size (microns)','FontSize',16)
xlabel('Turn Number','FontSize',16)
xlim([-10 1010])
title('Beam Size vs. Turn # with error bars','FontSize',18)
text(50,33,str1,'FontSize',16)
text(50,32,str2,'FontSize',16)