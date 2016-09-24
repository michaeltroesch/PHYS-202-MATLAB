function y=sawtooth(n,xv)
y = zeros(1,length(xv));
k=1;
while k <= n
    if mod(k,2) == 1
        y = y+(2*(sin(k.*xv))./k);
    else
        y = y+(-2*(sin(k.*xv))./k);
    end
    k=k+1;
end
end