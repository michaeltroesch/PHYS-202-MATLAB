function s=enxsum(num)
s=0;
k=0;
while k <= num
    s = s + (1/exp(k));
    k = k + 1;
end
end