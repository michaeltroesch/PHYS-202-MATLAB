function result=enxconverge()
x=0:20;
k=0;
result=0;
while k < length(x)
    if round(enxsum(k),7,'significant') == round(exp(1)/(exp(1)-1),7,'significant')
    if round(enxsum(k),7,'significant') == round(enxsum(k+1),7,'significant')
        result = k;
        k = k + 20;
    end
    end
        k = k + 1;
end
end

