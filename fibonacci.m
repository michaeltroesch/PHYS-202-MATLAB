function res=fibonacci(x)
res = [];
start = 1;
k=1;
while k <= x
    if length(res) <= 2
        res(k) = start;
    end
    if length(res) > 2
        res(k) = res(k-1)+ res(k-2);
    end
	k = k + 1;
end
end