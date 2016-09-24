function yq=m_linear(x,y,xq)
    %x,y: x and y data points to be interpolated from
    %xq: points to evaluate via interpolation
    %yq: returned value for line @ points given by xq  
    d=diff(y)./diff(x);
    n=length(x);
    k=ones(size(xq));
    for i=1:n-1
        if x(i+1) > x(i)
            k(x(i) <= xq) = i;
        else
            k(x(i) >= xq) = i;
        end
    end
    s=xq-x(k);
    yq=y(k) + s.*d(k);
end
