function quadSolver(a,b,c)  %quadSolver: returns the roots of a quadratic equation
                            %a,b,c: coefficients of the quadratic equation
det=b^2-(4*a*c);            %Calculating the determinant
root1=(-b+sqrt(det))/(2*a);  %Calculating the roots
root2=(-b-sqrt(det))/(2*a);

if det > 0
    disp(['There are two distinct real roots: ',num2str(root1),', ',num2str(root2)])
elseif det < 0
    disp(['There are two distinct imaginary roots: ',num2str(root1), ', ',num2str(root2)])
else
    disp(['There is one real root: ',num2str(root1)])
end
