function N=numtostr(n)
if n<=9
    N=['0' num2str(n)];
else
    N=num2str(n);
end