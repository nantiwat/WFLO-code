function N=strtostr(n)
N=[];
for i=1:length(n)
    if n(i)~= ' '
        N=[N n(i)];
    end
end

