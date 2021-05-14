function SODE(fun,foutput,nloop,nsol,nvar,nbit,a,b)
% Hybrid differential evolution & Nelder and Mead method
% fun = m-file for objective function
% foutput = mat-flie for output results
% nloop = no. of generations
% nsol = population size
% nvar = no. of design variables
% nbit = no. of binary bits for each design variable (in cases of the methdo using binary strings)
% a = vector of lower limits
% b = vector of upper limits

% predefined parameters
pc=0.7;% crossover probability
F=0.5;% Scaling factor for differentail evolution (DE) operator
CR=0.8;% probability of chosing element from offspring in crossover

rand('state',sum(100*clock));

[bin0,x0,fp0,f0,g0] = ea_initial(fun,nbit,nvar,nsol);
[fpmin,npmin]=min(fp0);
fpminhist=min(fp0);fpavghist=mean(fp0);fpmaxhist=max(fp0);
fhist=f0(npmin);ghist=g0(:,npmin);

iter=0;neval=0;
while neval < nloop*nsol
    iter=iter+1;
    [x1,fp1,f1,g1]=de_reproduct(fun,x0,fp0,f0,g0,a,b,pc,F,CR);% DE operator
    neval=neval+length(f1);
        
    for i=1:nsol
        if fp1(i) <= fp0(i)
            x0(:,i)=x1(:,i);
            fp0(i)=fp1(i);
            f0(i)=f1(i);
            g0(:,i)=g1(:,i);
        end
    end
    [fminI(iter),nminI]=min(fp0);gminI(:,iter)=g0(:,nminI);
    
    [fpmin,npmin]=min(fp0);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp0)];
    fpmaxhist=[fpmaxhist max(fp0)];
    fhist=[fhist f0(npmin)];ghist=[ghist g0(:,npmin)];
end
[fpmin,nmin]=min(fp0);
xmin=x0(:,nmin);fmin=f0(:,nmin);gmin=g0(:,nmin);
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
% figure(1),clf,hold on
% plot((1:length(besthist))*nloop*nsol/length(besthist),besthist,'r')
% plot((1:length(besthist))*nloop*nsol/length(besthist),avghist,'b')
% plot((1:length(besthist))*nloop*nsol/length(besthist),worsthist,'g')
%%%%%%%%%% Sub-programs %%%%%%%%%%%%%%
function [x1,fp1,f1,g1] = de_reproduct(fun,x0,fp0,f0,g0,a,b,pc,F,CR)
% best/2/bin strategy
[m0,n0]=size(x0);
[fbest,nmin]=min(fp0);
xbest=x0(:,nmin);
% DE point generation for n1 parents (x1)
% x0 = population from the previous generation
% f0 = corresponding objective values

for i=1:n0
    nr=randperm(n0);
    i1=nr(1);i2=nr(2);i3=nr(3);i4=nr(4);
    xr1=x0(:,i1);% Randomly seletced individual 1
    xr2=x0(:,i2);% Randomly seletced individual 2
    xr3=x0(:,i3);% Randomly seletced individual 3
    xr4=x0(:,i4);% Randomly seletced individual 4

    Ci=xbest+F*(xr1+xr2-xr3-xr4);

    for j=1:m0
        if Ci(j) < a(j)
            Ci(j)=a(j);
        elseif Ci(j) > b(j)
            Ci(j)=b(j);
        end
    end
    
    if rand < pc
        for j=1:m0 % binary crossover
            if rand < CR
                Vi(j,1)=Ci(j);
            else
                Vi(j,1)=x0(j,i);
            end
        end
        Ci=Vi;
    end
    [fpCi,fCi,gCi]=feval(fun,Ci);
    if fpCi <= fp0(i)
        x1(:,i)=Ci;
        fp1(i)=fpCi;
        f1(i)=fCi;
        g1(:,i)=gCi;
    else
        x1(:,i)=x0(:,i);
        fp1(i)=fp0(i);
        f1(i)=f0(i);
        g1(:,i)=g0(:,i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

