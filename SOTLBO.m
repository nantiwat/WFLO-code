function SODE_NMM(fun,foutput,nloop,nsol,nvar,nbit,a,b)
% Learning-Teaching Based optimisation
% fun = m-file for objective function
% foutput = mat-flie for output results
% nloop = no. of generations
% nsol = population size
% nvar = no. of design variables
% nbit = no. of binary bits for each design variable (in cases of the methdo using binary strings)
% a = vector of lower limits
% b = vector of upper limits


rand('state',sum(100*clock));
randn('state',sum(100*clock));

[bin0,x0,fp0,f0,g0] = ea_initial(fun,nbit,nvar,nsol);
[fpmin,npmin]=min(fp0);
fpminhist=min(fp0);fpavghist=mean(fp0);fpmaxhist=max(fp0);
fhist=f0(npmin);ghist=g0(:,npmin);

iter=0;neval=0;
while neval < nloop*nsol
    iter=iter+1;
    [x1,fp1,f1,g1]=tlbo_teaching(fun,x0,fp0,f0,g0,a,b);% teaching phase
    neval=neval+length(f1);
    
    [x2,fp2,f2,g2]=tlbo_learning(fun,x1,fp1,f1,g1,a,b);% learning phase
    neval=neval+length(f2);
       
    x0=x2;fp0=fp2;f0=f2;g0=g2;
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
function [xbest,fpbest,fbest,gbest]=tlbo_selection(x0,fp0,f0,g0,xbest,fpbest,fbest,gbest);
[m,n]=size(x0);
x=[xbest x0];
fp=[fpbest fp0];
f=[fbest f0];
g=[gbest g0];
[fpbest,ns]=sort(fp);
xbest=x(:,ns);
fbest=f(ns);
gbest=g(:,ns);

xbest=xbest(:,1:n);
fpbest=fpbest(1:n);
fbest=fbest(1:n);
gbest=gbest(:,1:n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,fp1,f1,g1]=tlbo_teaching(fun,x0,fp0,f0,g0,a,b);% teaching phase
[m,n]=size(x0);
[fmin,nmin]=min(fp0);
Teacher=x0(:,nmin);
ClassAvg=mean(x0,2);

for i=1:n
    x1(:,i)=rand(m,1).*(Teacher-ceil(2*rand)*ClassAvg);
    x1(:,i)=fTrustReg(x1(:,i),a,b);
    [fp1(i),f1(i),g1(:,i)]=feval(fun,x1(:,i));
end
[x1,fp1,f1,g1]=tlbo_selection(x0,fp0,f0,g0,x1,fp1,f1,g1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,fp2,f2,g2]=tlbo_learning(fun,x1,fp1,f1,g1,a,b);% learning phase
[m,n]=size(x1);

for i=1:n
    isl=randperm(n);j1=isl(1);j2=isl(2);
    X1=x1(:,j1);X2=x1(:,j2);
    F1=fp1(:,j1);F2=fp1(:,j2);
    if F1<F2
        x2(:,i)=x1(:,i)+rand*(X1-X2);
    else
        x2(:,i)=x1(:,i)+rand*(X2-X1);
    end
    x2(:,i)=fTrustReg(x2(:,i),a,b);
    [fp2(i),f2(i),g2(:,i)]=feval(fun,x2(:,i));
end
[x2,fp2,f2,g2]=tlbo_selection(x1,fp1,f1,g1,x2,fp2,f2,g2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function x=fTrustReg(x,a,b)
for i=1:length(x)
    if x(i) < a(i)
        x(i)=a(i);
    elseif x(i) > b(i)
        x(i)=b(i);
    end
end
%%%%%%%%%%%%% end of file %%%%%%%%%%%%%%%

