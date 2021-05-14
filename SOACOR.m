function SOACOR(fun,foutput,nloop,nsol,nvar,nbit,a,b)
% Real code ant colony optimisation

% fun = m-file for objective function
% foutput = mat-flie for output results
% nloop = no. of generations
% nsol = population size
% nvar = no. of design variables
% nbit = no. of binary bits for each design variable (in cases of the methdo using binary strings)
% a = vector of lower limits
% b = vector of upper limits

AchSize=nsol;
q=0.2;
xi=1.0;

rand('state',sum(100*clock));
randn('state',sum(100*clock));

[bin0,x0,fp0,f0,g0] = ea_initial(fun,nbit,nvar,nsol);
[fpmin,npmin]=min(fp0);
fpminhist=min(fp0);fpavghist=mean(fp0);fpmaxhist=max(fp0);
fhist=f0(npmin);ghist=g0(:,npmin);
    
iter=0;
while iter < nloop
    iter=iter+1;
    [ff,l]=sort(fp0);
    xmin=x0(:,l(1));fmin=fp0(l(1));
    for i=1:nsol
        for j=1:nvar
            meanval=x0(j,:);
            w=(1/nsol/q/sqrt(2*pi))*exp(-(l-1).^2/2/q^2/nsol^2);
            prob=w/sum(w);
            RW=[0 cumsum(prob)];% establish a roulette wheel
            sp=rand;% uniform random number [0,1]
            kk=0;sw=1;
            while sw==1
                kk=kk+1;
                if RW(kk)<=sp & sp < RW(kk+1)
                    sw=0;
                    isl=kk;
                elseif kk==nsol
                    sw=0;
                    isl=nsol;
                end
            end
            xl=x0(:,isl);
            sig=0;
            for e=1:nsol
                sig=sig+abs(xl(j)-x0(j,e));
            end
            sigma=xi*sig/(nsol-1);
            x(j,i)=meanval(isl)+sigma*randn;
            if x(j,i)<a(j);x(j,i)=a(j);end
            if x(j,i)>b(j);x(j,i)=b(j);end
        end
        [fp(i),f(i),g(:,i)]=feval(fun,x(:,i));
    end
    
    % selection
    XX=[x0 x];FP=[fp0 fp];F=[f0 f];G=[g0 g];
    [FFP,ns]=sort(FP);F=F(ns);XX=XX(:,ns);FP=FP(ns);G=G(:,ns);
    f0=F(1:nsol);x0=XX(:,1:nsol);fp0=FP(1:nsol);g0=G(:,1:nsol);
    nrp=randperm(nsol);
    f0=f0(nrp);x0=x0(:,nrp);fp0=fp0(nrp);g0=g0(:,nrp);
    
    [fpmin,npmin]=min(fp0);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp0)];
    fpmaxhist=[fpmaxhist max(fp0)];
    fhist=[fhist f0(npmin)];ghist=[ghist g0(:,npmin)];
end
[fpmin,nmin]=min(fp0);
xmin=x0(:,nmin);fmin=f0(nmin);
gmin=g0(:,nmin);
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%%%%%%% Sub-programs %%%%%%%%%%%%%%
function [bin,x,f] = aco_initial(fun,nvar,nbit,nsol,a,b)
%
% Randomly initiate the population, design variables
% nvar=no. of variables
% nbit is the number of cell in each variable
% nsol is a number of gene
%
bin = round(rand(nvar*nbit,nsol));
for i=1:nsol
    for j=1:nvar
        binn=bin((j-1)*nbit+1:j*nbit,i);
        x(j,i)=bin2dec(binn,a(j),b(j));
    end
    f(i)=feval(fun,x(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=fdistant(x1,x2)
x12=x1-x2;
d=sqrt(sum(x12.*x12));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=fTrustReg(x,a,b)
for i=1:length(x)
    if x(i) < a(i)
        x(i)=a(i);
    elseif x(i) > b(i)
        x(i)=b(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=bin2dec(bin,a,b)
%
% Transformation from binary string to real number
% with lowr limit a and upper limit b
n=max(size(bin));
trans=cumprod(2*ones(size(bin)))/2;
real1=sum(bin.*trans);
x=a+(real1*(b-a))/(2^n-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
