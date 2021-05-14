% -------------------------------------------------
% Citation details:
% Alireza Askarzadeh, Anovel metaheuristic method for solving constrained
% engineering optimization problems: Crow search algorithm, Computers &
% Structures, Vol. 169, 1-12, 2016.

% Programmed by Alireza Askarzadeh at Kerman Graduate %
% University of Advanced Technology (KGUT) %
% Date of programming: September 2015 %
% -------------------------------------------------
% This demo only implements a standard version of CSA for minimization of
% a standard test function (Sphere) on MATLAB 7.6.0 (R2008a).
% -------------------------------------------------
% Note:
% Due to the stochastic nature of meta-heuristc algorithms, different runs
% may lead to slightly different results.
% -------------------------------------------------

function SOCSA(fun,foutput,nloop,nsol,nvar,nbit,a,b)

pd=nvar; % Problem dimension (number of decision variables)
N=nsol; % Flock (population) size
AP=0.1; % Awareness probability
fl=2; % Flight length (fl)

[x l u]=init(N,pd,a(1),b(1)); % Function for initialization

xn=x;
[ft,f,g]=fitness(xn,N,pd,fun,a',b'); % Function for fitness evaluation
fp=ft;
[fpmin,npmin]=min(fp);
fpminhist=min(fp);fpavghist=mean(fp);fpmaxhist=max(fp);
fhist=f(npmin);ghist=g(:,npmin);

mem=x; % Memory initialization
fit_mem=ft; % Fitness of memory positions

tmax=nloop; % Maximum number of iterations (itermax)
for t=1:tmax

    num=ceil(N*rand(1,N)); % Generation of random candidate crows for following (chasing)
    for i=1:N
        if rand>AP
            xnew(i,:)= x(i,:)+fl*rand*(mem(num(i),:)-x(i,:)); % Generation of a new position for crow i (state 1)
        else
            for j=1:pd
                xnew(i,j)=l-(l-u)*rand; % Generation of a new position for crow i (state 2)
            end
        end
    end

    xn=xnew;
    
    [ft,f,g]=fitness(xn,N,pd,fun,a',b'); % Function for fitness evaluation of new solutions
    fp=ft;
    [fpmin,npmin]=min(fp);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp)];
    fpmaxhist=[fpmaxhist max(fp)];
    fhist=[fhist f(npmin)];ghist=[ghist g(:,npmin)];
    
    for i=1:N % Update position and memory
        if xnew(i,:)>=l & xnew(i,:)<=u
            x(i,:)=xnew(i,:); % Update position
            if ft(i)<fit_mem(i)
                mem(i,:)=xnew(i,:); % Update memory
                fit_mem(i)=ft(i);
            end
        end
    end

    ffit(t)=min(fit_mem); % Best found value until iteration t
    min(fit_mem);
end

ngbest=find(fit_mem== min(fit_mem));
g_best=mem(ngbest(1),:); % Solutin of the problem
xmin=g_best';
[fpmin,fmin,gmin]=feval(fun,xmin);
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x l u]=init(N,pd,l,u) % Function for initialization

%l=0; u=1; % Lower and upper bounds

for i=1:N % Generation of initial solutions (position of crows)
    for j=1:pd
        x(i,j)=l-(l-u)*rand; % Position of the crows in the space
    end
end

function [ft,f,g]=fitness(xn,N,pd,fun,lb,ub) % Function for fitness evaluation

for i=1:N
    xn(i,:)=max(min(xn(i,:),ub),lb);
    [ft(i),f(i),g(:,i)]=feval(fun,xn(i,:)'); % Sphere function
end

%%%%%%%%%%%%%%%%%%%%
