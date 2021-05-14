%_________________________________________________________________________________
%  Salp Swarm Algorithm (SSA) source codes version 1.0
%
%  Developed in MATLAB R2016a
%
%  Author and programmer: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
%
%   Main paper:
%   S. Mirjalili, A.H. Gandomi, S.Z. Mirjalili, S. Saremi, H. Faris, S.M. Mirjalili,
%   Salp Swarm Algorithm: A bio-inspired optimizer for engineering design problems
%   Advances in Engineering Software
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.07.002
%____________________________________________________________________________________

function SSA(fun,foutput,nloop,nsol,nvar,nbit,a,b)
%[FoodFitness,FoodPosition,Convergence_curve]=SSA(N,Max_iter,lb,ub,dim,fobj)
N=nsol;
fobj=fun;
dim=nvar;
Max_iter=nloop;
lb=a;
ub=b;
%%%%% Mod
fpminhist=zeros(1,Max_iter);
fpavghist=zeros(1,Max_iter);
fpmaxhist=zeros(1,Max_iter);
fhist=zeros(1,Max_iter);
ghist=[];
%%%%% End Mod

if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

Convergence_curve = zeros(1,Max_iter);

%Initialize the positions of salps
SalpPositions=initialization(N,dim,ub,lb);


FoodPosition=zeros(1,dim);
FoodFitness=inf;


%calculate the fitness of initial salps

for i=1:size(SalpPositions,1) %1-nsol
    [SalpFitness(1,i),f0(:,i),g0(:,i)]=feval(fobj,SalpPositions(i,:)');
end
fp0=SalpFitness;
[sorted_salps_fitness,sorted_indexes]=sort(SalpFitness);
for newindex=1:N
    Sorted_salps(newindex,:)=SalpPositions(sorted_indexes(newindex),:);
end
FoodPosition=Sorted_salps(1,:);
FoodFitness=sorted_salps_fitness(1);
%%%%% Mod
fp0=fp0(:,sorted_indexes);
g0=g0(:,sorted_indexes);
fpmin0=fp0(:,1);
gmin0=g0(:,1);
%
fpminhist(:,1)=fpmin0;
fpavghist(:,1)=mean(fp0);
fpmaxhist(:,1)=fp0(:,end);
fhist(:,1)=FoodFitness;
ghist=[ghist gmin0];
%%%%% End Mod
%Main loop
l=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
while l<Max_iter+1
    
    c1 = 2*exp(-(4*l/Max_iter)^2); % Eq. (3.2) in the paper
    
    for i=1:size(SalpPositions,1)
        
        SalpPositions= SalpPositions';
        
        if i<=N/2
            for j=1:1:dim
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=FoodPosition(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=FoodPosition(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif i>N/2 && i<N+1
            point1=SalpPositions(:,i-1);
            point2=SalpPositions(:,i);
            
            SalpPositions(:,i)=(point2+point1)/2; % % Eq. (3.4) in the paper
        end
        
        SalpPositions= SalpPositions';
    end
    
    for i=1:size(SalpPositions,1)
        
        Tp=SalpPositions(i,:)>ub';Tm=SalpPositions(i,:)<lb';SalpPositions(i,:)=(SalpPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        
        [SalpFitness(1,i),f(:,i),g(:,i)]=feval(fobj,SalpPositions(i,:)');
        fp=SalpFitness;
        if SalpFitness(1,i)<FoodFitness
            FoodPosition=SalpPositions(i,:);
            FoodFitness=SalpFitness(1,i);
            fpmin0=fp(:,i);
            gmin0=g(:,i);
        end
    end
    
    Convergence_curve(l)=FoodFitness;
    %%%%% Mod
        fpminhist(:,l)=fpmin0;
        fpavghist(:,l)=mean(fp);
        fpmaxhist(:,l)=max(fp);
        fhist(:,l)=FoodFitness;
        ghist=[ghist gmin0];
    %%%%% End Mod
    
    l = l + 1;
end

%%%%% Mod
fp0=fp;
f0=SalpFitness;
x0=SalpPositions';
g0=g;
[fpmin,nmin]=min(fp0);
xmin=x0(:,nmin);
fmin=f0(:,nmin);
gmin=g0(:,nmin);
maxeval=(l-1)*size(fp0,2);
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%% End Mod


