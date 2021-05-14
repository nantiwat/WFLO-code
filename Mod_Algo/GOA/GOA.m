%_________________________________________________________________________%
%  Grasshopper Optimization Algorithm (GOA) source codes demo V1.0        %
%                                                                         %
%  Developed in MATLAB R2016a                                             %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%  Main paper: S. Saremi, S. Mirjalili, A. Lewis                          %
%              Grasshopper Optimisation Algorithm: Theory and Application %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004   %
%                                                                         %
%_________________________________________________________________________%

% The Grasshopper Optimization Algorithm
function GOA(fun,foutput,nloop,nsol,nvar,nbit,a,b)
%[TargetFitness,TargetPosition,Convergence_curve,Trajectories,fitness_history, position_history]=
%GOA(N, Max_iter, lb,ub, dim, fobj)

% disp('GOA is now estimating the global optimum for your problem....')

fobj=fun;
Max_iter=nloop;
N=nsol;
dim=nvar;
nbit;
lb=0;
ub=1;

flag=0;
if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end

if (rem(dim,2)~=0) % this algorithm should be run with a even number of variables. This line is to handle odd number of variables
    dim = dim+1;
    ub = [ub; 100];
    lb = [lb; -100];
    flag=1;
end

%Initialize the population of grasshoppers
GrassHopperPositions=initialization(N,dim,ub,lb);
GrassHopperFitness = zeros(1,N);

fitness_history=zeros(N,Max_iter);
position_history=zeros(N,Max_iter,dim);
Convergence_curve=zeros(1,Max_iter);
Trajectories=zeros(N,Max_iter);
%%%%% Mod
fpminhist=zeros(1,Max_iter);
fpavghist=zeros(1,Max_iter);
fpmaxhist=zeros(1,Max_iter);
fhist=zeros(1,Max_iter);
ghist=[];
%%%%% End Mod
cMax=1;
cMin=0.00004;
%Calculate the fitness of initial grasshoppers

for i=1:size(GrassHopperPositions,1)
    %%%%% End Mod
    if flag == 1
        [GrassHopperFitness(1,i),f(:,i),g(:,i)]=feval(fobj,GrassHopperPositions(i,1:end-1)');
    else
        [GrassHopperFitness(1,i),f(:,i),g(:,i)]=feval(fobj,GrassHopperPositions(i,:)');
    end
    fp(i)=GrassHopperFitness(1,i);
    %%%%% End Mod
    fitness_history(i,1)=GrassHopperFitness(1,i);
    position_history(i,1,:)=GrassHopperPositions(i,:);
    Trajectories(:,1)=GrassHopperPositions(:,1);
end

[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);

% Find the best grasshopper (target) in the first population
Sorted_grasshopper=zeros(size(GrassHopperPositions));
for newindex=1:N
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end
% Sorted_grasshopper=GrassHopperPositions(sorted_indexes,:);
TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);
%%%%% Mod
fp=fp(:,sorted_indexes);
g=g(:,sorted_indexes);
fpmin0=fp(:,1);
gmin0=g(:,1);
%
fpminhist(:,1)=fpmin0;
fpavghist(:,1)=mean(fp);
fpmaxhist(:,1)=fp(:,end);
fhist(:,1)=TargetFitness;
ghist=[ghist gmin0];
%%%%% End Mod
% Main loop
l=2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
while l<Max_iter+1 %% Why ?????
    
    c=cMax-l*((cMax-cMin)/Max_iter); % Eq. (2.8) in the paper
    
    for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
        for k=1:2:dim
            S_i=zeros(2,1);
            for j=1:N
                if i~=j
                    Dist=distance(temp(k:k+1,j), temp(k:k+1,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(k:k+1,j)-temp(k:k+1,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((ub(k:k+1) - lb(k:k+1))*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total(k:k+1, :) = S_i;
            
        end
        
        X_new = c * S_i_total'+ (TargetPosition); % Eq. (2.7) in the paper
        GrassHopperPositions_temp(i,:)=X_new'; 
    end
    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp;
    
    for i=1:size(GrassHopperPositions,1)
        % Relocate grasshoppers that go outside the search space 
        Tp=GrassHopperPositions(i,:)>ub';Tm=GrassHopperPositions(i,:)<lb';GrassHopperPositions(i,:)=(GrassHopperPositions(i,:).*(~(Tp+Tm)))+ub'.*Tp+lb'.*Tm;
        
        % Calculating the objective values for all grasshoppers
        %%%%% End Mod
        if flag == 1
            [GrassHopperFitness(1,i),f(:,i),g(:,i)]=feval(fobj,GrassHopperPositions(i,1:end-1)');
        else
            [GrassHopperFitness(1,i),f(:,i),g(:,i)]=feval(fobj,GrassHopperPositions(i,:)');
        end
        fp(i)=GrassHopperFitness(1,i);
        %%%%% End Mod
        fitness_history(i,l)=GrassHopperFitness(1,i);
        position_history(i,l,:)=GrassHopperPositions(i,:);
        
        Trajectories(:,l)=GrassHopperPositions(:,1);
        
        % Update the target
        if GrassHopperFitness(1,i)<TargetFitness
            TargetPosition=GrassHopperPositions(i,:);
            TargetFitness=GrassHopperFitness(1,i);
            fpmin0=fp(:,i);
            gmin0=g(:,i);
        end
    end
        
    Convergence_curve(l)=TargetFitness;
%     disp(['In iteration #', num2str(l), ' , target''s objective = ', num2str(TargetFitness)])
%%%%% Mod
    fpminhist(:,l)=fpmin0;
    fpavghist(:,l)=mean(fp);
    fpmaxhist(:,l)=max(fp);
    fhist(:,l)=TargetFitness;
    ghist=[ghist gmin0];
%%%%% End Mod
    l = l + 1;
    
end


if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end
%%%%% Mod
fp0=fp;
f0=GrassHopperFitness;
x0=GrassHopperPositions';
g0=g;
[fpmin,nmin]=min(fp0);
xmin=x0(:,nmin);
fmin=f0(:,nmin);
gmin=g0(:,nmin);
maxeval=(l-1)*size(fp0,2);
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%% End Mod














