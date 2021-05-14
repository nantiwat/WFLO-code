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
N=nsol;
Max_iter=nloop;
lb=a';
ub=b';
dim=nvar;
fobj=fun;

if size(ub,1)==1
    ub=ones(dim,1)*ub;
    lb=ones(dim,1)*lb;
end


%Initialize the population of grasshoppers
GrassHopperPositions=initialization(N,dim,ub,lb);
GrassHopperFitness = zeros(1,N);

fitness_history=zeros(N,Max_iter);
position_history=zeros(N,Max_iter,dim);
Convergence_curve=zeros(1,Max_iter);
Trajectories=zeros(N,Max_iter);

cMax=1;
cMin=0.00004;
%Calculate the fitness of initial grasshoppers

for i=1:size(GrassHopperPositions,1)
    [fp0(i),f0(i),g0(:,i)]=feval(fobj,GrassHopperPositions(i,:)');
    GrassHopperFitness(1,i)=fp0(i);
    fitness_history(i,1)=GrassHopperFitness(1,i);
    position_history(i,1,:)=GrassHopperPositions(i,:);
    Trajectories(:,1)=GrassHopperPositions(:,1);
end

[fpmin,npmin]=min(fp0);
fpminhist=min(fp0);fpavghist=mean(fp0);fpmaxhist=max(fp0);
fhist=f0(npmin);ghist=g0(:,npmin);

[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);

% Find the best grasshopper (target) in the first population 
for newindex=1:N
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end

TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);

% Main loop
l=2; % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
while l<Max_iter+1
    
    c=cMax-l*((cMax-cMin)/Max_iter); % Eq. (2.8) in the paper
    
     for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
%         for k=1:dim-1
%             S_i=zeros(2,1);
            S_i_total=zeros(dim,1);
            for j=1:N
                if i~=j
                    Dist=sqrt((temp(:,j)-temp(:,i))'*(temp(:,j)-temp(:,i))); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(:,j)-temp(:,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((b - a)*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i_total=S_i_total+s_ij;
                end
            end
%             S_i_total(k:k+1,1) = S_i;
            
%         end
        X_new = c * S_i_total'+ (TargetPosition); % Eq. (2.7) in the paper  
        GrassHopperPositions_temp(i,:)=X_new'; 
    end
    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp;
    
    for i=1:size(GrassHopperPositions,1)
        % Relocate grasshoppers that go outside the search space 
        for j=1:dim
            if GrassHopperPositions(i,j)>b(j)
                GrassHopperPositions(i,j)=b(j);
            elseif GrassHopperPositions(i,j) < a(j)
                GrassHopperPositions(i,j)=a(j);
            end
        end
        
        % Calculating the objective values for all grasshoppers
        [fp0(i),f0(i),g0(:,i)]=feval(fobj,GrassHopperPositions(i,:)');
        GrassHopperFitness(1,i)=fp0(i);
        fitness_history(i,l)=GrassHopperFitness(1,i);
        position_history(i,l,:)=GrassHopperPositions(i,:);
        
        Trajectories(:,l)=GrassHopperPositions(:,1);
        
        % Update the target
        if GrassHopperFitness(1,i)<TargetFitness
            TargetPosition=GrassHopperPositions(i,:);
            TargetFitness=GrassHopperFitness(1,i);
        end
    end
        
    [fpmin,npmin]=min(fp0);
    fpminhist=[fpminhist fpmin];fpavghist=[fpavghist mean(fp0)];
    fpmaxhist=[fpmaxhist max(fp0)];
    fhist=[fhist f0(npmin)];ghist=[ghist g0(:,npmin)];
    
    Convergence_curve(l)=TargetFitness;    
    l = l + 1;
end

if (flag==1)
    TargetPosition = TargetPosition(1:dim-1);
end

xmin=TargetPosition';
[fpmin,fmin,gmin]=feval(fun,xmin);
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X]=initialization(N,dim,up,down)

if size(up,1)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
function d = distance(a,b)
d=sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
function o=S_func(r)
f=0.5;
l=1.5;
o=f*exp(-r/l)-exp(-r);  % Eq. (2.3) in the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

