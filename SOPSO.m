function SOPSO(fun,foutput,nloop,nsol,nvar,nbit,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single objective Particle Swarm Optimization %
% Created by Sujin Bureerat %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DEFINITION
% a is lower limit
% b is upper limit
% nvar is the number of design variables
% nloop is the number of iteration
% nsol is the number of individuals(genes) in population
% opt = the optimum solution
% fitest = the fitness of the optimum
% avrg= the avrage of the fitness on each generation
% bestfit=the fitest on each generation
%
% lsol is the last population as the real desing parameters
% lpop the population as bit strings in the last genmeration
% lfit the fitnesses of the last generation
%%%%%%%%%%Initialisation%%%%%%%%%%%%
rand('state',sum(100*clock));
Wst = 0.5;%Starting inertia weight
Wen = 0.01;%Ending inertia weight
C1 = 2.8;%cognitive learning factor
C2 = 1.3;%social learning factor
a=0;% lower bound
b=1;% upper bound

[bin0,x0,fp0,f0,g0] = ea_initial(fun,nbit,nvar,nsol);
pop0=x0;
[fpmin,npmin]=min(fp0);
fpminhist=min(fp0);fpavghist=mean(fp0);fpmaxhist=max(fp0);
fhist=f0(npmin);ghist=g0(:,npmin);

[p_leader,fp_leader,f_leader,g_leader]=pso_selection(pop0,fp0,f0,g0,[],[],[],[]);
Vel0=zeros(size(pop0));% initial velocities of the initial population pop0
pbest=pop0;fpbest=fp0;% initial personal best of pop0i
%%%%%%%%% Start MOPSO search %%%%%%%%
tic %
iter = 0;
while iter < nloop
    iter = iter+1;
    Wi=Wst+((Wen-Wst)/nloop)*iter;
    for i=1:nsol
        % selecting p_leader randomly from the external archive
        popi=pop0(:,i);
        pbesti=pbest(:,i);
        % update particls' velocities
        Vel1(:,i)=Wi*Vel0(:,i)+C1*rand*(pbesti-popi)+C2*rand*(p_leader-popi);
        % update particles
        pop1(:,i)=max(0,min(1,pop0(:,i)+Vel1(:,i)));
        %[f1(:,i),g1(:,i)]=feval(fun,bin0(:,i),level,geofile);
        % function evaluation
        [fp1(i),f1(i),g1(:,i)]=feval(fun,pop1(:,i));
        % update personal best
        if fp1(i)< fpbest(i);
            pbest(:,i)=pop1(:,i);
            fpbest(i)=fp1(i);
        end
    end
    % selection
    [p_leader,fp_leader,f_leader,g_leader]=pso_selection(pop1,fp1,f1,g1,...
        p_leader,fp_leader,f_leader,g_leader);
    pop00=pop0;
    pop0=pop1;fp0=fp1;f0=f1;g0=g1;Vel0=Vel1;
    
    %%%%%%%%
    
    fpminhist=[fpminhist fp_leader];fpavghist=[fpavghist mean(fp0)];
    fpmaxhist=[fpmaxhist max(fp0)];
    fhist=[fhist f_leader];ghist=[ghist g_leader];
end

% figure(1),clf,hold on
% plot((1:length(besthist))*nloop*nsol/length(besthist),besthist,'r')
% plot((1:length(besthist))*nloop*nsol/length(besthist),avghist,'b')
% plot((1:length(besthist))*nloop*nsol/length(besthist),worsthist,'g')

fpmin=fp_leader;
fmin=f_leader;
xmin=p_leader;
gmin=g_leader;
Time=toc % second
maxeval=nloop*nsol;
save(foutput,'xmin','fpmin','fmin','gmin','maxeval',...
    'fpminhist','fpavghist','fpmaxhist','fhist','ghist')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-programs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pop0,f0,g0] = soearsminitial(fun,nsol,nvar)
%
% Randomly initiate the population, design variables
%
pop0 = (rand(nvar,nsol));
for i=1:nsol
[f0(:,i),g0(:,i)]=feval(fun,pop0(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest,fpbest,fbest,gbest]=pso_selection(x0,fp0,f0,g0,...
    xbest,fpbest,fbest,gbest);
x=[xbest x0];
fp=[fpbest fp0];
f=[fbest f0];
g=[gbest g0];
[fpmin,nmin]=min(fp);
xbest=x(:,nmin);
fpbest=fp(:,nmin);
gbest=g(:,nmin);
fbest=f(:,nmin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%