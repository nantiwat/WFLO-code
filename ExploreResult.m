%%% Explore the restuls
%%% Code from the paper;
%%% Kunakote, T., Sabangban, N., Kumar, S., Tejani, G. G., Panagant, N., Pholdee, N., S.Bureerat & Yildiz, A. R. (2021), 
%%% Comparative Performance of Twelve Metaheuristics for Wind Farm Layout Optimisation, Archives of Computational Methods in Engineering, 1-14.

clc
clear all

fobj=[{'wflo_partialRotor01'}
    {'wflo_fullRotor01'}
    {'wflo_partialRotor02'}
    {'wflo_fullRotor02'}];

algo=[{'SOABC'}  %  Artificial bee colony method
    {'SOACOR'}        %  Real-code ant colony optimisation
    {'SODE'}          %  Differential evolution
    {'SOPSO'}         %  Particle swarm optimisation
    {'SOTLBO'}        %  Teaching-learning based optimisation
    {'SOCMAES'}       %  Evolution strategy with covarience matrix adatation
    {'SOMFO'}         % Moth-flame Optimization Algorithm
    {'SOSCA'}         % The Sine Cosine Algorithm
    {'SOWOA'}         %35Whale Optimization Algorithm-
    {'SOCSA'}         %crow search Optimization Algorithm
    {'SOSSA'}         % Salp Swarm Optimizaer
    {'SOGOA'}         %   Grassopper optimization algorithm
  ];

nrun=10;  % no. of optimisation runs
for j=1:4
    coutAl=0;
    for ii=1:12 % evolutionary algorithms
        coutAl=coutAl+1;
        algoi=char(algo(ii,:));
        xlswrite('CostWindlayout',{algoi},['sheet' num2str(j)],['A' num2str(coutAl+1)])
 
        for k=1:nrun % no. of optimisation runs
            load(['rst' numtostr(ii) numtostr(j) numtostr(k)])
            fout(coutAl,k)=fpmin;
            clear fpmin
        end
    end
    fmean=mean(fout,2);
    fstd=std(fout,[],2);
    fmin=min(fout,[],2);
    fmax=max(fout,[],2);
    xlswrite('CostWindlayout',[fmean,fstd,fmin,fmax],['sheet' num2str(j)],'B2')
    [~,nmin]=min(fmean);
    sal(j)=nmin;
    
end



%%% plot

clear fout
for j=1:4   % Problem no.
    al=sal(j);
    funj=char(fobj(j,:));
    for ii=al  % evolutionary algorithms
        
        algoi=char(algo(ii,:));
 
        for k=1:nrun % no. of optimisation runs
            load(['rst' numtostr(ii) numtostr(j) numtostr(k)])
            fout(k)=fpmin;
            clear fpmin
        end
        [fmin,nmin]=min(fout);
        load(['rst' numtostr(ii) numtostr(j) numtostr(nmin(1))])
        fp=feval(funj,xmin,2);
        
        
   end
end
   

