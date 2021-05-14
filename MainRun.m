%%% Main run for Comparative Performance of Twelve Metaheuristics for Wind Farm Layout Optimisation
%%% Code from the paper;
%%% Kunakote, T., Sabangban, N., Kumar, S., Tejani, G. G., Panagant, N., Pholdee, N., S.Bureerat & Yildiz, A. R. (2021), 
%%% Comparative Performance of Twelve Metaheuristics for Wind Farm Layout Optimisation, Archives of Computational Methods in Engineering, 1-14.


clc
clear all

% run single objective evolutionary algorithms
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


addpath([pwd '\Mod_Algo\SSA'])%add the folder of SSA
addpath([pwd '\Mod_Algo\GOA'])%add the folder of GOA
nvar=[100 100 39 39];% no. of design variables of the test problems

for i=1:12  % evolutionary algorithms
    algoi=char(algo(i,:));
    for j=1:4  % WFLO problem
        funj=char(fobj(j,:));
        nvarj=nvar(j);
        nsol=50;  %%population size
        nloop=50; %% Number of generation
        nrun=10;  % no. of optimisation runs
        nbit=10;
        a=zeros(nvarj,1);b=ones(nvarj,1);
        parfor k=1:nrun
            [i j k]
            feval(algoi,funj,['rst' numtostr(i) numtostr(j) numtostr(k)],nloop,nsol,nvarj,nbit,a,b);

        end
    end
end




