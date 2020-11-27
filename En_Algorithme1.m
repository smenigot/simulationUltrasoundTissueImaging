% Optimization Algorithm
%
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr
function resultats = En_Algorithme(optim, dir)

% Initialization
Nitermax = optim.Nitermax;
medium.flag = false;    % Medium initialization

%% I Setup the GA
npar=optim.Npoint;                % number of optimization variables

varhi=-5;               % variable limits
varlo=5;     
%% II Stopping criteria
maxit=Nitermax;         % max number of iterations
mincost=-9999999;       % minimum cost
%% III GA parameters
popsize=12;                         % set population size
mutrate=.4;                         % set mutation rate
selection=0.5;                      % fraction of population kept
Nt=npar;                            % continuous parameter GA Nt=#variables
keep=floor(selection*popsize);      % population members that survive
nmut=ceil((popsize-1)*Nt*mutrate);  % total number of mutations
M=ceil((popsize-keep)/2);           % number of matings
%% Create the initial population
fprintf('Initialisation\n');
iga=0;                              % generation counter initialized

par=(varhi-varlo)*rand(popsize,npar)+varlo;  % random

cost = zeros(popsize,1);
for k2 = 1:popsize
    parameter(k2).A          = optim.Amp;	% [Pa]  Pulse amplitude
    parameter(k2).Niter      = k2;
    parameter(k2).Nc         = optim.Nc;
    parameter(k2).bruit      = par(k2,:);
    if k2 == 1
        [Energie,medium] = En_FindEnergy(parameter(k2), medium);
        cost(k2) = -Energie.optim;
        medium.flag = true;
    end
end

for k2 = 2:popsize
    [Energie,medium_temp] = En_FindEnergy(parameter(k2), medium);
    cost(k2) = -Energie.optim;
end

[cost,ind]=sort(cost);              % min cost in element 1
par=par(ind,:);                     % sort continuous
minc(1)=min(cost);                  % minc contains min of
meanc(1)=mean(cost);                % meanc contains mean of population


resultats.A         = optim.Amp;
resultats.Niter     = 1;
resultats.Eiter     = minc(1);
resultats.EiterMean = meanc(1);
resultats.bruit     = par(1,:);

save([dir '/resFundamentalFT-f'],'resultats');

fprintf('\n');

%% Iterate through generations
while iga<maxit
    iga=iga+1;                              % increments generation counter
    fprintf('Generation n°%i\n',iga);
    %% Pair and mate
    M=ceil((popsize-keep)/2);               % number of matings
    prob=flipud([1:keep]'/sum([1:keep]));   % weights chromosomes
    odds=[0 cumsum(prob(1:keep))'];         % probability distribution function
    pick1=rand(1,M);                        % mate #1
    pick2=rand(1,M);                        % mate #2
    % ma and pa contain the indicies of the chromosomes that will mate
    ic=1;
    while ic<=M
        for id=2:keep+1
            if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
                ma(ic)=id-1;
            end
            if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
                pa(ic)=id-1;
            end
        end
        ic=ic+1;
    end
    %% Performs mating using single point crossover
    ix=1:2:keep;                            % index of mate #1
    xp=ceil(rand(1,M)*Nt);                  % crossover point
    r=rand(1,M);                            % mixing parameter
    for ic=1:M
        xy=par(ma(ic),xp(ic))-par(pa(ic),xp(ic)); % ma and pa mate
        par(keep+ix(ic),:)=par(ma(ic),:);         % 1st offspring
        par(keep+ix(ic)+1,:)=par(pa(ic),:);       % 2nd offspring
        par(keep+ix(ic),xp(ic))=par(ma(ic),xp(ic))-r(ic).*xy;   % 1st
        par(keep+ix(ic)+1,xp(ic))=par(pa(ic),xp(ic))+r(ic).*xy; % 2nd
        if xp(ic)<npar % crossover when last variable not selected
            par(keep+ix(ic),:)=[par(keep+ix(ic),1:xp(ic)) par(keep+ix(ic)+1,xp(ic)+1:npar)];
            par(keep+ix(ic)+1,:)=[par(keep+ix(ic)+1,1:xp(ic)) par(keep+ix(ic),xp(ic)+1:npar)];
        end
    end
    %%     Mutate the population
    mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);
    mcol=ceil(rand(1,nmut)*Nt);
    for ii=1:nmut
        par(mrow(ii),mcol(ii))=(varhi-varlo)*rand+varlo; % mutation
    end
    %%    The new offspring and mutated chromosomes are  evaluated
    for k2 = 2:popsize
        parameter(k2).A          = optim.Amp;	% [Pa]  Pulse amplitude
        parameter(k2).Nc         = optim.Nc;
        parameter(k2).bruit      = par(k2,:);
        cost(k2)                 = 0;
    end
    
    for k2 = 2:popsize
        [Energie,medium_temp] = En_FindEnergy(parameter(k2), medium);
        cost(k2) = -Energie.optim;
    end

    %%     Sort the costs and associated parameters
    [cost,ind]=sort(cost);
    par=par(ind,:);
    parc(iga+1,:) = par(1,:);
    
    %%    Do statistics for a single nonaveraging run
    minc(iga+1)=min(cost);
    meanc(iga+1)=mean(cost);
    fprintf('Contrast = %2.2f\n\n',minc(iga+1));
    
    %%
    resultats.A(iga+1)         = optim.Amp;
    resultats.Niter(iga+1)     = iga+1;
    resultats.Eiter(iga+1)     = minc(iga+1);
    resultats.EiterMean(iga+1) = meanc(iga+1);
    resultats.bruit(iga+1,:)   = parc(iga+1,:);
    
    save([dir '/resFundamentalFT-f'],'resultats');
end

return