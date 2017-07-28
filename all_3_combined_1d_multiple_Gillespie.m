% combination of Gillespie (40 realisations), logistic and Runge-Kutta in 1D

clear all

% Parameters

% lattice dimension
dim = 1;

% length of the domain
L = 1000;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 1; % proliferation rate per unit time of giving rise to another agent

Pd = 0; % death rate per unit time


cells = zeros(L,dim); % array that will store the sites, this is just for 1d


% initially choose randomly if each state is occupied or no

% choose a random number for each state
rand_initial = rand(L,dim)/1.9;

for i = 1:L      
   % if the random number is greater than 0.5, then assume that initially
   % there is a cell at that site. 1 corresponds to that cell
   if rand_initial(i)>0.5
       cells(i) = 1; 
   end
end

cells_initial = cells;

Q = sum(cells); % store the number of cells at different times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for logistic growth

Q_initial = Q;

t_final = 50; % time interval for Gillespie algorithm

%%%%%%%
% Gillespie will be done mutliple times
    
Gillespie_iterations = 40;

%%%%%
% Preallocate size approximately to save time, not sure how much time it saves 
a0 = (Pm + Pp + Pd)*Q_initial;
tau = log (1/rand) / a0; % first time step

len = round(t_final/tau); % assume that there will be len simulations

NumberOfCells = zeros (len,Gillespie_iterations); % approximately len iterations for each Gillespie iteration

store_time = zeros(1,len); 

for j =1:Gillespie_iterations

    a0 = (Pm + Pp + Pd)*Q_initial; % total propensity function 
    Q = Q_initial;
    cells = cells_initial;
        
    
    t = 0; % initialise time
    iter = 0; % count the number of iterations

    %Cells = zeros(L,30);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while t < t_final

        iter = iter + 1;
        NumberOfCells(iter,j) = Q; % keep track of the changes in the number of cells
        %Cells(:,iter) = cells; % check how it looks like


        tau = log (1/rand) / a0; % choose time till the next event

        % choose which event occurs (movement, proliferation or death)#

        R = a0*rand; % a0 is updated when the event affecting it occurs

        % Movement event
        if R < Pm*Q

           random_permutation = randperm(L);
           agent_ID = random_permutation(1); % choose a random agent 

           % check if it is actually an agent and not zero at that point and if
           % not choose until it is an agent
           while cells(agent_ID) == 0
            random_permutation = randperm(L);
            agent_ID = random_permutation(1); % choose a random agent   
           end



           % choose one of its neighbours
           if rand > 0.5
              neighbour_ID = agent_ID + 1; 

              %periodic boundary conditions
              if neighbour_ID == L+1
                 neighbour_ID = 1;  
              end

           else
              neighbour_ID = agent_ID -1;

              %periodic boundary conditions
              if neighbour_ID == 0
                 neighbour_ID = L;  
              end        
           end

           % check if the neighbour is empty, if yes, move the cell to that grid
           if cells(neighbour_ID) == 0
              cells(neighbour_ID) = 1;
              cells(agent_ID) = 0;

           end

        end

        % Proliferation event
        if R >= Pm*Q && R < (Pm+Pp)*Q

           random_permutation = randperm(L);
           agent_ID = random_permutation(1); % choose a random agent 

           % check if it is actually an agent and not zero at that point and if
           % not choose until it is an agent
           while cells(agent_ID) == 0
            random_permutation = randperm(L);
            agent_ID = random_permutation(1); % choose a random agent   
           end

           % choose one of its neighbours
           if rand > 0.5
              neighbour_ID = agent_ID + 1; 

              %periodic boundary conditions
              if neighbour_ID == L+1
                 neighbour_ID = 1;  
              end

           else
              neighbour_ID = agent_ID -1;

              %periodic boundary conditions
              if neighbour_ID == 0
                 neighbour_ID = L;  
              end  
           end

           % check if the neighbour is empty, if yes, agent proliferate by placing
           % a new agent in the target site
           if cells(neighbour_ID) == 0
              cells(neighbour_ID) = 1;
              Q = Q + 1; % increase in total population
              a0 = a0 + (Pm + Pp + Pd); % change in total propensity function

           end

        end   

        % Death event

        if R >= (Pm+Pp)*Q  && R < (Pm+Pp + Pd)*Q

           random_permutation = randperm(L);
           agent_ID = random_permutation(1); % choose a random agent 

           % check if it is actually an agent and not zero at that point and if
           % not choose until it is an agent
           while cells(agent_ID) == 0
            random_permutation = randperm(L);
            agent_ID = random_permutation(1); % choose a random agent   
           end


           % cell dies

           cells(agent_ID) = 0;
           Q = Q - 1; % decrease in total population
           a0 = a0 - (Pm + Pp + Pd); % change in total propensity function
        end

        t = t + tau;

        store_time(iter) = t;
    end


    % if a0=0, tau infinity, for further usage of timeset maximal time to the
    % previous one
    if store_time(length(store_time))==Inf
        store_time(length(store_time)) = store_time(length(store_time)-1);
    end
    

end



rescaled_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison

TimeReductionFactor = 2; % when death rate is increase, this has to be increased as well

rescaled_time = rescaled_time(1:round(length(rescaled_time)/TimeReductionFactor)); % make sure that all 40 simulations were running till the end time

Sum_NumberOfCells = zeros(size(NumberOfCells,1),1);

for j=1:Gillespie_iterations
    Sum_NumberOfCells = Sum_NumberOfCells + NumberOfCells(:,j);
end
average_NumberOfCells = Sum_NumberOfCells/Gillespie_iterations;

rescaled_density = (Pp - Pd)/Pp * average_NumberOfCells/L;

% plot the changes in the number of cells
figure

plot(rescaled_time,rescaled_density(1:length(rescaled_time)),'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
%h_legend = legend('1st component','2nd component','Total');
%set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)



%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

cA0 = Q_initial/L; % initial density
cA0_bar = (Pp - Pd)/(Pp)*cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(rescaled_time)
       y(i) = f(rescaled_time(i)); 
end

%%% plot logistic solution 

hold on
plot(rescaled_time,y,'LineWidth', 3)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%
% adding KSA


cA0 = Q_initial/L; % initial density of cells

dt = 0.01; % time step

N = round(store_time(length(1:round(length(store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie

% Initial conditions 

y  = zeros (L+1,N);

y(:,1) = cA0^2 * ones(L+1,1); % initial density of two-point distribution functions 
y(1,1) = cA0; % initial condition for 1-point distribution function


t = 0; % start time


for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd);
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd);
    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    
    store_time_Runge(i) = t;
end


rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = (Pp-Pd)/Pp * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
h_legend = legend('Gillespie 1D','Logistic growth','KSA');
set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['L = ' num2str(L) ', P_m = ' num2str(Pm) ', P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)


function deriv = dynamics(t,y,L,Pm,Pp,Pd)
    

    deriv = zeros(L+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*(y(3) - y(2)) - 2* Pd * y(2) + ...
        Pp * (y(1) - y(2)) + Pp * (y(3)*(y(3)*(y(1)-y(2))^2)/(y(1)^2*(1-y(1))))); % dynamics of two-point distribution function (distance 1)
    
 % dynamics of all two-point functions from distance 2 to L-1
    for i = 3:L
        deriv(i) = (Pm*(y(i-1) + y(i+1) - 2* y(i)) - 2 * Pd * y(i) +...
            Pp * ((y(1)-y(i))*(y(1)-y(2))*(y(i-1)+y(i+1)))/(y(1)^2*(1-y(1))));
    end

    deriv(L+1) = (Pm*(y(L) + y(2) - 2* y(L+1)) - 2 * Pd * y(L+1) +...
            Pp * ((y(1)-y(L+1))*(y(1)-y(2))*(y(L)+y(2)))/(y(1)^2*(1-y(1)))); % for the last one, since we have
        % periodic boundary condition L+1 neighbour is 1 neighbour, i.e. y(2)
    
end

