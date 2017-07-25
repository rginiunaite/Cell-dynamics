% Gillespie algorithm for 1d

clear all

% Parameters

% lattice dimension
dim = 1;

% length of the side of square
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.05; % proliferation rate per unit time of giving rise to another agent

Pd = 0; % death rate per unit time


cells = zeros(L,dim); % array that will store the sites, this is just for 1d

% initially choose randomly if each state is occupied or no

% choose a random number for each state
rand_initial = rand(L,dim);

for i = 1:L      
   % if the random number is greater than 0.5, then assume that initially
   % there is a cell at that site. 1 corresponds to that cell
   if rand_initial(i)>0.5
       cells(i) = 1; 
   end
end


Q = sum(cells); % store the number of cells at different times


a0 = (Pm + Pp + Pd)*Q; % total propensity function 


t_final = 20; % time interval
t = 0; % initialise time
iter = 0; % count the number of iterations

%Cells = zeros(L,30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for logistic growth

Q_initial = Q;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

while t < t_final

    iter = iter + 1;
    NumberOfCells(iter) = Q; % keep track of the changes in the number of cells
    %Cells(:,iter) = cells; % check how it looks like
    
    tau = log (1/rand) / a0; % choose time till the next event

    % choose which event occurs (movement, proliferation or death)#

    R = a0*rand;

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
    if R > Pm*Q && R < (Pm+Pp)*Q

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

    if R > (Pm+Pp)*Q  && R < (Pm+Pp + Pd)*Q

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


% plot the changes in the number of cells
figure

store_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison

plot(store_time,NumberOfCells/L,'LineWidth', 2);set(gca,'FontSize',14); % number of cells in the system
%h_legend = legend('1st component','2nd component','Total');
%set(h_legend,'FontSize',14)
xlabel('time','FontSize',14)
ylabel('density','FontSize',14)



%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

cA0 = Q_initial/L; % initial density
cA0_bar = (Pp - Pd)/(Pp)*cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(store_time)
    %store_time(i) = (Pp - Pd) * store_time(i); % rescale time
    y(i) = f(store_time(i)); 
end

%%% plot logistic solution 

hold on
plot(store_time,y,'LineWidth', 2)
h_legend = legend('Gillespie 1D','Logistic growth','Total');
set(h_legend,'FontSize',14)
xlabel('time','FontSize',14)
ylabel('density','FontSize',14)

