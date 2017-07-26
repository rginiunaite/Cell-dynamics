
% Gillespie algorithm for 1d

clear all

% Parameters

% lattice dimension
dim = 2;

% length of the side of square
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.1; % proliferation rate per unit time of giving rise to another agent

Pd = 0.05; % death rate per unit time


cells = zeros(L,L); % array that will store the sites

% initially choose randomly if each state is occupied or no

% choose a random number for each state
rand_initial = rand(L,L);

for i = 1:L      
   % if the random number is greater than 0.5, then assume that initially
   % there is a cell at that site. 1 corresponds to that cell
   for j = 1:L
      if rand_initial(i,j)>0.5
        cells(i,j) = 1; 
      end
   end
end



Q = sum(sum(cells)); % store the number of cells at different times


a0 = (Pm + Pp + Pd)*Q; % total propensity function 


t_final = 20; % time interval
t = 0; % initialise time
iter = 0; % count the number of iterations

%Cells = zeros(L,30);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for logistic growth

Q_initial = Q;

%%%%%%%
%Gillespie simulation

while t < t_final

    iter = iter + 1;
    NumberOfCells(iter) = Q;
    %Cells(:,iter) = cells; % check how it looks like
    
    tau = log (1/rand) / a0; % choose time till the next event

    % choose which event occurs (movement, proliferation or death)#

    R = a0*rand;

    % Movement event
    if R < Pm*Q
        
       random_permutation_row = randperm(L);
       random_permutation_column = randperm(L);
       agent_ID_row = random_permutation_row(1); % choose a random agent 
       agent_ID_column = random_permutation_column(1); % choose a random agent 
       
       % check if it is actually an agent and not zero at that point and if
       % not choose until it is an agent
       while cells(agent_ID_row,agent_ID_column) == 0
           random_permutation_row = randperm(L);
           random_permutation_column = randperm(L);
           agent_ID_row = random_permutation_row(1); % choose a random agent 
           agent_ID_column = random_permutation_column(1); % choose a random agent  
       end
       
       

       % choose one of its 4 neighbours
       random_number = rand;
       % neighbour up
       if random_number < 0.25
          neighbour_ID_row = agent_ID_row + 1; 
          neighbour_ID_column = agent_ID_column; 
          %periodic boundary conditions
          if neighbour_ID_row == L+1
             neighbour_ID_row = 1;  
          end
          
       % neighbour down        
       elseif random_number > 0.25 && random_number < 0.5
          neighbour_ID_row = agent_ID_row -1;
          neighbour_ID_column = agent_ID_column;       
          %periodic boundary conditions
          if neighbour_ID_row == 0
             neighbour_ID_row = L;  
          end
          
       % neighbour left
       elseif random_number > 0.5 && random_number < 0.75
          neighbour_ID_column = agent_ID_column -1;
          neighbour_ID_row = agent_ID_row;
          %periodic boundary conditions
          if neighbour_ID_column == 0
             neighbour_ID_column = L;  
          end        
          
       % neighbour right
       elseif random_number > 0.75
          neighbour_ID_column = agent_ID_column +1;
          neighbour_ID_row = agent_ID_row;
          %periodic boundary conditions
          if neighbour_ID_column == L+1
             neighbour_ID_column = 1;  
          end    
          
       end

       % check if the neighbour is empty, if yes, move the cell to that grid
       if cells(neighbour_ID_row,neighbour_ID_column) == 0
          cells(neighbour_ID_row,neighbour_ID_column) = 1;
          cells(agent_ID_row,agent_ID_column) = 0;

       end

    end

    % Proliferation event
    if R > Pm*Q && R < (Pm+Pp)*Q

       random_permutation_row = randperm(L);
       random_permutation_column = randperm(L);
       agent_ID_row = random_permutation_row(1); % choose a random agent 
       agent_ID_column = random_permutation_column(1); % choose a random agent 
       
       % check if it is actually an agent and not zero at that point and if
       % not choose until it is an agent
       while cells(agent_ID_row,agent_ID_column) == 0
           random_permutation_row = randperm(L);
           random_permutation_column = randperm(L);
           agent_ID_row = random_permutation_row(1); % choose a random agent 
           agent_ID_column = random_permutation_column(1); % choose a random agent  
       end
       
       

       % choose one of its 4 neighbours
       random_number = rand;
       % neighbour up
       if random_number < 0.25
          neighbour_ID_row = agent_ID_row + 1; 
          neighbour_ID_column = agent_ID_column; 
          %periodic boundary conditions
          if neighbour_ID_row == L+1
             neighbour_ID_row = 1;  
          end
          
       % neighbour down        
       elseif random_number > 0.25 && random_number < 0.5
          neighbour_ID_row = agent_ID_row -1;
          neighbour_ID_column = agent_ID_column;       
          %periodic boundary conditions
          if neighbour_ID_row == 0
             neighbour_ID_row = L;  
          end
          
       % neighbour left
       elseif random_number > 0.5 && random_number < 0.75
          neighbour_ID_column = agent_ID_column -1;
          neighbour_ID_row = agent_ID_row;
          %periodic boundary conditions
          if neighbour_ID_column == 0
             neighbour_ID_column = L;  
          end        
          
       % neighbour right
       elseif random_number > 0.75
          neighbour_ID_column = agent_ID_column +1;
          neighbour_ID_row = agent_ID_row;
          %periodic boundary conditions
          if neighbour_ID_column == L+1
             neighbour_ID_column = 1;  
          end    
          
       end

       % check if the neighbour is empty, if yes, move the cell to that grid
       if cells(neighbour_ID_row,neighbour_ID_column) == 0
          cells(neighbour_ID_row,neighbour_ID_column) = 1;
          Q = Q + 1; % increase in total population
          a0 = a0 + (Pm + Pp + Pd); % change in total propensity function
       end

    end
    
    % Death event

    if R > (Pm+Pp)*Q  && R < (Pm+Pp + Pd)*Q

       random_permutation_row = randperm(L);
       random_permutation_column = randperm(L);
       agent_ID_row = random_permutation_row(1); % choose a random agent 
       agent_ID_column = random_permutation_column(1); % choose a random agent 
       
       % check if it is actually an agent and not zero at that point and if
       % not choose until it is an agent
       while cells(agent_ID_row,agent_ID_column) == 0
           random_permutation_row = randperm(L);
           random_permutation_column = randperm(L);
           agent_ID_row = random_permutation_row(1); % choose a random agent 
           agent_ID_column = random_permutation_column(1); % choose a random agent  
       end
       
       % cell dies

       cells(agent_ID_row, agent_ID_column) = 0;
       Q = Q - 1; % decrease in total population
       a0 = a0 - (Pm + Pp + Pd); % change in total propensity function
    end

    t = t + tau;

    store_time(iter) = t;
end

rescaled_time  = (Pp - Pd)*store_time; % rescale time to allow parameter comparison

rescaled_density = ((Pp-Pd)/(Pp)) * NumberOfCells/L^2;


% plot the changes in the number of cells
figure

%plot(rescaled_time,NumberOfCells/L^2,'LineWidth', 2);set(gca,'FontSize',14); % number of cells in the system
plot(rescaled_time,rescaled_density,'LineWidth', 2);set(gca,'FontSize',14); % number of cells in the system


xlabel('time','FontSize',14)
ylabel('density','FontSize',14)



%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

cA0 = Q_initial/L^2; % initial density
cA0_bar = (Pp - Pd)/(Pp)*cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(store_time)
   % store_time(i) = (Pp - Pd) * store_time(i); % rescale time
    y(i) = f(rescaled_time(i)); 
end

%%% plot logistic solution 



hold on
plot(rescaled_time,y,'LineWidth', 2)
h_legend = legend('Gillespie 2D','Logistic growth');
set(h_legend,'FontSize',14)
xlabel('time','FontSize',14)
ylabel('density','FontSize',14)

