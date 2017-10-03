% combination of Gillespie (40 realisations) in 1D
% save data for different parameter values, no scaling

for Number =2
    
        % Parameters

        % lattice dimension
        dim = 1;

        % half length of the side of square
        L = 1000;

        Pm = 1; % transition rate per unit time of moving to another lattice site

        %%%%%%%%%%%%%%%%%%%% Save simulations for three different parameter values

        if Number == 1
                Pp = 0.6; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.4; % death rate per unit time
                initial_percentage = 1;
        end

        if Number == 2
                Pp = 0.4; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.6; % death rate per unit time
                initial_percentage = 1;
        end


    cells = zeros(L,dim); % array that will store the sites, this is just for 1d


    % initially choose randomly if each state is occupied or no

    % choose a random number for each state
    rand_initial = rand(L,dim)/initial_percentage;

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

    t_final = 100; % time interval for Gillespie algorithm

    %%%%%%%
    % Gillespie will be done mutliple times

    Gillespie_iterations = 1;

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



    %rescaled_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison
    rescaled_time = store_time; % no scaling

    TimeReductionFactor = 10; % when death rate is increase, this has to be increased as well, %%%%%%%% even bigger when death rate is bigger

    rescaled_time = rescaled_time(1:round(length(rescaled_time)/TimeReductionFactor)); % make sure that all 40 simulations were running till the end time

    Sum_NumberOfCells = zeros(size(NumberOfCells,1),1);

    for j=1:Gillespie_iterations
        Sum_NumberOfCells = Sum_NumberOfCells + NumberOfCells(:,j);
    end
    average_NumberOfCells = Sum_NumberOfCells/Gillespie_iterations;

    %rescaled_density = Pp/(Pp - Pd) * average_NumberOfCells/L;
    rescaled_density =  average_NumberOfCells/L; % no scaling


    % save data

    save_time = rescaled_time;
    save_density = rescaled_density(1:length(rescaled_time));
    save_cA0 = Q_initial/L; % initial density
    save_store_time = store_time; % these have to be saved

    if Number == 1 && rescaled_time(length(rescaled_time)) ~= 0 % ignore the solutions where time jumps back to zero and starts moving again
        save('gillespie1D_large_initial_Pp06_Pd04_40G_long.mat','save_time','save_density','save_cA0','save_store_time');

    end
    
    if Number == 2 && rescaled_time(length(rescaled_time)) ~= 0
        save('gillespie1D_large_initial_Pp04_Pd06_40G_long.mat','save_time','save_density','save_cA0','save_store_time');

    end    
    

end