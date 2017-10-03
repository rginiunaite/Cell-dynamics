% combination of Gillespie on a hexagonal lattice, in "D with hexagonal
% lattice, save data for different parameter valeus


for Number =4
    
        % Parameters

        % lattice dimension
        dim = 2;

        % half length of the side of square
        L = 100;

        Pm = 1; % transition rate per unit time of moving to another lattice site

        %%%%%%%%%%%%%%%%%%%% Save simulations for three different parameter values

        if Number == 1
                Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.1; % death rate per unit time
                initial_percentage = 1;
        end

        if Number == 2
                Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.1; % death rate per unit time
                initial_percentage = 1.9;
        end

        if Number == 3
                Pp = 1; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.5; % death rate per unit time
                initial_percentage = 1.9;
        end

        if Number == 4
                Pp = 0.6; % proliferation rate per unit time of giving rise to another agent

                Pd = 0.4; % death rate per unit time
                initial_percentage = 1.9;
        end


    cells = zeros(L,L); % array that will store the sites

    % initially choose randomly if each state is occupied or no

    % choose a random number for each state
    
    rand_initial = rand(L,L)/initial_percentage; 

    for i = 1:L      
       % if the random number is greater than 0.5, then assume that initially
       % there is a cell at that site. 1 corresponds to that cell
       for j = 1:L
          if rand_initial(i,j)>0.5
            cells(i,j) = 1; 
          end
       end
    end

    cells_initial = cells;

    Q = sum(sum(cells)); % store the number of cells at different times


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

    %store_time = zeros(1,len); 

    for j =1:Gillespie_iterations

        a0 = (Pm + Pp + Pd)*Q_initial; % total propensity function 
        Q = Q_initial; % counts the number of cells
        cells = cells_initial; % keeps track of cells at different times 


        t = 0; % initialise time
        iter = 0; % count the number of iterations


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        while t < t_final

            iter = iter + 1;
            NumberOfCells(iter,j) = Q;
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


                % choose one of its 6 neighbours (a,b)
                random_number = rand; 

                % neighbours are different based on the column ID, if it is odd
                % or even 

                if mod(agent_ID_column,2) == 1 % if it is odd

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5     
                       % neighbour 1, (a,b-1)
                       if random_number <= 0.17
                          neighbour_ID_row = agent_ID_row - 1; 
                          neighbour_ID_column = agent_ID_column; 
                          %periodic boundary conditions
                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end

                       % neighbour 2, (a+1,b)       
                       elseif random_number > 0.17 && random_number <= 0.34
                          neighbour_ID_row = agent_ID_row;
                          neighbour_ID_column = agent_ID_column+1;       
                          %periodic boundary conditions
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end

                       % neighbour 3, (a,b+1)
                       elseif random_number > 0.34 && random_number <= 0.51
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row+1;
                          %periodic boundary conditions
                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end        

                       % neighbour 4, (a-1,b+1)
                       elseif random_number > 0.51 && random_number <= 0.68
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row +1;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    

                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end

                        % neighbour 5, (a-1,b)
                       elseif random_number > 0.68 && random_number <= 0.85
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    


                       % neighbour 6, (a-1,b-1)
                       elseif random_number > 0.85 && random_number <= 1
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row-1;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    

                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end    

                       end                
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                   else % if b is even

                       % neighbour 1, (a+1,b-1)
                       if random_number <= 0.17
                          neighbour_ID_row = agent_ID_row - 1; 
                          neighbour_ID_column = agent_ID_column +1; 
                          %periodic boundary conditions
                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end

                           if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end  


                       % neighbour 2, (a+1,b)       
                       elseif random_number > 0.17 && random_number <= 0.34
                          neighbour_ID_row = agent_ID_row;
                          neighbour_ID_column = agent_ID_column+1;       
                          %periodic boundary conditions
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end

                       % neighbour 3, (a+1,b+1)
                       elseif random_number > 0.34 && random_number <= 0.51
                          neighbour_ID_column = agent_ID_column+1;
                          neighbour_ID_row = agent_ID_row+1;
                          %periodic boundary conditions
                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end        
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end   

                       % neighbour 4, (a,b+1)
                       elseif random_number > 0.51 && random_number <= 0.68
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row +1;
                          %periodic boundary conditions

                           if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                           end


                        % neighbour 5, (a-1,b)
                       elseif random_number > 0.68 && random_number <= 0.85
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end

                        % beighbour 6, (a,b-1)
                       elseif random_number > 0.85 && random_number <= 1
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row-1;
                          %periodic boundary conditions

                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end    

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



                % choose one of its 6 neighbours (a,b)
                random_number = rand; 

                % neighbours are different based on the column ID, if it is odd
                % or even 

                if mod(agent_ID_column,2) == 1 % if it is odd

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5     
                       % neighbour 1, (a,b-1)
                       if random_number <= 0.17
                          neighbour_ID_row = agent_ID_row - 1; 
                          neighbour_ID_column = agent_ID_column; 
                          %periodic boundary conditions
                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end

                       % neighbour 2, (a+1,b)       
                       elseif random_number > 0.17 && random_number <= 0.34
                          neighbour_ID_row = agent_ID_row;
                          neighbour_ID_column = agent_ID_column+1;       
                          %periodic boundary conditions
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end

                       % neighbour 3, (a,b+1)
                       elseif random_number > 0.34 && random_number <= 0.51
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row+1;
                          %periodic boundary conditions
                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end        

                       % neighbour 4, (a-1,b+1)
                       elseif random_number > 0.51 && random_number <= 0.68
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row +1;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    

                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end

                        % neighbour 5, (a-1,b)
                       elseif random_number > 0.68 && random_number <= 0.85
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    


                       % neighbour 6, (a-1,b-1)
                       elseif random_number > 0.85 && random_number <= 1
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row-1;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end    

                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end    

                       end                
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                   else % if b is even

                       % neighbour 1, (a+1,b-1)
                       if random_number <= 0.17
                          neighbour_ID_row = agent_ID_row - 1; 
                          neighbour_ID_column = agent_ID_column +1; 
                          %periodic boundary conditions
                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end

                           if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end  


                       % neighbour 2, (a+1,b)       
                       elseif random_number > 0.17 && random_number <= 0.34
                          neighbour_ID_row = agent_ID_row;
                          neighbour_ID_column = agent_ID_column+1;       
                          %periodic boundary conditions
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end

                       % neighbour 3, (a+1,b+1)
                       elseif random_number > 0.34 && random_number <= 0.51
                          neighbour_ID_column = agent_ID_column+1;
                          neighbour_ID_row = agent_ID_row+1;
                          %periodic boundary conditions
                          if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                          end        
                          if neighbour_ID_column == L+1
                             neighbour_ID_column = 1;  
                          end   

                       % neighbour 4, (a,b+1)
                       elseif random_number > 0.51 && random_number <= 0.68
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row +1;
                          %periodic boundary conditions

                           if neighbour_ID_row == L+1
                             neighbour_ID_row = 1;  
                           end


                        % neighbour 5, (a-1,b)
                       elseif random_number > 0.68 && random_number <= 0.85
                          neighbour_ID_column = agent_ID_column -1;
                          neighbour_ID_row = agent_ID_row;
                          %periodic boundary conditions
                          if neighbour_ID_column == 0
                             neighbour_ID_column = L;  
                          end

                        % beighbour 6, (a,b-1)
                       elseif random_number > 0.85 && random_number <= 1
                          neighbour_ID_column = agent_ID_column;
                          neighbour_ID_row = agent_ID_row-1;
                          %periodic boundary conditions

                          if neighbour_ID_row == 0
                             neighbour_ID_row = L;  
                          end    

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


        % if a0=0, tau infinity, for further usage of timeset maximal time to the
        % previous one
        if store_time(length(store_time))==Inf
            store_time(length(store_time)) = store_time(length(store_time)-1);
        end


    end



    rescaled_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison
    %rescaled_time = store_time; % rescale time to allow parameter comparison, temporarily for big death rate


    TimeReductionFactor = 1.5; % when death rate is increased, this has to be increased as well

    rescaled_time = rescaled_time(1:round(length(rescaled_time)/TimeReductionFactor)); % make sure that all 40 simulations were running till the end time

    Sum_NumberOfCells = zeros(size(NumberOfCells,1),1);

    for j=1:Gillespie_iterations
        Sum_NumberOfCells = Sum_NumberOfCells + NumberOfCells(:,j);
    end
    average_NumberOfCells = Sum_NumberOfCells/Gillespie_iterations;

    %rescaled_density = (Pp - Pd)/Pp * average_NumberOfCells/L^2; %wrong
    rescaled_density = Pp/(Pp - Pd) * average_NumberOfCells/L^2;
    %rescaled_density =  average_NumberOfCells/L^2; % for big death rate



    % plot(rescaled_time,rescaled_density(1:length(rescaled_time)),'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
    % %h_legend = legend('1st component','2nd component','Total');
    % %set(h_legend,'FontSize',14)
    % xlabel('rescaled time','FontSize',14)
    % ylabel('rescaled density','FontSize',14)

    % save data

    save_time = rescaled_time;
    save_density = rescaled_density(1:length(rescaled_time));
    save_cA0 = Q_initial/L^2; % initial density
    save_store_time = store_time; % these have to be saved



    if Number == 1 && rescaled_time(length(rescaled_time)) ~= 0 % ignore the solutions where time jumps back to zero and starts moving again
        save('gillespie2D_large_initial_Pp05_Pd01_40G_hexa.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_large_initial_Pp05_Pd01_40G_hexa.mat');
    end
    
    if Number == 2 && rescaled_time(length(rescaled_time)) ~= 0
        save('gillespie2D_Pp05_Pd01_40G_hexa.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_Pp05_Pd01_40G_hexa.mat');
    end    
    
    if Number == 3 && rescaled_time(length(rescaled_time)) ~= 0
        save('gillespie2D_Pp1_Pd05_40G_hexa.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_Pp1_Pd05_40G_hexa.mat');
    end
    
    if Number == 4 && rescaled_time(length(rescaled_time)) ~= 0
        save('gillespie2D_Pp06_Pd04_40G_hexa_rescaled.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_Pp06_Pd04_40G_hexa_rescaled.mat');
    end

end


