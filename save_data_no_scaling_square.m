% combination of Gillespie (40 realisations) i n"D square lattice
% no scaling, allowing Pd > Pp, save data
clear all

for Number =1:2
    
    % Parameters

    % lattice dimension
    dim = 2;

    % half length of the side of square
    L = 100;

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
    

    cells = zeros(L,L); % array that will store the sites

    % initially choose randomly if each state is occupied or no

    % choose a random number for each state
    rand_initial = rand(L,L);%/1.9;

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

    t_final = 70; % time interval for Gillespie algorithm

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
        Q = Q_initial; % set this for each new iteration
        cells = cells_initial;


        t = 0; % initialise time
        iter = 0; % count the number of iterations

        %Cells = zeros(L,30);



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
        b(iter)=tau;
        a0_iter(iter) = a0;

        store_time(iter) = t;
    end


        % if a0=0, tau infinity, for further usage of timeset maximal time to the
        % previous one
        if store_time(length(store_time))==Inf
            store_time(length(store_time)) = store_time(length(store_time)-1);
        end


    end



    %rescaled_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison


    TimeReductionFactor = 1.5; % when death rate is increase, this has to be increased as well

    store_time = store_time(1:round(length(store_time)/TimeReductionFactor)); % make sure that all 40 simulations were running till the end time

    Sum_NumberOfCells = zeros(size(NumberOfCells,1),1);

    for j=1:Gillespie_iterations
        Sum_NumberOfCells = Sum_NumberOfCells + NumberOfCells(:,j);
    end
    average_NumberOfCells = Sum_NumberOfCells/Gillespie_iterations;

    %rescaled_density = (Pp - Pd)/Pp * average_NumberOfCells/L^2; %wrong
    density = average_NumberOfCells/L^2;


    % plot the changes in the number of cells
    %figure

    % plot(rescaled_time,rescaled_density(1:length(rescaled_time)),'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
    % %h_legend = legend('1st component','2nd component','Total');
    % %set(h_legend,'FontSize',14)
    % xlabel('rescaled time','FontSize',14)
    % ylabel('rescaled density','FontSize',14)

    % save data

    save_time = store_time;
    save_density = density(1:length(store_time));
    save_cA0 = Q_initial/L^2; % initial density
    save_store_time = store_time; % these have to be saved, same as save_time in this case



    if Number == 1 && store_time(length(store_time)) ~= 0 % ignore the solutions where time jumps back to zero and starts moving again
        save('gillespie2D_large_initial_Pp06_Pd04_square_40G.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_large_initial_Pp06_Pd04_square_40G.mat');
    end
    
    if Number == 2 && store_time(length(store_time)) ~= 0
        save('gillespie2D_large_initial_Pp04_Pd06_square_40G.mat','save_time','save_density','save_cA0','save_store_time');
        S = load('gillespie2D_large_initial_Pp04_Pd06_square_40G.mat');
    end    

end


function deriv = dynamics(t,y,Pm,Pp,Pd,max_rad,dr)
    
    

    numberDer = (round(max_rad) - 5)/dr ; % this is to cound the number of additional derivatives
     
    deriv = zeros(14 + numberDer,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*0.5* (y(4) + 2 * y(3) - 3* y(2)) - 2 * Pd*y(2) + ...
        Pp * 0.5 * (y(1) - y(2)) + 0.5 * Pp * (y(3) + y(4)) *(y(1) - y(2))^2 /(y(1)^2*(1-y(1)))); % distance 1
    
    deriv(3) = 0.5 * Pm * (y(2) + y(5) + y(2) + y(5) - 4*y(3)) - ...
        2 * Pd * y(3) + 0.5 * Pp * (y(1) - y(3))* (y(1) - y(2))*(y(2)+y(5)+y(2)+ y(5))/(y(1)^2*(1-y(1))); % distance srt(2)
    
    deriv(4) = 0.5 * Pm * (y(7) + y(2) + y(5) + y(5) - 4*y(4)) - ...
        2 * Pd * y(4) + 0.5 * Pp * (y(1) - y(4))* (y(1) - y(2))*(y(7)+y(2)+y(5)+ y(5))/(y(1)^2*(1-y(1))); % distance 2
    
    deriv(5) = 0.5 * Pm * (y(8) + y(3) + y(6) + y(4) - 4*y(5)) - ...
        2 * Pd * y(5) + 0.5 * Pp * (y(1) - y(5))* (y(1) - y(2))*(y(8)+y(3)+y(6)+ y(4))/(y(1)^2*(1-y(1))); % distance sqrt(5)
    
    deriv(6) = 0.5 * Pm * (y(9) + y(5) + y(9) + y(5) - 4*y(6)) - ...
        2 * Pd * y(6) + 0.5 * Pp * (y(1) - y(6))* (y(1) - y(2))*(y(9)+y(5)+y(9)+ y(5))/(y(1)^2*(1-y(1))); % distance sqrt(8)
    
    deriv(7) = 0.5 * Pm * (y(10) + y(4) + y(8) + y(8) - 4*y(7)) - ...
        2 * Pd * y(7) + 0.5 * Pp * (y(1) - y(7))* (y(1) - y(2))*(y(10)+y(4)+y(8)+ y(8))/(y(1)^2*(1-y(1))); % distance 3     
    
    deriv(8) = 0.5 * Pm * (y(11) + y(9) + y(5) + y(7) - 4*y(8)) - ...
        2 * Pd * y(8) + 0.5 * Pp * (y(1) - y(8))* (y(1) - y(2))*(y(11)+y(9)+y(5)+ y(7))/(y(1)^2*(1-y(1)));
    
    deriv(9) = 0.5 * Pm * (y(12) + y(8) + y(6) + y(13) - 4*y(9)) - ...
        2 * Pd * y(9) + 0.5 * Pp * (y(1) - y(9))* (y(1) - y(2))*(y(12)+y(8)+y(6)+ y(13))/(y(1)^2*(1-y(1)));   
    
    deriv(10) = 0.5 * Pm * (y(7) + y(14) + y(11) + y(11) - 4*y(10)) - ...
        2 * Pd * y(10) + 0.5 * Pp * (y(1) - y(10))* (y(1) - y(2))*(y(7)+y(14)+y(11)+ y(11))/(y(1)^2*(1-y(1)));
    
    deriv(11) = 0.5 * Pm * (y(13) + y(10) + y(8) + y(15) - 4*y(11)) - ...
        2 * Pd * y(11) + 0.5 * Pp * (y(1) - y(11))* (y(1) - y(2))*(y(13)+y(10)+y(8)+ y(15))/(y(1)^2*(1-y(1)));
        
    deriv(12) = 0.5 * Pm * (y(9) + y(9) + y(14) + y(14) - 4*y(12)) - ...
        2 * Pd * y(12) + 0.5 * Pp * (y(1) - y(12))* (y(1) - y(2))*(y(9)+y(9)+y(14)+ y(14))/(y(1)^2*(1-y(1)));  
    
    deriv(13) = 0.5 * Pm * (y(11) + y(14) + y(9) + y(15) - 4*y(13)) - ...
        2 * Pd * y(13) + 0.5 * Pp * (y(1) - y(13))* (y(1) - y(2))*(y(11)+y(14)+y(9)+ y(15))/(y(1)^2*(1-y(1))); % distance sqrt(20)
    
    
    %%%% If dr =0.5 fixed, I will not vary it, then
    
    deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(16) + y(16) - 4*y(14)) - ...
        2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(16)+ y(16))/(y(1)^2*(1-y(1))); % distance 5
    
    
    %%%% be careful because actually y(16) depend on dr    FALSE!!!!!!!!!
%     if sqrt(32) < 5 + dr && sqrt(34) < 5 +dr
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(15) + y(15) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(15)+ y(15))/(y(1)^2*(1-y(1)));
%     elseif  sqrt(32) > 5 +dr
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(16) + y(16) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(16)+ y(16))/(y(1)^2*(1-y(1)));
%     else 
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(15) + y(16) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(15)+ y(16))/(y(1)^2*(1-y(1)));       
%     end       
    % now the distance is more than 5, so all will be approximated
    % uniformly 
    
   
    
    current_radius = 5;
    
    i = 15; % will start considering regular grid from this derivative
    
    % update derivatives until 1 before maximal radius
    while current_radius +dr < max_rad -dr % this is very approximate
        
        current_radius = current_radius + dr; % increase radius by dr
        
        deriv(i) = Pm * (y(i-1)+ y(i+1) - 2 *y(i))/dr^2 + 1/current_radius * (y(i+1) - y(i-1))/(2*dr) ...
            - 2 * Pd * y(i) + Pp * (y(1) - y(i))* (y(1) - y(2))*(y(i+1)+ y(i-1))/(y(1)^2*(1-y(1)));
        
        i = i+1;
    end
    
    last = 14 + numberDer;
    % periodic boundary conditions
%     deriv(last) =  Pm * (y(last-1)+ y(2) - 2 *y(last))/dr^2 + 1/current_radius * (y(2) - y(last-1))/(2*dr) ...
%             - 2 * Pd * y(last) + Pp * (y(1) - y(last))* (y(1) - y(2))*(y(2)+ y(last-1))/(y(1)^2*(1-y(1)));
    % boundary conditions correctly
    deriv(last) =  Pm * (y(last-1)+ y(last-1) - 2 *y(last))/dr^2 + 1/current_radius * (y(last-1) - y(last-1))/(2*dr) ...
            - 2 * Pd * y(last) + Pp * (y(1) - y(last))* (y(1) - y(2))*(y(last-1)+ y(last-1))/(y(1)^2*(1-y(1)));
        
end

