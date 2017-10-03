%%%% ABC for logistic

 
S = load('gillespie2D_Pp05_Pd01_40G.mat'); % load the data for a particular choice of Pp and Pd, this is our data

cA0 = S.save_density(1);

density_length = length(S.save_density); % length of the size of density

meas = 8; % assume there are meas measurements of real data

values = round(density_length/5);

% number of simulations to apply ABC
K = 50;

sum_stat = zeros(K,K);
Pp = zeros(1,K);
Pd = zeros(1,K);
% sample Pp(k) and Pd(k) from uniform distribution
for k=1:K
   Pp(k) = rand; 
   Pd(k) = rand/2;
end

% M = 3; % the number of replicates, should have different data sets, and most likely different initial conditions
% %%% which would be used for 

data_log = zeros(meas +1,M);
    
    for k = 1:K
        for j = 1:K
            %%%%%%%%%%%%%%%%%%%
           
                %%% simulate D(k) with the model parameters Pp(k) and Pd(k)
                cA0_bar = (Pp(k))/(Pp(k) - Pd(j))*cA0;

                f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

                for i = 1 : length(S.save_time)
                       y_log(i) = f(S.save_time(i)); 
                end
            
                %%%%%%%%%% compute summary statistics
                sum_stat(k,j) = norm (S.save_density(1:values:density_length) - y_log(1:values:density_length));   

        end
    end
    
    




%%%%

sum_stat = sum_stat'; % transpose, so that it would be easier to keep track of indices later

sum_stat = sum_stat(:); % convert it into a vector

[sorted_sum_stat,indices] = sort(sum_stat,'ascend');

%%%% eps will be 1% of the smallest distances

one_percent_final = round(K*K/100);

small_indices = find(indices<one_percent_final); % find the indices that are below the threshold 1%


% extract the matrix entry number which corresponds to different values of
% Pp and Pd
for j = 1 : length(small_indices)
   if mod(small_indices(j),K)==0
      Pp_ind(j) = small_indices(j)/K;
      Pd_ind(j) = K;
     
   else 
     Pp_ind(j) = floor(small_indices(j)/K)+1;
     Pd_ind(j) = mod(small_indices(j),K);
   end
   
end

sample = sorted_sum_stat(1:one_percent_final);

figure
histogram(Pp(Pp_ind),5)

figure 
histogram(Pd(Pd_ind),5)
