%%%% ABC for logistic

 
S = load('gillespie2D_Pp05_Pd01_40G_hexa.mat'); % load the data for a particular choice of Pp and Pd, this is our data

cA0 = S.save_density(1);

density_length = length(S.save_density); % length of the size of density

meas = 5; % assume there are meas measurements of real data

values = round(density_length/5);

% number of simulations to apply ABC
K = 10^4;

sum_stat = zeros(meas+1,1);
Pp = zeros(1,K);
Pd = zeros(1,K);
% sample Pp(k) and Pd(k) from uniform distribution
for k=1:K
   Pp(k) = rand; 
   Pd(k) = rand/2;
end



    
    for k = 1:K
        %for j = 1:K
            %%%%%%%%%%%%%%%%%%%
           
                %%% simulate D(k) with the model parameters Pp(k) and Pd(k)
                cA0_bar = (Pp(k))/(Pp(k) - Pd(k))*cA0;

                f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

                for i = 1 : length(S.save_time)
                       y_log(i) = f(S.save_time(i)); 
                end
            
                %%%%%%%%%% compute summary statistics
                sum_stat(k) = norm (S.save_density(1:values:density_length) - y_log(1:values:density_length)');   

        %end
    end
    
    

[sorted_sum_stat,indices] = sort(sum_stat,'ascend'); % sort in ascending order for further selection of 1%

%%%% eps will be 1% of the smallest distances

one_percent_final = round(K/100);

%small_indices = find(indices<one_percent_final); % find the indices that are below the threshold 1%
small_indices = indices(1:one_percent_final);
%small_indices = indices(1:50);

%%%%%%%
% Nonsense
% extract the matrix entry number which corresponds to different values
% of%%%%%%%%%%% 
% Pp and Pd

% 
% for j = 1 : length(small_indices)
%    if mod(small_indices(j),dim)==0
%       Pp_ind(j) = small_indices(j)/dim;
%       Pd_ind(j) = dim;
%      
%    else 
%      Pp_ind(j) = floor(small_indices(j)/dim)+1;
%      Pd_ind(j) = mod(small_indices(j),dim);
%    end
%    
% end

% sample = sorted_sum_stat(1:one_percent_final);

figure
scatter(Pp(small_indices),Pd(small_indices),'LineWidth', 2)
xlabel('P_p')
ylabel('P_d')
set(gca,'linewidth',3)
%set(gca,'FontWeight','bold')
set(gca,'FontSize',36)
ax = gca;
ax.YAxis.TickLabelFormat = '%,.2f';
ax.XAxis.TickLabelFormat = '%,.1f';
