%%%% ABC for the KSA on a hexagonal lattice
clear all
 
S = load('gillespie2D_large_initial_Pp05_Pd01_40G.mat'); % load the data for a particular choice of Pp and Pd, this is our data



cA0 = S.save_density(1);

density_length = length(S.save_density); % length of the size of density

TimeReductionFactor = 1; % when death rate is increase, this has to be increased as well


dt = 0.01; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie

t = 0; % initial time

Pm = 1;

L = 3; % truncated length of the domain

n_meas = 5; % number of measurements

% assume that there are 6 measurements
values = round(density_length/(n_meas));

values_for_KSA = round(N/(n_meas-1))-1; % -1 so that the final one would be included in case N=5000

% number of simulations to apply ABC
dim = 10; % square matrix size

K = dim^4;

sum_stat = zeros(K,1);
Pp = zeros(1,K);
Pd = zeros(1,K);
% sample Pp(k) and Pd(k) from uniform distribution
for k=1:K
   Pp(k) = rand; 
   Pd(k) = rand/2;
end
%    Pp=0.1;
%    Pd=0.2;

% M = 3; % the number of replicates, should have different data sets, and most likely different initial conditions
% %%% which would be used for 

%data_log = zeros(meas +1,M);
    
    for k = 1:K
        %for j = 1:K
            %%%%%%%%%%%%%%%%%%%
           
                %%% simulate D(k) with the model parameters Pp(k) and Pd(k)
                
                y(:,1) = S.save_cA0^2 * ones(L+1,1); % initial density of two-point distribution functions 
                y(1,1) = S.save_cA0; % initial condition for 1-point distribution function

                for i = 2 : N

                    t = t + dt;

                    k1 = dynamics (t,y(:,i-1),Pm,Pp(k),Pd(k),L);
                    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,Pm,Pp(k),Pd(k),L);
                    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,Pm,Pp(k),Pd(k),L);
                    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,Pm,Pp(k),Pd(k),L);
                    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  

                    store_time_Runge(i) = t;
                end
                y_rescaled = Pp(k)/(Pp(k)-Pd(k)) * y;

                dens_KSA = y_rescaled(1,:);
                %%%%%%%%%% compute summary statistics
            sum_stat(k) = norm (S.save_density(1:values:density_length) - dens_KSA(1:values_for_KSA:N)');   

      %  end
    end
    


%%%%

% If I had a matrix, I would need this
% sum_stat = sum_stat'; % transpose, so that it would be easier to keep track of indices later when converted to a vector
% 
% sum_stat = sum_stat(:); % convert it into a vector

[sorted_sum_stat,indices] = sort(sum_stat,'ascend'); % sort in ascending order for further selection of 1%

%%%% eps will be 1% of the smallest distances

one_percent_final = round(K/100);

%small_indices = find(indices<one_percent_final); % find the indices that are below the threshold 1%
small_indices = indices(1:one_percent_final);
small_indices = indices(1:50);

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
%title(['Percentage of accepted samples = ' num2str()],'FontSize',14,'FontWeight','Normal')

function deriv = dynamics(t,y,Pm,Pp,Pd,L)
    
 
     
    deriv = zeros(L+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*1/3* (2*y(2) + 3 * y(3) - 5* y(2)) - 2 * Pd*y(2) + ...
        Pp * 1/3 * (y(1) - y(2)) + 1/3 * Pp * (2*y(2) + 3*y(3)) *(y(1) - y(2))^2 /(y(1)^2*(1-y(1)))); % distance 1
    
    
    for n=3:L
    
        deriv(n) = 1/3 * Pm * (y(n-1) + 2 * y(n) +3 * y(n+1)- 6*y(n)) - ...
            2 * Pd * y(n) + 1/3 * Pp * (y(1) - y(n))* (y(1) - y(2))*(y(n-1)+2*y(n)+3*y(n+1))/(y(1)^2*(1-y(1))); % 
    
    end
    

   deriv(L+1) = 1/3 * Pm * (y(L) + 2 * y(L+1) +3 * y(L)- 6*y(L+1)) - ...
           2 * Pd * y(L+1) + 1/3 * Pp * (y(1) - y(L+1))* (y(1) - y(2))*(y(L)+2*y(L+1)+3*y(L))/(y(1)^2*(1-y(1))); % Periodic boundary conditions
     
%   deriv(L+1) = 1/3 * Pm * (y(L) + 2 * y(L+1) +3 * y(1)^2- 6*y(L+1)) - ...
%           2 * Pd * y(L+1) + 1/3 * Pp * (y(1) - y(L+1))* (y(1) -...
%           y(2))*(y(L)+2*y(L+1)+3*y(1)^2)/(y(1)^2*(1-y(1))); % cells far
%           % apart independent
    
end