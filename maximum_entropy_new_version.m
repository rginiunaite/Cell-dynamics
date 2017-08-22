% combination of Gillespie, logistic and Runge-Kutta with
% maximum entropy approximation in 1D, using fsolve for step 3 in maximal
% entropy scheme

clear all

% Parameters

% lattice dimension
dim = 1;

% half of the length of the domain
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

Pd = 0.1; % death rate per unit time

S = load('gillespie2D_large_initial_Pp05_Pd01.mat'); % load the data for a particular choice of Pp and Pd, note that L=100 is fixed

TimeReductionFactor = 2; % when death rate is increase, this has to be increased as well

figure 

plot(S.save_time,S.save_density,'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
%h_legend = legend('1st component','2nd component','Total');
%set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)


%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

%cA0 = Q_initial/L; % initial density
cA0_bar = (Pp)/(Pp - Pd)*S.save_cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(S.save_time)
       y_log(i) = f(S.save_time(i)); 
end

%%% plot logistic solution 

hold on
plot(S.save_time,y_log,'LineWidth', 3)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)

dt = 0.01; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie



%%%%%%%%%%%%%%%%%%%%

%MAXIMUM ENTROPY METHOD

%%%%%%%%%%%%%%%%%%%%

max_iter_gamma = 10; % guessed number of iterations to find gamma that converges

gammaAA = ones(L,max_iter_gamma); % size of gamma AA
gammaA0 = ones(L,max_iter_gamma); % size of gamma A0
% 
% gamma = 5*ones(2*L,max_iter_gamma);

% initial conditions for gamma, first L/2 correspond to gamma(A,A), the other L/2, gamma(A,0)


    for i = 1:L-1
         gammaAA(i) = 3.89*S.save_cA0^3; % worked for these, old data
         gammaA0(i) = 3.89*S.save_cA0*(1-S.save_cA0)^2; 
%          gammaAA(i) = 4.13*S.save_cA0^3; % worked for these
%          gammaA0(i) = 4.13*S.save_cA0*(1-S.save_cA0)^2;          
 
          gammaAA(i) = 1.0122*(1-S.save_cA0); % worked for these
          gammaA0(i) = 1.0122*S.save_cA0;   

%          gammaAA(i) = 1-S.save_cA0^1; % worked for these
%          gammaA0(i) = S.save_cA0;%(1-S.save_cA0)^1; 


%%%%%%%%%%%%%% try to choose initial conditions based on the KSA
%%%%%%%%%%%%%% approximation of triple, not a good idea

%     %%% doubles from Power 1 approx, they are similar to Gillespie
%     Data = load('gillespie1D_large_initial_Pp05_Pd01_doubles.mat');
%     pairs = Data.save_doubles;
%     den = S.save_density(length(S.save_density));
%     
%     gammaAA(i) = (pairs(i)+pairs(i+2))/(1- den);
% 
%     gammaA0(i) = sqrt((den-pairs(i+1))*(den-pairs(2)))/den;


    end
    
%     gammaAA(L) = (pairs(L+1)+pairs(L))/(1- den);
% 
%     gammaA0(L) = sqrt((den-pairs(L+1))*(den-pairs(2)))/den;  
    
    
    

%     gamma(1:L,1) = gammaAA(:,1);
%     gamma(L+1:2*L,1) = gammaA0(:,1);

% calculate the second one separately from the rest, so that could one
% could keep track of changes in the norm correctly


% Runge Kutta method to solve ODEs

%%%%%%%Initial conditions
y(:,1) = S.save_cA0^2 * ones(L+1,1); % initial density of two-point distribution functions
y(1,1) = S.save_cA0; % initial condition for 1-point distribution function

count = 1;


%%%%%%%%%%%%
% calculate the second gamma separately
%%%%%%%%%


t = 0; % start time


for i = 2 : N

     t = t + dt;

     k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));

     y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  

     store_time_Runge(i) = t;
end


% use fsolve to find new gamma
final = size(y(1,:),2);

options = optimoptions('fsolve','MaxIterations',10000, 'FunctionTolerance',1e-3,'MaxFunctionEvaluations',10000,'StepTolerance',1.0000e-06);


[xAA,fval] = fsolve(@(gamma_unk)step3AA(gamma_unk,y(2:L+1,final),L,y(1,final)),gammaAA(:,count),options);
[xA0,fval] = fsolve(@(gamma_unk)step3A0(gamma_unk,y(2:L+1,final),L,y(1,final)),gammaA0(:,count),options);

gammaAA(:,count+1) = xAA;
gammaA0(:,count+1) = xA0;   


% gamma(1:L,count+1) = gammaAA(:,count+1);
% gamma(L+1:2*L,count+1) = gammaA0(:,count+1);


count = count +1; % update count, it is 2 now


for i = 2 : N

     t = t + dt;

     k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));

     y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  

     store_time_Runge(i) = t;
end



tol = 3; % tolerance for the convergence of gamma


while norm (gammaA0(:,count)-gammaA0(:,count-1))>tol 
    

    % Runge Kutta method to solve ODEs to find density and pair densities

    t = 0; % start time


    for i = 2 : N

        t = t + dt;

        k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
        k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
        k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
        k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
        y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  

        store_time_Runge(i) = t;
    end


    % use fsolve to find new gamma
    final = size(y(1,:),2);
    
     
    [xAA,fval] = fsolve(@(gamma_unk)step3AA(gamma_unk,y(2:L+1,final),L,y(1,final)),gammaAA(:,count),options);
    [xA0,fval] = fsolve(@(gamma_unk)step3A0(gamma_unk,y(2:L+1,final),L,y(1,final)),gammaA0(:,count),options);
   
    gammaAA(:,count+1) = xAA; % update new value of gamma
    gammaA0(:,count+1) = xA0; % update new value of gamma

    
%     gamma(1:L,count+1) = gammaAA(:,count+1);
%     gamma(L+1:2*L,count+1) = gammaA0(:,count+1);
    
    count = count + 1; % update new step
end



rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = Pp/(Pp-Pd) * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
h_legend = legend('Gillespie 2D','Logistic growth','Maximum entropy');
set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['L = ' num2str(L) ', P_m = ' num2str(Pm) ', P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)


function deriv = dynamics(t,y,Ln,Pm,Pp,Pd,gammaAA,gammaA0)


    deriv  = zeros (Ln+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*(y(3) - y(2)) - 2* Pd * y(2) + ...
        Pp * (y(1) - y(2)) + Pp * (gammaAA(2)*gammaA0(1)^2)); % dynamics of two-point distribution function (distance 1)
    
 % dynamics of all two-point functions from distance 2 to L-1, note that
 % indix i+1 correspond to distance i, therefore, gamma indices are smaller
    for i = 3:Ln
        deriv(i) = (Pm*(y(i-1) + y(i+1) - 2* y(i)) - 2 * Pd * y(i) +...
            Pp * (gammaA0(i-1)*gammaA0(1)*(gammaAA(i-2) + gammaAA(i)))); 
    end

    deriv(Ln+1) = (Pm*(y(Ln) + y(2) - 2* y(Ln+1)) - 2 * Pd * y(Ln+1) +...
            Pp * (gammaA0(Ln)*gammaA0(1)*(gammaAA(Ln-1)+gammaAA(1)))); % for the last one, since we have
        % periodic boundary condition L+1 neighbour is L-1 neighbour, i.e. y(2)
    
end


% this function corresponds to the non-linear function in step 3 that
% needs to be solved for gamma^{(i+1)}

function y = step3AA(gamma,p2,Ln,cA)
    
    add = zeros(1,Ln);
    y = zeros(1,Ln);
        
    for l=1:Ln
      
        for j = 1:Ln
            if j ~= l
                add(l) = add(l) + gamma(j)*gamma(abs(j-l));
            else % if distance is zero, gamma is zero
                add(l) = add(l); 
            end
        end
        
        y(l) = p2(l) - gamma(l)*add(l);     
    end

    
end


function y = step3A0(gamma,p2,Ln,cA)
    
    add = zeros(1,Ln);
    y = zeros(1,Ln);
        
    for l=1:Ln
      
        for j = 1:Ln
            if j ~= l
                add(l) = add(l) + gamma(j)*gamma(abs(j-l));
            else % if distance is zero, gamma is zero
                add(l) = add(l); 
            end
        end
        
        y(l) = cA - p2(l) - gamma(l)*add(l);     
    end

    
end




