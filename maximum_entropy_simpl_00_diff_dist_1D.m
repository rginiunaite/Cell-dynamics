% combination of Gillespie, logistic and Runge-Kutta with
% maximum entropy approximation in 1D, using fsolve for step 3 in maximal
% entropy scheme

clear all

% Parameters

% lattice dimension
dim = 1;

%  length of the domain
L = 1000;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

Pd = 0.1; % death rate per unit time

S = load('gillespie1D_Pp05_Pd01_40G_long.mat'); % load the data for a particular choice of Pp and Pd, note that L=100 is fixed

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

dt = 0.1; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie


%%%%%%%%%%%%%%%%%%%%

%MAXIMUM ENTROPY METHOD

%%%%%%%%%%%%%%%%%%%%

max_iter_gamma = 10; % guessed number of iterations to find gamma that converges
tol = 0.005;

gammaAA = ones(L,max_iter_gamma); % size of gamma AA
gammaA0 = ones(L,max_iter_gamma); % size of gamma A0
gamma00 = ones(L,max_iter_gamma); % size of gamma A0
% 
% gamma = 5*ones(2*L,max_iter_gamma);

% initial conditions for gamma, first L/2 correspond to gamma(A,A), the other L/2, gamma(A,0)


    for i = 1:L-1

         gammaAA(i) = (1-S.save_cA0^1)^1;%  % worked for these
         gammaA0(i) = S.save_cA0^1;
         gamma00(i) = S.save_cA0^1; % chosen without any good reason
    end
    
% Runge Kutta method to solve ODEs

%%%%%%%Initial conditions
y(:,1) = S.save_cA0^2 * ones(L+1,1); % initial density of two-point distribution functions
y(1,1) = S.save_cA0; % initial condition for 1-point distribution function

count = 1;


t = 0; % start time
kiekis = zeros(1,N); % keep track how many times go into the while loop


for i = 2 : N

     t = t + dt;

     count = 1; % always start from the first gamma
     
     k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));
     k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd,gammaAA(:,count),gammaA0(:,count));

     y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  

     % use fsolve to find new gamma
     final = size(y(1,:),2);

     options = optimoptions('fsolve','MaxIterations',10000, 'FunctionTolerance',1e-2,'MaxFunctionEvaluations',10000,'StepTolerance',1.0e-06);

     [xAA,fval] = fsolve(@(gamma_unk)step3AA(gamma_unk,gammaA0(:,count),y(2:L+1,final),L,y(1,final)),gammaAA(:,count),options);
     [xA0,fval] = fsolve(@(gamma_unk)step3A0(gamma_unk,gammaAA(:,count),gamma00(:,count),y(2:L+1,final),L,y(1,final)),gammaA0(:,count),options);
     [x00,fval] = fsolve(@(gamma_unk)step300(gamma_unk,gammaA0(:,count),gammaAA(:,count),y(2:L+1,final),L,y(1,final)),gamma00(:,count),options);

     gammaAA(:,count+1) = xAA;
     gammaA0(:,count+1) = xA0;   
     gamma00(:,count+1) = x00;

    while norm (gammaA0(:,count)-gammaA0(:,count+1)) > tol || norm(gammaAA(:,count)-gammaAA(:,count+1)) > tol
        count = count + 1; % update new step
                
        [xAA,fval] = fsolve(@(gamma_unk)step3AA(gamma_unk,gammaA0(:,count),y(2:L+1,final),L,y(1,final)),gammaAA(:,count),options);
        [xA0,fval] = fsolve(@(gamma_unk)step3A0(gamma_unk,gammaAA(:,count),gamma00(:,count),y(2:L+1,final),L,y(1,final)),gammaA0(:,count),options);
        [x00,fval] = fsolve(@(gamma_unk)step300(gamma_unk,gammaA0(:,count),gammaAA(:,count),y(2:L+1,final),L,y(1,final)),gamma00(:,count),options);
        
        % update new value of gamma
        gammaAA(:,count+1) = xAA; 
        gammaA0(:,count+1) = xA0;     
        gamma00(:,count+1) = x00;
        
        kiekis(i) = i; % just check when it goes into the loop
    end
    
    % set gamma to the final value
    gammaA0(:,1) = gammaA0(:,count); 
    gammaAA(:,1) = gammaAA(:,count); 
    gamma00(:,1) = gamma00(:,count);
  
    store_time_Runge(i) = t;

end


rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = Pp/(Pp-Pd) * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
% h_legend = legend('Gillespie 1D','Logistic growth','Maximum entropy');
% set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)
ax = gca;
ax.YAxis.TickLabelFormat = '%,.1f';

function deriv = dynamics(t,y,Ln,Pm,Pp,Pd,gammaAA,gammaA0)


    deriv  = zeros (Ln+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*(y(3) - y(2)) - 2* Pd * y(2) + ...
        Pp * (y(1) - y(2)) + Pp * (gammaAA(2)*gammaA0(1)^2)); % dynamics of two-point distribution function (distance 1)
    
 % dynamics of all two-point functions from distance 2 to L-1, note that
 % indix i+1 correspond to distance i, therefore, gamma indices are smaller

     deriv(3:Ln) = Pm * (y(2:Ln-1) + y(4:Ln+1) - 2*y(3:Ln)) - 2 * Pd * y(3:Ln)+ ...
        Pp*(gammaA0(2:Ln-1).*(gammaA0(1)*ones(Ln-2,1)).*(gammaAA(1:Ln-2)+gammaAA(3:Ln)));

    
   deriv(Ln+1) = (Pm*(y(Ln) + y(1) - 2* y(Ln+1)) - 2 * Pd * y(Ln+1) +...
           Pp * (gammaA0(Ln)*gammaA0(1)*(gammaAA(Ln-1)+gammaAA(Ln-1)))); % for the last one,
%         deriv(Ln+1) = (Pm*(y(Ln) + y(1) - 2* y(Ln+1)) - 2 * Pd * y(Ln+1) +...
%             Pp * (gammaA0(Ln)*gammaA0(1)*(gammaAA(Ln-1)+gammaAA(1)))); % for the last one,
end


% this function corresponds to the non-linear function in step 3 that
% needs to be solved for gamma^{(i+1)}

%% assume the only constraint is for (A_L, 0_{L+n}, X_{L+1})

function y = step3AA(gammaAA,gammaA0,p2,Ln,cA)
    
    y = zeros(1,Ln);
       
    %y(1:Ln-1) = p2(1:Ln-1) - gammaAA(1:Ln-1).*(gammaAA(2:Ln).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gammaA0(2:Ln));   
    %y(Ln) = p2(Ln) - gammaAA(Ln)*(gammaAA(1)*gammaAA(1)+gammaA0(1)*gammaA0(1));
    
    %y(1:Ln) = p2(1:Ln) - gammaAA(1:Ln).*(gammaAA(1:Ln).*(gammaAA(1)*ones(Ln,1))+gammaA0(1)*ones(Ln,1).*gammaA0(1:Ln));   
   
    % works for this
    y(1:Ln) = p2(1:Ln) - gammaAA(1:Ln).*(gammaAA(1:Ln).*(gammaAA(1)*ones(Ln,1))+gammaA0(1)*ones(Ln,1).*gammaA0(1:Ln));   
    
    
    % smaller distance to the neighbour
%     y(2:Ln) = p2(2:Ln) - gammaAA(2:Ln).*(gammaAA(1:Ln-1).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gammaA0(1:Ln-1));   
%     y(1) = p2(1) - gammaAA(1)*(gammaAA(Ln-1)*gammaAA(1)+gammaA0(Ln-1)*gammaA0(1));
%     
%     y(Ln+1:2*Ln-1) = p2(1:Ln-1) - gammaAA(1:Ln-1).*(gammaAA(2:Ln).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gammaA0(2:Ln));
%     y(2*Ln) = p2(Ln) - gammaAA(Ln)*(gammaAA(Ln-1)*gammaAA(1)+gammaA0(1)*gammaA0(Ln-1));
    
    
    % another attempt to make it more correct
    % y(2:Ln) = p2(2:Ln) - gammaAA(2:Ln).*(gammaAA(1:Ln-1).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gammaA0(1:Ln-1));   
    % y(1) = p2(1) - gammaAA(1)*(gammaAA(Ln-1)*gammaAA(1)+gammaA0(Ln-1)*gammaA0(1));   
   
    
end


function y = step3A0(gammaA0,gammaAA,gamma00,p2,Ln,cA)
    
    y = zeros(1,Ln);
           
   
   % y(1:Ln-1) = cA*ones(Ln-1,1) - p2(1:Ln-1) - gammaA0(1:Ln-1).*(gammaAA(2:Ln).*(gammaA0(1)*ones(Ln-1,1))+gamma00(1)*ones(Ln-1,1).*gammaA0(2:Ln));   
   % y(Ln) = cA - p2(Ln) - gammaA0(Ln)*(gammaAA(1)*gammaA0(1)+gamma00(1)*gammaA0(1));
  
    %y(1:Ln) = cA*ones(Ln,1) - p2(1:Ln) - gammaA0(1:Ln).*(gammaAA(1:Ln).*(gammaA0(1)*ones(Ln,1))+gamma00(1)*ones(Ln,1).*gammaA0(1:Ln));   
 
    
    %  works for this 
    y(1:Ln) = cA*ones(Ln,1) - p2(1:Ln) - gammaA0(1:Ln).*(gammaA0(1:Ln).*(gammaAA(1)*ones(Ln,1))+gammaA0(1)*ones(Ln,1).*gamma00(1:Ln));

   % smaller distance to the neighbour  
%    y(2:Ln) = cA*ones(Ln-1,1) - p2(2:Ln) - gammaA0(2:Ln).*(gammaA0(1:Ln-1).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gamma00(1:Ln-1));   
%    y(1) = cA - p2(1) - gammaA0(1)*(gammaA0(Ln-1)*gammaAA(1)+gammaA0(1)*gamma00(Ln-1));   
%    
%     y(Ln+1:2*Ln-1) = cA*ones(Ln-1,1) - p2(1:Ln-1) - gammaA0(1:Ln-1).*(gammaA0(2:Ln).*(gammaAA(1)*ones(Ln-1,1))+gammaA0(1)*ones(Ln-1,1).*gamma00(2:Ln));  
%     y(2*Ln) = cA - p2(Ln) - gammaA0(Ln)*(gammaA0(Ln-1)*gammaAA(1)+gammaA0(1)*gamma00(Ln-1));

end


function y = step300(gamma00,gammaA0,gammaAA,p2,Ln,cA)

    y = zeros(1,Ln);


 %   y(1:Ln-1) = ones(Ln-1,1) - 2 * cA*ones(Ln-1,1) + p2(1:Ln-1) - gamma00(1:Ln-1).*(gammaA0(2:Ln).*(gammaA0(1)*ones(Ln-1,1))+gamma00(1)*ones(Ln-1,1).*gamma00(2:Ln));   
 %   y(Ln) = 1 - 2*cA + p2(Ln) - gamma00(Ln)*(gammaA0(1)*gammaA0(1)+gamma00(1)*gamma00(1));
    
 %works for this
    y(1:Ln) = ones(Ln,1) - 2 * cA*ones(Ln,1) + p2(1:Ln) - gamma00(1:Ln).*(gammaA0(1:Ln).*(gammaA0(1)*ones(Ln,1))+gamma00(1)*ones(Ln,1).*gamma00(1:Ln)); 
    
   %smaller distance neighbour 
%     y(2:Ln) = ones(Ln-1,1) - 2 * cA*ones(Ln-1,1) + p2(2:Ln) - gamma00(2:Ln).*(gammaA0(1:Ln-1).*(gammaA0(1)*ones(Ln-1,1))+gamma00(1)*ones(Ln-1,1).*gamma00(1:Ln-1));   
%     y(1) = 1 - 2*cA + p2(1) - gamma00(1)*(gammaA0(Ln-1)*gammaA0(1)+gamma00(1)*gamma00(Ln-1)); % WRONG, it does not exist in this setting, I always subtract
%     
%     y(Ln+1:2*Ln-1) = ones(Ln-1,1) - 2 * cA*ones(Ln-1,1) + p2(1:Ln-1) - gamma00(1:Ln-1).*(gammaA0(2:Ln).*(gammaA0(1)*ones(Ln-1,1))+gamma00(1)*ones(Ln-1,1).*gamma00(2:Ln));
%     y(2*Ln) = 1 - 2*cA + p2(Ln) - gamma00(Ln)*(gammaA0(Ln-1)*gammaA0(1)+gamma00(1)*gamma00(Ln-1));

end
