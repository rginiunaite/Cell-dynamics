clear all;


t(1) = 0;
tfinal = 1000;
X = [100, 0];

A(1) = 100;
B(1) = 0;
k1 = 0.01;
k2 = 1;
v = 1;
N=1000;

for i =1:N
    
a1 = k1/v*A(i)*(A(i)-1);

a2 = k2 * B(i);

a0 = a1 + a2;

a1 = a1;
a2 = a2;

r1 = rand;

tau = 1/a0 * log(1/r1);

if t(i) + tau < tfinal
    r2 = rand;
    
    if r2*a0 < a1
     A(i+1) = A(i) - 2;
     B(i+1) = B(i) + 1;
    else
     A(i+1) = A(i) + 2;
     B(i+1) = B(i) - 1;
    end
    t(i+1) = t(i)+tau;
else 
    break;
end

allB(i) = B(size(B,2))

end
 
histogram(allB)
%kazkas negerai

plot(t,A)
hold 
plot(t,B)
    