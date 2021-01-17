function [theta, d] = Dutter(S,p,theta0)
N = length(S) - p;
d0 = (median(abs(S-median(S))))/0.6745;
H = zeros(N,p);

for m = 1: N
   H(m,:) = - S(m + p - 1: -1 : m);
end

n= zeros(N,1);
delta = zeros(N,1);
ksi_D=zeros(N,N);
k = 1.5; a=5.6; a2=3;
q = min(1/(2*erf(k)),1.9);
epsilon = 1e-4;

contionue_func=1;
num_iter=0;
while (contionue_func && num_iter<10)
    num_iter=num_iter+1;
    
    % Step 1
    for m = 1:N
        n(m) = S(m+p) - H(m,:)*theta0;
    end
    
    % Step 2
    suma = 0;
    for m = 1: N
        eps = n(m)/d0;
        if(abs(eps) <= k) % Dutter
            ksi1=eps;
        else
            ksi1=k*sign(eps);
        end
        if abs(eps)<a % Tukey
                ksi2=eps*(1-(eps/a)^2);
            else
                ksi2=0;
        end
        if(abs(eps) <= a2) % truncation
            ksi3=eps;
        else
            ksi3=0;
        end
        ksi=ksi1;
        suma = suma + ksi^2;
    end
    d1 = sqrt(d0^2*suma/N/0.7785);
    
    % Step 3
    for m = 1:N
        eps = n(m)/d1;
        if(eps > k)
            delta(m) = k*d1;
        elseif(eps < -k)
            delta(m) = -k*d1;
        else
            delta(m) = n(m);
        end
    end
    
    % Step 4
    d_theta = inv(H'*H)*H'*delta;
    
    % Step 5
    theta = theta0 + q*d_theta;
    
    % Step 6  
    contionue_func = 0;
    for m=1:p
       if ~(abs(theta(m)-theta0(m))<abs(epsilon*theta0(m)) && abs(d1-d0)<epsilon*abs(d0))
           contionue_func=1;           
       end
    end
    d0=d1;
    theta0=theta;
end
d=d0;

