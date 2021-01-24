function theta = WLSQ( samples, p,theta0, d0 )
N = length(samples) - p;
H = zeros(N,p);
for m = 1:N
    H(m,:) = - samples(m + p - 1: -1 : m);
end
W = zeros(N,N);
a4 = 0.5;
epsilon = 1e-6;
continue_func = 1;
num_iter = 0;
while continue_func && num_iter<4
    num_iter = num_iter+1;
    suma = 0;
    
    for m = 1: N
        eps = (samples(m+p)-H(m,:)*theta0)/d0;
        if samples(m+p)~=H(m,:)*theta0
            if ( abs(eps) <= a4*pi ) % Andrews nonlinearity
                ksi = sin(eps/a4);
            else
                ksi = 0;
            end
            suma = suma +  d0^2*ksi^2;
            W(m,m) =ksi/eps;
        else
            W(m,m) = 1;
        end
    end
    
    d1 = sqrt(suma/N/0.561722); % Expectation is different than in Dutter agorithm because of the use of Andrews nonlinearity
    theta = inv(H'*W*H)*H'*W*samples(1+p:N+p);
    
    continue_func = 0;
    for m=1:p
        if ~(abs(theta(m)-theta0(m))<abs(epsilon*theta0(m)))
            continue_func = 1;
        end
    end
    theta0 =theta;
    d0 = d1;
end