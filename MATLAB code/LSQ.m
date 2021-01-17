function theta0 = LSQ(samples,p) 
N = length(samples) - p;
H = zeros(N,p);
for m = 1: N
   H(m,:) = - samples(m + p - 1: -1 : m);
end
theta0 = (H'*H)\H'*samples(p+1:N+p);