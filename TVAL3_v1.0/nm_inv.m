function nm_inv

clear; close all;
path(path,genpath(pwd));
fullscreen = get(0,'ScreenSize');

% problem size
ratio = 1;
sidelength = 32;
n = sidelength^2;
m = round(ratio*n);
k = floor(m/10);

% original staircase signal
xs = zeros(n,1);
p = randperm(n); p = sort(p(1:k-1)); p = [1 p n];
for sct = 1:k
    xs(p(sct):p(sct+1)) = rand-0.5 + (rand-0.5)*1i;
end
nrmxs = norm(xs,'fro');

H = hadamard(n);
A = H(1:m,:);

% observation
b = A*xs;
bavg = mean(abs(b));

% add noise
sigma = 0.05;  % noise std
noise = randn(m,1);
b = b + sigma*bavg*noise;

% reconstruction
t = cputime;
x = 1/n*A*b;
t = cputime - t;
rerr = norm(x-xs,'fro')/nrmxs;

% plotting
figure('Name','TVAL3','Position',...
    [fullscreen(1) fullscreen(2) fullscreen(3) fullscreen(4)]);
subplot(211); set(gca,'fontsize',16)
plot(1:n,real(xs),'r.-',1:n,real(x),'b-');
title(sprintf('Real Part         Noise: %2.1f%%,   Rel-Err: %4.2f%%,   CPU: %4.2fs',sigma*100,rerr*100,t))
subplot(212); set(gca,'fontsize',16)
plot(1:n,imag(xs),'r.-',1:n,imag(x),'b-');
title('Image Part                                                                             ')



