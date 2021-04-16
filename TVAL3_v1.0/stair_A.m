function stair_wht

clear; close all;
path(path,genpath(pwd));
fullscreen = get(0,'ScreenSize');

% problem size
ratio = .6;
sidelength = 16;
n = sidelength^2;
m = round(ratio*n);
k = floor(m/10);

% original staircase signal
xs = zeros(n,1);
p = randperm(n); p = sort(p(1:k-1)); p = [1 p n];
for sct = 1:k
    xs(p(sct):p(sct+1)) = rand + rand*1i;
end
nrmxs = norm(xs,'fro');

% generate measurement matrix
p = randperm(n);
picks = p(1:m);
for ii = 1:m
    if picks(ii) == 1
        picks(ii) = p(m+1);
        break;
    end
end
perm = randperm(n); % column permutations allowable
A = @(x,mode) dfA(x,picks,perm,mode);

% observation
b = A(xs,1);
bavg = mean(abs(b));

% add noise
sigma = 0.04;  % noise std
noise = randn(m,1);
b = b + sigma*bavg*noise;

% set the optional paramaters
clear opts
opts.mu = 2^5;
opts.beta = 2^5;
opts.mu0 = 2^1;
opts.beta0 = 2^1;
opts.tol = 1E-4;
opts.maxit = 600;
opts.TVnorm = 1;

% reconstruction
t = cputime;
[x, out] = TVAL3(A,b,n,1,opts);
x = x - min(x(:));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dfA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = dfA(x,picks,perm,mode)
switch mode
    case 1
        y = A_fWH(x,picks,perm);
    case 2
        y = At_fWH(x,picks,perm);
    otherwise
        error('Unknown mode passed to f_handleA!');
end