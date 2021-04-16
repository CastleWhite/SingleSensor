function stair_wht

clear; close all;
path(path,genpath(pwd));
fullscreen = get(0,'ScreenSize');

% problem size
ratio = 0.7;
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

% generate measurement matrix
% p = randperm(n);
% picks = p(1:m);
% for ii = 1:m
%     if picks(ii) == 1
%         picks(ii) = p(m+1);
%         break;
%     end
% end
% perm = randperm(n); % column permutations allowable
H = hadamard(n);
A = H(1:m,:);

% observation
b = A*xs;
bavg = mean(abs(b));

% add noise
sigma = 0.02;  % noise std
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


plotting = 0;
if plotting
    figure(2);
    subplot(241); plot(out.lam1); title('\_al: ||w||');
    subplot(242); plot(out.lam2); title('\_al: ||Du-w||^2');
    subplot(243); plot(out.lam3); title('\_al: ||Au-f||^2');
    subplot(244); plot(abs(out.obj),'b-'); title('\_al: objective values');
    subplot(245); plot(out.res); title('\_al: residue');
    subplot(246); plot(abs(out.tau)); title('\_al: steplenths');
    subplot(247); plot(out.itrs); title('\_al: inner iterations');
    subplot(248); plot(abs(out.C),'r-'); title('\_al: reference vlaues');
    
    figure(3);
        semilogy(1:length(out.lam1),out.lam1,'b*:',1:length(out.lam2),sqrt(out.lam2),'rx:',...
        1:length(out.lam3),sqrt(out.lam3),'g.--', 1:length(out.f),sqrt(out.f),'m+-');
    legend('lam1(||w||_1)','lam2(||D(d_tu)-w||_2)','lam3(||Au-b||_2)','obj function');
end
