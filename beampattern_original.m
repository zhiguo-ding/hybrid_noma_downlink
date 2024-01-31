%courtesy to Zhaolin Wang @ QMUL
clear all
%close all

N =513; % number of antennas
f = 28e9; % 28 GHz
c = 3e8; % speed of light
lambda = c/f; % wavelength
d = lambda / 2; % antenna spacing
D = (N-1)*d; %apture or size
dis_Ray = 2*D^2/lambda
    
r = 5; % distance
theta = 45 / 180 * pi; % direction
w = beamfocusing(r, theta, N, d, f)/sqrt(N); % beamforming vector

        
% calculate beampattern
m = 500;
area = 100; area2=0;%area/2;
X = linspace(area2, area, m);
Y = linspace(area2, area, m);
[X, Y] = meshgrid(X, Y);
[theta_all, r_all] = cart2pol(X, Y);

P = zeros(length(r_all), length(theta_all));
P_ff = zeros(length(r_all), length(theta_all));
for i = 1:m
    for j = 1:m
        a = beamfocusing(r_all(i,j), theta_all(i,j), N, d, f)/sqrt(N);
        P(i,j) = abs(a' * w)^2;
 
    end
    
end
 

% plot beampattern
[X, Y] = pol2cart(theta_all, r_all);
figure; hold on; colormap jet; colorbar;
mesh(X,Y,P); 
view([90,-90]);
xlim([area2,area]); ylim([area2,area]);
xlabel('x (m)'); ylabel('y (m)');

 
