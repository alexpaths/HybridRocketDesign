% Alexander Athougies
% Senior Design
% 
% plots Mohr's Circle
% 
% INPUTS
%  1) sigma_max
%  2) sigma_min
%  3) c = color
%  
function Mcirc(sigma_max, sigma_min, c)

hold on
center = (sigma_max + sigma_min)/2;
radius = sigma_max - center;
N = 100;
circle([center,0],radius,N,c);
grid on
axis square

end