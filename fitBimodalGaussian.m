function [fitCoeff, gof] = fitBimodalGaussian(respAmp, baseAmp)

% input directions, average response during drift grating
x = 0 : 45 : 335;
y = respAmp;


% initialize start point of all the parameters for the optimzation process
pref_d = find(respAmp == max(respAmp));
startP(1) = x(pref_d); % presented direction with the maximum amplitude as the initialized point for the optimal preferred direction
startP(2) = baseAmp; % baseline firing rate estimated from the lowest 10% of dF for each cell
oppo_d = mod(pref_d + 4, 8); % 
if oppo_d == 0
    oppo_d = 8;
end
startP(3) = y(oppo_d);
startP(4) = y(pref_d);
startP(5) = 40;


% set lower limit for all the parameters
lowerLimit(1) = startP(1) - 25;
lowerLimit(2) = startP(2) - 0.02;
lowerLimit(3) = startP(3);
lowerLimit(4) = startP(4);
lowerLimit(5) = 10;


% set upper limit for all the parameters
upperLimit(1) = startP(1) + 25;
upperLimit(2) = 0.1;
upperLimit(3) = startP(3) + 0.03;
upperLimit(4) = startP(4) + 0.06;
upperLimit(5) = 100;


% bimodal gaussian with same standard deviation
gaussEqn = 'rc + rp * exp(-(min([abs(x-pref), abs(360-x+pref), abs(x-pref+360)], [], 2)/(sqrt(2)*sigma))^2) + ro * exp(-(min([abs(x-pref+180), abs(180-x+pref), abs(x-pref+540)], [], 2)/(sqrt(2)*sigma))^2)';
[fitCoeff, gof] = fit(x', y', gaussEqn, 'StartPoint', startP, 'Lower', lowerLimit, 'Upper', upperLimit);


