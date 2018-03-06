function [wei] = weibull(x,xdata) 
%LET is greater than the threshold in LET therefore we use the second equation
%wei = x(1).*(1-exp(-((BLAH-x(4)))./x(2)).^x(3)); %weibull function equation where 
% x(1) is C_s, x(2) is L_0, x(3) is W and x(4) is S

E_b = [0.6, 0.72, 9.6, 4.8, 20, 56, 84, 786]; %energy of beam from appendix
F = [25e6, 25e6, 25e6, 17e6, 23e6, 10e6, 7e6, 4e6]; % Flux from appendix
et = [5, 5, 5, 5, 5, 10, 10, 15]; %exposure time from appendix
flip = [8, 7497, 22514, 23986, 33810, 29991, 21022, 18043]; %number of flips from appendix
u = [12,12, 12, 16, 40, 56, 84 ,131];  %atomic number from appendix

sigma = zeros (8 ,1);
E = zeros (8 ,1);
for i=1:8 %to calculate for each instancec in the table in the appendix
    sigma(i) = flip(i)./(et(i).*60.*F(i)); %cross section
    E(i) = E_b(i)./u(i); %Energy
end

LET = [2.5e3, 3e3, 4.7e3, 6e3, 1.7e4, 2.7e4, 3.8e4, 2.4e4]
x= abs(lsqcurvefit(@weibull ,[3.01e-4, 1, 1, 2] , LET , sigma ));
clf reset ;
plot ( LET , sigma ,'m +');
hold on ;
xdata = linspace (0 ,60 ,100);
plot ( xdata , weibull (x , xdata ));
xlabel (' LET ( MeV cm ^2/ mg )')
ylabel ('\ sigma ( LET ) ( cm ^2)')
grid on
end

