% Rossby wave dispersion relation
% twnh Feb 18

% Housekeeping
clear
close all

% Parameters
kstar     = [-10:0.02:10] ;
beta      = 1 ;
ellstar   = [0:0.5:4] ;

% Compute
for ll = 1 :length(ellstar)
    Kstar2 = ellstar(ll).^2 + kstar.^2 ;
omegastar = -kstar./(Kstar2 + 1) ;
hold on
plot(kstar,omegastar,'linewidth',2)
end % ll

grid on
set(gca,'YLim',[0 0.5]) ;
xlabel('$k$','interpreter','latex') ;
ylabel('$\omega$','interpreter','latex') ;