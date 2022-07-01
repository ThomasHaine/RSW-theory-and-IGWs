% RSW theory for N layer f-plane modes
% No symbolic variables in this version
% twnh Feb '18

% Housekeeping
close all
clear
more off

%% Define problem
N           = 10 ;
gravity     = 9.81 ;
rho0        = 1025 ;
latitude    = 15.0 ;
Omega       = 7.292e-5 ;
fprintf(1,' RSW theory for N=%d layer f-plane modes.\n twnh Feb ''18\n\n',N) ;
fprintf(1,' This version is purely numerical (no symbolic calculations, unlike earlier versions), which is much faster.\n\n') ;


%% Define stratification.
if(1)
    % Eq. Pac.
    rho_profile = 1e3.*[
        0.0223583           0
        0.0223947  -0.0100000
        0.0224207  -0.0200000
        0.0224462  -0.0300000
        0.0225148  -0.0500000
        0.0226095  -0.0750000
        0.0228001  -0.1000000
        0.0232596  -0.1250000
        0.0240957  -0.1500000
        0.0257601  -0.2000000
        0.0263685  -0.2500000
        0.0265599  -0.3000000
        0.0267711  -0.4000000
        0.0269425  -0.5000000
        0.0270778  -0.6000000
        0.0271781  -0.7000000
        0.0272535  -0.8000000
        0.0273197  -0.9000000
        0.0273782  -1.0000000
        0.0274229  -1.1000000
        0.0274652  -1.2000000
        0.0275006  -1.3000000
        0.0275369  -1.4000000
        0.0275670  -1.5000000
        0.0276246  -1.7500000
        0.0276595  -2.0000000
        0.0277113  -2.5000000
        0.0277332  -3.0000000
        0.0277509  -3.5000000
        0.0277665  -4.0000000
        0.0277777  -4.5000000
        0.0277888  -5.0000000 ] ;
    rho_profile = [rho_profile(:,2),rho_profile(:,1)] ;
    total_H = 5000.0 ;
elseif(0)
    % Doctored tropical Atlantic data.
    rho_profile = [
        25.568  0.0   ;
        25.568  3.000 ;
        26.430 125.000 ;
        26.926 350.000 ;
        27.296 750.000 ;
        27.605 1200.000 ;
        27.713 1400.000 ;
        27.763 1600.000 ;
        27.798 1800.000 ;
        27.844 2400.000 ;
        27.862 3000.000 ;
        27.877 3201.000 ;
        27.878 3601.000 ;
        27.888 3901.000 ;
        27.890 4100.000 ;
        27.894 4400.000 ;
        27.896 4800.000 ;
        27.897 5000.0
        ] ;
    rho_profile = [-rho_profile(:,2),rho_profile(:,1)] ;
    total_H = 5000.0 ;
    
elseif(0)
    % Uniform stratification. Matches Emery & Thompson example, p344. Useful for testing code.
    rho_profile = [
        24.000  0.0   ;
        25.045  2500 ;
        ] ;
    total_H = 2500.0 ;
    rho_profile = [-rho_profile(:,2),rho_profile(:,1)] ;
end % if

% Set grid spacing and interpolate density profile to this grid.
delta =  total_H/N ;   % m
grid2 = -[delta/2:delta:total_H-delta/2]' ;
layer_densities = 1000+interp1(rho_profile(:,1),rho_profile(:,2),grid2) ;

%% Compute variables
fcoriolis   = 2*Omega*sind(latitude) ;
layer_thicknesses = (total_H/N).*ones(N,1) ;
reduced_gravities = gravity.*([layer_densities;0]-[0;layer_densities])./rho0 ;
reduced_gravities = reduced_gravities(1:N) ;
defrad            = sqrt(gravity*sum(layer_thicknesses))/fcoriolis ;  % Barotropic deformation radius.
kdefrad           = 1/defrad ;
% NB The first reduced gravity is the actual gravity!
%    This convention is not the same as usual in texts
%    (because I can't define "gp0" as the actual gravity).

%% Plot density profile
figure(3)
for ll = 1:N
    if ll > 1
        this_start_z = sum(layer_thicknesses(1:ll-1)) ;
    else
        this_start_z = 0 ;
    end
    this_end_z   = this_start_z + layer_thicknesses(ll) ;
    plot([layer_densities(ll) layer_densities(ll)],-[this_start_z this_end_z],'k-') ;
    hold on
    if ll < N
        plot([layer_densities(ll) layer_densities(ll+1)],-[this_end_z this_end_z],'k-') ;
    end % if
end % ll
grid on
xlabel('Density (kg/m^3)') ;
ylabel('Height (m)') ;

%% Build A matrix
fprintf(1,' Building matrix...') ;
tic
A = zeros(3*N) ;

for nn = 1:3:3*N          % Loop over rows: nn is the x-mom eqn for each layer in turn
    layer        = floor((nn-1)/3)+1 ;
    thisH        = layer_thicknesses(layer) ;
    A(nn,nn+1)   = -fcoriolis ;            % Coriolis u-mom
    A(nn+1,nn)   =  fcoriolis ;            % Coriolis v-mom
    A(nn+2,nn)   = 1i*thisH ;              % Continuity: multiply by k in loop below
    A(nn+2,nn+1) = 1i*thisH ;              % Continuity: multiply by l in loop below
end % nn

% Add pressure gradient terms. This is the slow part for large N
for ll = 1:N                       % Loop over layer
    thisrow = ((ll-1)*3)+1 ;
    for mm = 0:ll-1                % First  sum in pressure gradient
        for kk = mm+1:N            % Second sum in pressure gradient
            thisg   = reduced_gravities(mm+1) ;
            thiscol = ((kk-1)*3)+3 ;
            A(thisrow  ,thiscol) = A(thisrow  ,thiscol) + 1i*thisg ; % multiply by k in loop below
            A(thisrow+1,thiscol) = A(thisrow+1,thiscol) + 1i*thisg ; % multiply by l in loop below
        end % kk
    end % mm
end % ll
A = -A.*1i ;          % Accounts for -i*omega on diagonal, which has been replaced with eigenvalue analysis here rather than vanishing of the determinant (because it's quicker).
fprintf(1,'done in [%g]s\n',toc) ;

%% Find the modes
fprintf(1,' Finding the modes...') ;
tic
kks = 10.^[-2:0.1:5]' ;
omega_nondim = zeros(length(kks),N) ;
this_l       = 0 ;
for kk = 1:length(kks)
    this_k             = kks(kk)*kdefrad ;
    this_A             = A ;
    
    % Insert k and l into matrix for this value of k.
    for nn = 1:3:3*N          % Loop over rows: nn is the x-mom eqn for each layer in turn
        this_A(nn+2,nn)      = this_A(nn+2,nn)  *this_k ;          % Continuity
        this_A(nn+2,nn+1)    = this_A(nn+2,nn+1)*this_l;           % Continuity
        this_A(nn  ,3:3:3*N) = this_A(nn  ,3:3:3*N)*this_k ;    % Pressure gradient
        this_A(nn+1,3:3:3*N) = this_A(nn+1,3:3:3*N)*this_l ;    % Pressure gradient
    end % nn
    
    [modes,omegas]     = eig(this_A) ;
    [omegas,inds]      = sortrows(real(diag(omegas))) ;
    omega_nondim(kk,:) = omegas(end-N+1:end)./fcoriolis ;
end %kk

% Compute bottom pressure and surface elevation modes.
h_modes   = modes(3:3:3*N,inds(1:N)) ;
p_mode    = zeros(N,1) ;
surf_mode = zeros(N,1) ;
for nn = 1:N
    p_mode(nn)     = abs(sum(h_modes(:,nn).*layer_densities)) ;   % Hydrostatic pressure anomaly at the bottom
    surf_mode(nn)  = abs(sum(h_modes(:,nn))) ;                    % Surface deflection.
end % nn
p_mode    = flipud(p_mode) ;
surf_mode = flipud(surf_mode) ;

% Infer deformation radii
def_rads = sqrt(mean((omega_nondim.^2 - 1)./kks.^2))' ;
fprintf(1,'done in [%g]s\n',toc) ;
fprintf(1,' Non-dimensional deformation radii = [%g-%g].\n\n',min(def_rads),max(def_rads)) ;
fprintf(1,' Dimensional deformation radii:\n') ;
for nn = 1:min(16,N)
    fprintf(1,' %d : %g km.\n',nn,def_rads(end+1-nn)*defrad/1000) ;
end %nn


figure(1)
for ff = 1:N
   loglog(kks,omega_nondim(:,ff),'-') ;
   hold on
   fitted_omega = sqrt(1 + def_rads(ff)^2.*(kks.^2)) ;
   loglog(kks,fitted_omega,'o') ;
end % ff
xlabel('Frequency (non-dimensional)') ;
ylabel('Wavenumber (non-dimensional)') ;
grid on

figure(2)
subplot(2,1,1)
semilogy(flipud(def_rads),'o-')
ylabel('deformation radius')
xlabel('mode') ;
grid on
set(gca,'XScale','log')

subplot(2,1,2)
semilogy(flipud(p_mode),'o-')
xlabel('mode') ;
hold on
semilogy(flipud(surf_mode),'o-')
grid on
legend('bottom pressure','surface height') ;
set(gca,'XScale','log')