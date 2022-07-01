% RSW theory for N layer f-plane modes
% twnh Feb '18

% Housekeeping
close all
clear
more off

%% Define problem
N = 2 ;
fcoriolis         = 1e-4 ;
gravity           = 9.81 ;
rho0              = 1025 ;
Total_H           = 4000 ;
Delta_rho         = 5 ;
layer_densities   = rho0+(Delta_rho/(N-1)).*[0:N-1]' ;
fprintf(1,' RSW theory for N=%d layer f-plane modes.\n twnh Feb ''18\n\n',N) ;

% Compute variables
layer_thicknesses = (Total_H/N).*ones(N,1) ;
reduced_gravities = gravity.*([layer_densities;0]-[0;layer_densities])./rho0 ;
reduced_gravities = reduced_gravities(1:N) ;
defrad            = sqrt(gravity*sum(layer_thicknesses))/fcoriolis ;
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

%% Create symbolic variables
syms k l omega f
gp = sym('gp',[1 N+1]) ;
H  = sym('H', [1 N]) ;
A  = sym('A',3*N) ;

%% Build A matrix
fprintf(1,' Building matrix...') ;
tic
for nn = 1:3*N
    A(nn,:)  = 0 ;            % Zero all elements
end % nn

for nn = 1:3:3*N          % Loop over rows: nn is the x-mom eqn for each layer in turn
    layer        = floor((nn-1)/3)+1 ;
    thisH        = H(layer) ;
    A(nn,nn+1)   = -f ;            % Coriolis u-mom
    A(nn+1,nn)   =  f ;            % Coriolis v-mom
    A(nn+2,nn)   = 1i*k*thisH ;    % Continuity
    A(nn+2,nn+1) = 1i*l*thisH ;    % Continuity
end % nn

toc
tic
% Add pressure gradient terms. This is the slow part for large N
for ll = 1:N                       % Loop over layer
    thisrow = ((ll-1)*3)+1 ;
    for mm = 0:ll-1                % First  sum in pressure gradient
        for kk = mm+1:N            % Second sum in pressure gradient
            thisg   = gp(mm+1) ;
            thiscol = ((kk-1)*3)+3 ;
            A(thisrow  ,thiscol) = A(thisrow  ,thiscol) + 1i*k*thisg ;
            A(thisrow+1,thiscol) = A(thisrow+1,thiscol) + 1i*l*thisg ;
            fprintf(1,'.') ;
        end % kk
    end % mm
    fprintf(1,'\n') ;
end % ll

%% Substitute the parameters
tmp = subs(  A,f,fcoriolis) ;
tmp = subs(tmp,l,0) ;
tmp = subs(tmp,k,k*kdefrad) ;
for nn = 1:N
    thisgp = reduced_gravities(nn) ;
    thisH  = layer_thicknesses(nn) ;
    tmp    = subs( tmp,gp(nn),thisgp) ;
    tmp    = subs( tmp,H( nn),thisH ) ;
end % nn
fprintf(1,'done in [%g]s\n',toc) ;
A = -tmp.*1i ;          % Accounts for -i*omega on diagonal, which has been replaced with eigenvalue analysis here rather than vanishing of the determinant (because it's quicker).

%% Find the modes
fprintf(1,' Finding the modes...') ;
tic
kks = 10.^[-2:0.1:3]' ;
omega_nondim = zeros(length(kks),N) ;
this_l       = 0 ;
for kk = 1:length(kks)
    this_k             = kks(kk) ;
    tmp                = double(subs( A,k,this_k)) ;
    omegas             = sortrows(real(eig(tmp)))./fcoriolis ;
    omega_nondim(kk,:) = omegas(end-N+1:end) ;
end %kk

% Infer deformation radii
def_rads = sqrt(mean((omega_nondim.^2 - 1)./kks.^2))' ;
fprintf(1,'done in [%g]s\n',toc) ;
fprintf(1,' Deformation radii = [%g-%g].\n\n',min(def_rads),max(def_rads)) ;

figure(1)
for ff = 1:N
   loglog(kks,omega_nondim(:,ff),'-') ;
   hold on
   fitted_omega = sqrt(1 + def_rads(ff)^2.*(kks.^2)) ;
   loglog(kks,fitted_omega,'o') ;
end % ff
grid on

figure(2)
semilogy(flipud(def_rads),'o-')
ylabel('deformation radius')
xlabel('mode') ;
grid on
set(gca,'XScale','log')