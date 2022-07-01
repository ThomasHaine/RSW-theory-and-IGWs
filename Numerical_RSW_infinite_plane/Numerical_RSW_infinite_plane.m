% Script to solve the RSW adjustment problem numerically
% twnh Feb '21

%% Housekeeping
clear
close all
more off
clc
fprintf(1,' RSW_adjustment_numerical.m\n Script to solve the 1-layer linear RSW equations numerically from an arbitrary initial condition.\n twnh Feb ''21\n\n')

%% Parameters
fprintf(1,' Setup...') ;
% Physical parameters
f = 1 ;         % Coriolis parameter
g = 1 ;         % Gravitational acceleration
H = 1 ;         % Resting layer depth
L = 16 ;        % Length of domain
M = 16 ;        % Width of domain

% Numerical parameters
Nx = 512 ;      % Number of discretization points in x
Ny = 512 ;      % Number of discretization points in y
Nt = 256 ;      % Number of discretization points in t
Np = 32 ;       % Number of Lagrangian particles
Tf = (L/3)/sqrt(g*H) ;
times = linspace(0,Tf,Nt) ;

% Compute
NN   = Nx*Ny ;
dx   = L/Nx ;
dy   = M/Ny ;
xvec = -L/2:dx:L/2-dx ;
yvec = -M/2:dy:M/2-dy ;
[xgrid,ygrid] = ndgrid(xvec,yvec) ;
x    = xgrid(:) ;
y    = ygrid(:) ;
kvec = (2*pi/L)*(-Nx/2:Nx/2-1) ;
lvec = (2*pi/M)*(-Ny/2:Ny/2-1) ;
kvec = fftshift(kvec) ;
lvec = fftshift(lvec) ;
[kgrid,lgrid] = ndgrid(kvec,lvec) ;
k    = kgrid(:) ;
l    = lgrid(:) ;
fprintf(1,'done.\n Grid: Nx, Ny, Nt = [%d, %d, %d].\n',Nx,Ny,Nt) ;

%% Initial condition
fprintf(1,' Initial condition...') ;
%u0      = 0*xgrid ;
u0       = (exp(-xgrid.^2/((L/32)^2))).*(exp(-ygrid.^2/((M/32)^2))) ;
v0       = 0*xgrid ;
%eta0     = 0*xgrid ;
eta0    = 0.5*(exp(-xgrid.^2/((L/32)^2))).*(exp(-ygrid.^2/((M/32)^2))) ;
u0hat   = fft2(u0) ;
v0hat   = fft2(v0) ;
eta0hat = fft2(eta0) ;
y0hat   = [u0hat(:); v0hat(:); eta0hat(:)] ;
fprintf(1,'done.\n') ;

%% Simulate in Fourier frequency domain
fprintf(1,' Build propagator matrix...') ;
tic
A = make_A_matrix(k,l,g,H,f) ;
fprintf(1,'done\n') ;

% expm solution: This works, but scales worse with grid size.
if(0)
    fprintf(1,' Solve problem with expm...') ;
    tic
    Delta_t = times(2) - times(1) ;
    propagator = expm(A*Delta_t) ;
    %[V,Omega] = eigs(A,32) ;
    %propagator = V*exp(Omega*Delta_t)/V ;
    
    u = zeros(Nt,Nx,Ny) ;
    v = zeros(Nt,Nx,Ny) ;
    eta = zeros(Nt,Nx,Ny) ;
    
    for tt = 1:numel(times)
        if(tt == 1)
            this_yhat = y0hat ;
        else
            this_yhat   = propagator*this_yhat ;
        end % if
        this_uhat   = reshape(this_yhat(     1:  NN),  Nx,Ny) ;
        this_vhat   = reshape(this_yhat(  NN+1:2*NN),  Nx,Ny) ;
        this_etahat = reshape(this_yhat(2*NN+1:3*NN),  Nx,Ny) ;
        u(tt,:,:)   = ifft2(this_uhat) ;        % ifft2 to return to spatial domain
        v(tt,:,:)   = ifft2(this_vhat) ;
        eta(tt,:,:) = ifft2(this_etahat) ;
    end
    fprintf(1,'done in [%6.3f]s.\n',toc) ;
end % if

% Numerical ode solver:
fprintf(1,' Solve problem with ode45...') ;
tic
u    = zeros(Nt,Nx,Ny) ;
v    = zeros(Nt,Nx,Ny) ;
eta  = zeros(Nt,Nx,Ny) ;
zeta = zeros(Nt,Nx,Ny) ;
dive = zeros(Nt,Nx,Ny) ;
[~,yhat] = ode45(@(t,yhat)rhs_RSW(t,yhat,A),times,y0hat) ;

for tt = 1:numel(times)
    this_uhat    = reshape(yhat(tt,     1:  NN),  Nx,Ny) ;
    this_vhat    = reshape(yhat(tt,  NN+1:2*NN),  Nx,Ny) ;
    this_etahat  = reshape(yhat(tt,2*NN+1:3*NN),  Nx,Ny) ;
    u(tt,:,:)    = real(ifft2(this_uhat)) ;        % ifft2 to return to spatial domain
    v(tt,:,:)    = real(ifft2(this_vhat)) ;
    eta(tt,:,:)  = real(ifft2(this_etahat)) ;
    this_zetahat = -1i.*(kgrid.*this_vhat - lgrid.*this_uhat) ;   % Compute vorticity
    this_divehat = -1i.*(kgrid.*this_uhat + lgrid.*this_vhat) ;   % Compute divergence
    zeta(tt,:,:) = real(ifft2(this_zetahat)) ;
    dive(tt,:,:) = real(ifft2(this_divehat)) ;
end
PV = (f+zeta)./(H+eta) ;
fprintf(1,'done in [%6.3f]s.\n',toc) ;

%% Compute Lagrangian trajectories
fprintf(1,' Solve Lagrangian trajectories...') ;
part_x_0 = (L/16).*randn((Np-1)*2,1) ;      % Start them near the origin
part_x_0 = [0; part_x_0(1:Np); 0; part_x_0(Np+1:end)] ;     % Add one at the origin
[tmat,xmat,ymat] = ndgrid(times,xvec,yvec) ;
uInterp = griddedInterpolant(tmat,xmat,ymat,u) ;
vInterp = griddedInterpolant(tmat,xmat,ymat,v) ;
[~,part_x] = ode45(@(t,part_x)rhs_particle(t,part_x,uInterp,vInterp),times,part_x_0) ;
fprintf(1,'done in [%6.3f]s.\n',toc) ;

%% FIGURES (PRODUCTION)

vidfiles{1} = VideoWriter('eta.mp4','MPEG-4');
open(vidfiles{1}) ;
vidfiles{2} = VideoWriter('PV.mp4','MPEG-4');
open(vidfiles{2}) ;
vidfiles{3} = VideoWriter('zeta.mp4','MPEG-4');
open(vidfiles{3}) ;
vidfiles{4} = VideoWriter('divergence.mp4','MPEG-4');
open(vidfiles{4}) ;
for tt = 1:numel(times)
    plot_pcolor(1,xgrid,ygrid,tt,eta,'$\eta$',L,M,part_x,vidfiles) ;
    plot_pcolor(2,xgrid,ygrid,tt,PV,'$\frac{f+\zeta}{h}$',L,M,part_x,vidfiles) ;
    plot_pcolor(3,xgrid,ygrid,tt,zeta,'$\zeta$',L,M,part_x,vidfiles) ;
    plot_pcolor(4,xgrid,ygrid,tt,dive,'$\delta$',L,M,part_x,vidfiles) ;
end % tt
close(vidfiles{1})
close(vidfiles{2})
close(vidfiles{3})
close(vidfiles{4})

%% Local functions

function dyhatdt = rhs_RSW(~,yhat,A)
dyhatdt = A*yhat ;
end

function d_part_xdt = rhs_particle(t,part_x,uInterp,vInterp)
Np = numel(part_x)/2 ;
tvec = t.*ones(Np,1) ;
u = uInterp(tvec,part_x(1:Np),part_x(Np+1:end)) ;
v = vInterp(tvec,part_x(1:Np),part_x(Np+1:end)) ;
d_part_xdt = [u; v] ;
end

function A = make_A_matrix(k,l,g,H,f)

NN    = numel(k) ;
fvec  = f.*ones(NN,1) ;
gkvec = 1i*g.*k ;
glvec = 1i*g.*l ;
Hkvec = 1i*H.*k ;
Hlvec = 1i*H.*l ;
A1    = spdiags([ fvec,gkvec],[NN, 2*NN],NN,3*NN) ;
A2    = spdiags([-fvec,glvec],[ 0, 2*NN],NN,3*NN) ;
A3    = spdiags([Hkvec,Hlvec],[ 0,   NN],NN,3*NN) ;
A     = [A1; A2; A3] ;

end

function plot_pcolor(ff,xgrid,ygrid,tt,fld,tit_str,L,M,part_x,vidfiles)

% Prepare
figure(ff)
%set(gcf,'Renderer', 'painters', 'Position', [10 10 900 900]) ;
set(gcf,'Position', [10 10 800 800]) ;
%set(gca,'FontSize',48) ;
hold off

% Field
pcolor(xgrid,ygrid,squeeze(fld(tt,:,:))) ;
shading flat
set(gca,'XLim',0.5.*[min(xgrid(:)) max(xgrid(:))],'YLim',0.5.*[min(ygrid(:)) max(ygrid(:))]) ;

colorbar('southoutside') ;
clims = [min(fld(:)) max(fld(:))] ;
caxis(clims) ;
title(tit_str,'interpreter','latex')
daspect([L M 1]) ;
xlabel('$\frac{x f}{\sqrt{g H}}$','interpreter','latex') ;
ylabel('$\frac{y f}{\sqrt{g H}}$','interpreter','latex') ;

% Particles
NP = size(part_x,2)/2 ;
hold on
plot(part_x(tt,1:NP),part_x(tt,NP+1:end),'.') ;

% Tidy up
set(gca,'Box','on') ;
c_map = flipud(diverging_map(linspace(0,1,64),[0.230, 0.299, 0.754],[0.706, 0.016, 0.150])) ;
colormap(c_map) ;
drawnow

% Save frame
F = getframe(gcf);
writeVideo(vidfiles{ff},F);

end