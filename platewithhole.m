clear all;

% Define all constants

% Misc

rt2 = sqrt(2);
tol = 1.e-6;
maxiter = 100;

% Material properties - Aluminum (N/mm^2)

lambda = 60.49e3; % Lame constant
mu = 25.93e3; % Shear modulus

props = [lambda mu]; % Array of material properties
nint = 0; % Number of internal state variables

% Geometry and meshing

width = 1.e2; % Half the width of the plate (mm)
height = 1.e2; % Half the height of the plate (mm)
vvf = 0.01; % Volume fraction of the hole
Nr = 20; % Number of elements in the radial direction
Nt = 40; % Number of elements in the tangential direction

% Loading parameters

rdisp = 1.0; % Displacement of the right edge (mm). Left edge is fixed.
ninc = 20; % Number of increments
ddisp = rdisp/ninc; % Incremental displacement

% Call the mesh generator

[nels nnodes nodexy con bound left right bottom top void nvoidE nvoidN] = holemesh(width,height,vvf,Nr,Nt);

nleft = [left(1,:) left(2,length(left(1,:)))];
nright = [right(1,:) right(2,length(right(1,:)))];

%triangulation_plot_eps('original_mesh.eps', nnodes, nodexy, 2*nels, [[con(1,:);con(2,:);con(4,:)] [con(2,:);con(3,:);con(4,:)]]);

% Initialization

gauss = 1/sqrt(3);
xilist = gauss*[-1 1 1 -1];
etlist = gauss*[-1 -1 1 1];

u = zeros(2*nnodes,1); % Nodal displacements
sig = zeros(4*nels,4); % Stresses at the integration points
epsp = zeros(4*nels,4); % Plastic strains at the integration points
intvar = zeros(4*nels,nint); % Internal variables at the integration points

rdisp = 0; % Displacement of the right edge
rfrc = 0; % Force on the right edge
disphist = [rdisp];
frchist = [rfrc];

% Incremental loading loop

for inc=1:ninc
  inc

  % Save current material state

  uprev = u;
  sigprev = sig;
  epspprev = epsp;
  intvprev = intvar;

  % Setup the RHS vector

  Q = zeros(2*nnodes,1);
  Q(nright) = ddisp; % Dirichlet boundary conditions

  % Solve for the incremental nodal displacements by iteration

  du = zeros(2*nnodes,1);

  for iter=1:maxiter

    % Set up the tangent stiffness and R vector

    Ktan = zeros(2*nnodes);
    R = zeros(2*nnodes,1);

    for e=1:nels
      n = con(:,e);
      xe = nodexy(1,n);
      ye = nodexy(2,n);

      due = [du(n)'; du(nnodes+n)'];

      K11 = zeros(4); K22 = zeros(4); K12 = zeros(4);
      Re = zeros(4,2);

      for i=1:4
        xi = xilist(i);
        et = etlist(i);

        M = 1/4*[-(1-et) -(1-xi); (1-et) -(1+xi); (1+et) (1+xi); -(1+et) (1-xi)];

        J = [xe*M(:,1) xe*M(:,2); ye*M(:,1) ye*M(:,2)];
        dJ = det(J);
        Jinv = inv(J);

        N = M*Jinv;

        graddu = due*N;
        deps = 1/2*(graddu+graddu');
        deps = [deps(1,1); deps(2,2); 0; 0; 0; rt2*deps(1,2)]; % Convert to Voigt form

        sigvec = sigprev(4*(e-1)+i,:)';
        sigma = [sigvec(1); sigvec(2); sigvec(3); 0; 0; rt2*sigvec(4)]; % Voigt form
        qvec = intvprev(4*(e-1)+i,1:nint)';

        [dsig,depsp,dqvec,ddsdde] = stressinc(sigma,qvec,deps,props);

        stif11 = [ddsdde(1,1) ddsdde(1,6)/rt2; ddsdde(1,6)/rt2 ddsdde(6,6)/2];
        stif22 = [ddsdde(6,6)/2 ddsdde(2,6)/rt2; ddsdde(2,6)/rt2 ddsdde(2,2)];
        stif12 = [ddsdde(1,6)/rt2 ddsdde(1,2); ddsdde(6,6)/2 ddsdde(2,6)/rt2];

        K11 = K11 + N*stif11*N' * dJ;
        K22 = K22 + N*stif22*N' * dJ;
        K12 = K12 + N*stif12*N' * dJ;

        Re = Re + N*[dsig(1) dsig(6)/rt2; dsig(6)/rt2 dsig(2)] * dJ;

        sig(4*(e-1)+i,:) = sigprev(4*(e-1)+i,:)+[dsig(1) dsig(2) dsig(3) dsig(6)/rt2];
        epsp(4*(e-1)+i,:) = epspprev(4*(e-1)+i,:)+[depsp(1) depsp(2) depsp(3) depsp(6)/rt2];
        if (nint > 0)
          intvar(4*(e-1)+i,1:nint) = intvprev(4*(e-1)+i,1:nint)+dqvec';
        end;
      end;

      Ke = [K11 K12; K12' K22];

      % Assembly
      for i=1:4
      for j=1:4
        Ktan(n(i),n(j)) = Ktan(n(i),n(j)) + Ke(i,j);
        Ktan(n(i),nnodes+n(j)) = Ktan(n(i),nnodes+n(j)) + Ke(i,4+j);
        Ktan(nnodes+n(i),n(j)) = Ktan(nnodes+n(i),n(j)) + Ke(4+i,j);
        Ktan(nnodes+n(i),nnodes+n(j)) = Ktan(nnodes+n(i),nnodes+n(j)) + Ke(4+i,4+j);
      end;
      end;

      for i=1:4
        R(n(i)) = R(n(i)) + Re(i,1);
        R(nnodes+n(i)) = R(nnodes+n(i)) + Re(i,2);
      end;
    end;

    Ksave = Ktan; % Save the secent stiffness for force computations

    % Dirichlet boundary conditions

    for i=1:length(nleft)
      n = nleft(i);
      Ktan(n,:) = 0;
      Ktan(n,n) = 1;
      R(n) = Q(n);

      if (i==length(nleft)) % Fully fix bottom-left corner node
        Ktan(nnodes+n,:) = 0;
        Ktan(nnodes+n,nnodes+n) = 1;
        R(nnodes+n) = Q(nnodes+n);
      end
    end;

    for i=1:length(nright)
      n = nright(i);
      Ktan(n,:) = 0;
      Ktan(n,n) = 1;
      R(n) = Q(n);
    end;

    % Solve for the incremental nodal displacements

    if (iter==1)
      ddu = Ktan\Q;
    else
      ddu = Ktan\(Q-R);
    end
    du = du+ddu;

    if (max(abs(ddu)) <= tol) % Convergence check
      break;
    end;

  end; % Equilibrium iterations

  if (iter == maxiter)
    u = uprev;
    sig = sigprev;
    epsp = epspprev;
    intvar = intvprev;

    disp('Convergence failure');
    break;
  end

  u = u+du;

  % Calculate the force on the right edge

  dfrc = Ksave*du;
  drfrc = sum(dfrc(nright));
  rfrc = rfrc+drfrc;
  frchist = [frchist rfrc];

  rdisp = rdisp+ddisp;
  disphist = [disphist rdisp];

end; % Incremental loading loop


% Plot the force-displacement graph

plot(disphist,frchist,'-o');
title('Force-displacement curve')

% Extrapolate stresses and plastic strains to the nodes

sigma = zeros(4,nnodes);
eplas = zeros(4,nnodes);
ncount = zeros(nnodes,1);

for e=1:nels
  n = con(:,e);

  sip = sig(4*e-3:4*e,:);
  epip = epsp(4*e-3:4*e,:);

  xtrapl = [1+sqrt(3)/2, -1/2, 1-sqrt(3)/2, -1/2,
            -1/2, 1+sqrt(3)/2, -1/2, 1-sqrt(3)/2,
            1-sqrt(3)/2, -1/2, 1+sqrt(3)/2, -1/2,
            -1/2, 1-sqrt(3)/2, -1/2, 1+sqrt(3)/2];

  snod = xtrapl*sip;
  epnod = xtrapl*epip;

  sigma(:,n) = sigma(:,n) + snod';
  eplas(:,n) = eplas(:,n) + epnod';
  ncount(n) = ncount(n)+1;
end;

% Compute node averaged stresses and plastic strains

nodexydef = zeros(2,nnodes);
for i=1:nnodes
  sigma(:,i) = sigma(:,i)/ncount(i);
  eplas(:,i) = eplas(:,i)/ncount(i);

  nodexydef(1,i) = nodexy(1,i) + u(i);
  nodexydef(2,i) = nodexy(2,i) + u(nnodes+i);
end;


%triangulation_plot_eps('deformed_mesh.eps', nnodes, nodexydef, 2*nels, [[con(1,:);con(2,:);con(4,:)] [con(2,:);con(3,:);con(4,:)]]);


% Plot contour maps of the stresses

x1 = reshape(nodexy(1,con'),nels,4)';
y1 = reshape(nodexy(2,con'),nels,4)';

s11 = reshape(sigma(1,con'),nels,4)';
s22 = reshape(sigma(2,con'),nels,4)';
s33 = reshape(sigma(3,con'),nels,4)';
s12 = reshape(sigma(4,con'),nels,4)';

ep11 = reshape(eplas(1,con'),nels,4)';
ep22 = reshape(eplas(2,con'),nels,4)';
ep33 = reshape(eplas(3,con'),nels,4)';
ep12 = reshape(eplas(4,con'),nels,4)';


%Mean Stress and plastic strain
sigm= (s11 + s22 + s33)/3;
epm= (ep11 + ep22 + ep33)/3;

%Equivalent stress and plastic strain
seq= sqrt((s11-s22).^2 + (s22-s33).^2 + (s33-s11).^2 + 6*((s12).^2)/2);
epeq= sqrt((2/3)*(ep11.^2 + ep22.^2 + ep33.^2 + 2*(ep12.^2 )));
%%
%Plotting Hydrostatic
figure;
patch(x1,y1,sigm);
axis equal auto
shading interp
colorbar;
title('Hydro-static Stress using D-P model')

figure;
patch(x1,y1,epm);
axis equal auto
shading interp
colorbar;
title('Volumetric plastic strain using D-P model')

%Plotting equivalent
figure;
patch(x1,y1,seq);
axis equal auto
shading interp
colorbar;
title('Equivalent stress using D-P model')

figure;
patch(x1,y1,epeq);
axis equal auto
shading interp
colorbar;
title('Equivalent plastic strain using D-P model')

% figure;
% patch(x1,y1,s11);
% axis equal auto
% shading interp
% title('sigma_{11} using D-P model')
% 
% figure;
% patch(x1,y1,s22);
% axis equal auto
% shading interp
% title('sigma_{22} using D-P model')
% 
% 
% figure;
% patch(x1,y1,s12);
% axis equal auto
% shading interp
% title('sigma_{12} using D-P model')

% figure;
% patch(x1,y1,ep11);
% axis equal auto
% shading interp
% 
% figure;
% patch(x1,y1,ep22);
% axis equal auto
% shading interp
% 
% figure;
% patch(x1,y1,ep12);
% axis equal auto
% shading interp