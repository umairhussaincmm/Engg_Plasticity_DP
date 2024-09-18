%Drucker-Prager plasticity using Cutting Plane Algorithm
%initialisation
a=-1/2;
b=0;
n=0.1;
nu= 0.333;

E= 200;
phi= 35;
theta= 28;
sigo= E/500;
mu= E/(2*(1+nu));
lambda= (E*nu)/((1+nu)*(1-2*nu));

alpha=(2*sin(phi))/(sqrt(3)*(3-sin(phi)));
beta=(2*sin(theta))/(sqrt(3)*(3-sin(theta)));

steps=50;
eo= linspace(0,.01,steps);


iteration=100;

tol = 1.e-12;

h=0.00001;

ee= [1 a b 0 0 0];
e= zeros(6,steps);
sig= zeros(6,steps);
dev= zeros(6,steps);
ep= zeros(6,steps);
s= zeros(1,steps);
sigy= zeros(1,steps);
sigy(1)= sigo*(1+500*s(1)).^n;
I2 = [1; 1; 1; 0; 0; 0]; % 2nd order identity tensor in Voigt form
stif = 2*mu*eye(6)+lambda*(I2*I2'); % Elastic stiffness tensor in Voigt form

for i=1:steps
    e(:,i)= ee*eo(i);
end

%disp(stif);
%%
%increment loop
for i= 1:steps-1
    i
    de= e(:,i+1) - e(:,i);
    dlam=0;
    %disp(de);
    dsigt= stif*de;
    sigt= sig(:,i) + dsigt;
    %trace_sigt= sigt(1) + sigt(2) + sigt(3);
    %devt= sigt - trace_sigt*I2;
    %dotprod= dot(devt,devt);
    %deveq= sqrt((3/2)*dotprod);
    f= DPyieldf(sigt, alpha, sigo);
    
    if (f) <= tol
        sig(:,i+1)= sigt;
        ep(:,i+1)= ep(:,i);
    else
        dlamk= 0;
        depk= [0; 0; 0; 0; 0; 0];
        dsigk= dsigt;
        sigk= sig(:,i) + dsigk;
        
        
        for j=1:iteration
            j
            f= DPyieldf(sigk, alpha, sigo); %yield function value at current state
            if abs(f)<=tol %checking for convergence
                sig(:,i+1)= sigk ;%Updating values after convergence
                dep= depk;
                ep(:,i+1)= ep(:,i) + dep;
                break;
            end
            
            N= DPflowf(sigk, beta);%explicit estimation
            
            %Updating values now
            
            %Newton-Raphson- for updating dlambda
            x0=0;
            for k=1:iteration
                temp= x0;%storing previous iteration value
                x11= x0+h;%taking a point ahead current point
                x1= x0-h;%taking a point before current point
                
                dsigk= dsigt-(dlamk+x0)*stif*N;
                sigk= sig(:,i) + dsigk;
                f0= DPyieldf(sigk, alpha, sigo);%yield function value at x0
                
                dsigk= dsigt-(dlamk+x1)*stif*N;
                sigk= sig(:,i) + dsigk;
                f1= DPyieldf(sigk, alpha, sigo);%yield function value at x1
                
                dsigk= dsigt-(dlamk+x11)*stif*N;
                sigk= sig(:,i) + dsigk;
                f11= DPyieldf(sigk, alpha, sigo);%yield function value at x11
                
                df0= (f11-f1)/(2*h);%Derivative estimation numerically
                
                x0= x0- f0/df0;%Newton-Raphson formula
                
                err= abs(temp-x0);
                if err<5*10^-(8) %checking for convergence in Newton-raphson
                    break;
                end
                if k== iteration
                    disp('Newton-Raphson failure!');
                    return;
                end
            end
            
            ddlam=x0;
            
            %Updating
            dlamk= dlamk + ddlam;
            %N= (3/(2*deveqk))*devk + (beta/3)*I2;
            depk= dlamk*N;
            dsigk= stif*(de - depk);
            sigk= sig(:,i) + dsigk;
            
            
            
            if j==iteration
                disp('Convergence Failure');
                return;
            end
        end
            
    
    end
    
end
%disp(sig);
%disp(ep);

%%
sig11= zeros(1,steps);
sig22= zeros(1,steps);
sig33= zeros(1,steps);
sigm= zeros(1,steps);
sigeq= zeros(1,steps);

for i=1:steps
    sig11(1,i)= sig(1,i);
    sig22(1,i)= sig(2,i);
    sig33(1,i)= sig(3,i);
    sigm(1,i)= (sig(1,i) + sig(2,i) + sig(3,i))/3;
    dev(:,i)= sig(:,i) - sigm(1,i)*I2;
    dotpro= dot(dev(:,i),dev(:,i));
    sigeq(1,i)= sqrt((3/2)*dotpro);
end
figure;
plot(sigm,sigeq)
title('V-M equivalent stress as a function of Hydro-static stress(D-P)')

% figure;
% plot(eo,sigeq)
% title('V-M equivalent stress as a function of total strain(D-P)')

%%

figure;
plot(eo,sigm)
title('Hydro-static stress as a function of total strain')

hold on

ep11= zeros(1,steps);
ep22= zeros(1,steps);
ep33= zeros(1,steps);
epm= zeros(1,steps);

for i=1:steps
    ep11(1,i)= ep(1,i);
    ep22(1,i)= ep(2,i);
    ep33(1,i)= ep(3,i);
    epm(1,i)= (ep(1,i) + ep(2,i) + ep(3,i))/3;
end

figure;
plot(eo,epm)
title('Volumetric plastic strain as a function of total strain(D-P)')
% 
% figure;
% plot(epm,sigm)
% title('Mean stress as a function of Mean plastic strain')