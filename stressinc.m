function [dsig,depsp,dqvec,ddsdde] = stressinc(sigma,qvec,deps,props)

% Input
% -----
% sigma - Current stress in Voigt form
% qvec - Vector of internal variables (current values)
% props - Array of material properties
%
% Output
% ------
% dsig - Stress increment in Voigt form
% depsp - Plastic strain increment in Voigt form
% dqvec - Increments to qvec
% ddsdde - Consistent tangent modulus in Voigt form
%
% IMPORTANT - If deps=0, ddsdde must be set to the elastic stiffness tensor in Voigt form


  tol = 1.e-12;

  lambda = props(1);
  mu = props(2);
  
  E= mu*(3*lambda + 2*mu)/(mu + lambda);
  phi= 35;
  theta= 28;
  

  I2 = [1; 1; 1; 0; 0; 0]; % 2nd order identity tensor in Voigt form
  stif = 2*mu*eye(6)+lambda*I2*I2'; % Elastic stiffness tensor in Voigt form

  nint = length(qvec);
  dqvec = zeros(nint,1);
  depsp = zeros(6,1);
  ddsdde = stif;
  
  %D-P parameters
  alpha=(2*sin(phi))/(sqrt(3)*(3-sin(phi)));
  beta=(2*sin(theta))/(sqrt(3)*(3-sin(theta)));
  sigmay= E/500;
  
  h=0.001;%for newton-raphson
  
  iteration=100;
  
  dlam=0;%elastic predictor
  dsigt= stif*deps;%trial stress increment
  sigmaft= sigma + dsigt;%trial final stress
  f= DPyieldf(sigmaft, alpha, sigmay);%Calling yield function
  
  %Checking for Yielding
  if (f) <= tol
      dsig= dsigt;
      depsp= [0; 0; 0; 0; 0; 0];
      
  else
      %Initialising the iteration
      dlamk= 0;
      depk= [0; 0; 0; 0; 0; 0];
      dsigk= dsigt;
      sigk= sigma + dsigk;
      
      for j=1:iteration
          f= DPyieldf(sigk, alpha, sigmay); %yield function value at current state
          
          if abs(f)<=tol %checking for convergence
              dsig= dsigk ;%Updating values after convergence
              sigmaf= sigma+dsig;
              depsp= depk;
              ddsdde= consis_mod(stif,sigmaf,alpha,beta); %Calling consistent tangent modulus
              break;
          end
          
          N= DPflowf(sigk, beta);%explicit estimation of plastic flow direction
          
          %Updating values now
          
          %Newton-Raphson- for updating dlambda
          x0=0;
          for k=1:iteration
              temp= x0;%storing previous iteration value
              x11= x0+h;%taking a point ahead current point
              x1= x0-h;%taking a point before current point
              
              dsigk= dsigt-(dlamk+x0)*stif*N;
              sigk= sigma + dsigk;
              f0= DPyieldf(sigk, alpha, sigmay);%yield function value at x0
              
              dsigk= dsigt-(dlamk+x1)*stif*N;
              sigk= sigma + dsigk;
              f1= DPyieldf(sigk, alpha, sigmay);%yield function value at x1
              
              dsigk= dsigt-(dlamk+x11)*stif*N;
              sigk= sigma + dsigk;
              f11= DPyieldf(sigk, alpha, sigmay);%yield function value at x11
              
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
          depk= dlamk*N;
          dsigk= stif*(deps - depk);
          sigk= sigma + dsigk;
          %Values ready for next iteration
          
          if j==iteration
              disp('Convergence Failure in D-P');
              return;
          end
      end
      
   end
          
           
      
end
