function [consis_mod] = consis_mod(stif,sigmaf,alpha,beta)
%UNTITLED2 Summary of this function goes here
    I2 = [1; 1; 1; 0; 0; 0]; % 2nd order identity tensor in Voigt form
    sigma_m= (sigmaf(1)+sigmaf(2)+sigmaf(3))/3;
    dev= sigmaf - sigma_m*I2;
    dotprod= dot(dev,dev);
    deveq= sqrt((3/2)*dotprod);
    
    A= (3/(2*deveq))*dev + (beta/3)*I2;
    B= (3/(2*deveq))*dev + (alpha/3)*I2;
    
    denom= A'*stif*B;
    
    C= stif*B;
    D= (stif*A)';
    
    consis_mod= stif - (C*D)/denom;
    
end

