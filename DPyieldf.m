function [f] = DPyieldf(sigma, alpha, sigy)
    I2 = [1; 1; 1; 0; 0; 0];
    tracesig= sigma(1) + sigma(2) + sigma(3);
    dev= sigma - (tracesig/3)*I2;
    dotprod= dot(dev,dev);
    deveq= sqrt((3/2)*dotprod);
    f= deveq + (alpha*tracesig/3) -sigy;

end