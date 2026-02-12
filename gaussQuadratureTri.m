function integral = gaussQuadratureTri(integrand,p)
% Ben Sarfati 2/2026
    switch p
        case 0 %Not actual quadrature; just constant val*area
            integral = integrand(NaN,NaN)/2;
            return
        case 1 %might be sufficiently low error order introduced 
            u = [1/3 1/3];
            w = 1;
        case 3
            u = [1/6 1/6; 2/3 1/6; 1/6 2/3];
            w = [1 1 1]/3;
        case 4
            u = [1/3 1/3; 1/5 1/5; 3/5 1/5; 1/5 3/5];
            w = [-27 25 25 25]/48;
    end    
    integral = 0;
    for k = 1:p
        integral = integral+0.5*w(k)*integrand(u(k,1),u(k,2));
    end
end