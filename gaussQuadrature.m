function integral = gaussQuadrature(integrand,p)
    switch p
        case 0 %Not actual quadrature; just constant val*area
            integral = integrand(NaN,NaN)/2;
            return
        case 1 %might be sufficiently low error order introduced 
            u = [1/3 1/3];
            % w = 1;
            w = 1/2;
        case 3
            u = [1/6 1/6; 2/3 1/6; 1/6 2/3];
            w = [1/6 1/6 1/6];
    end
    % switch p
    %     case 1
    %         u = 0;
    %         w = 2; 
    %     case 2
    %         u = [-sqrt(1/3); sqrt(1/3)];
    %         w = [1 1];
    %     case 3
    %         u = [-sqrt(3/5); 0; sqrt(3/5)];
    %         w = [5/9 8/9 5/9];    
    %     case 4
    %         u = [-sqrt(3/7-2/7*sqrt(5/6)) -sqrt(3/7+2/7*sqrt(5/6)) sqrt(3/7-2/7*sqrt(5/6)) sqrt(3/7+2/7*sqrt(5/6))];
    %         w = [1/2+sqrt(30)/36 1/2-sqrt(30)/36 1/2+sqrt(30)/36 1/2-sqrt(30)/36];
    % end
    
    integral = 0;
    for k = 1:p
        integral = integral+w(k)*integrand(u(k,1),u(k,2));
    end
end