function [Mprop, Mdeg] = MatLib(MatKind, iprop, xyg, Lg, ifield, evalflg, invflag)
% Material library
% prop(1) = material relative permittivity [Adim]
% prop(2) = material relative permeability [Adim]
% prop(3) = material electric conductivity [S/m]
% prop(4) = material thermal conductivity [W/(m K)]
% Mdeg = grado della funzione che descrive la proprietà
% Mgeg = 0 -> funzione costante
% Mgeg > 0 -> grado della funzione espressa in termini di variabili
% spaziali
% Mgeg = -1 -> grado della funzione pari al grado delle funzioni di forma
% utilizzate




switch MatKind
    case "air"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 1, 0, 2.62e-2];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "aluminium"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 1, 3.54e7, 237];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "semicond"
        Mdegs = [0, 0, 0, 0];
        Mprops = [2.30, 1, 6e3, 0.2857];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "iron"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 2000, 1.03e+007, 5.2e1];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "copper"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 0.999991, 5.8e+007, 3.8e2];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "porcelain"
        Mdegs = [0, 0, 0, 0];
        Mprops = [5.7, 1, 1.e-15, 1.5];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "Teflon"
        Mdegs = [0, 0, 0, 0];
        Mprops = [2.5, 1, 0, 2.5e-1]; %
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "MagSteel"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 2000, 2.e+006, 5.4e1];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "Test1"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1e8, 1, 1, 1];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
    case "Test2"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1e-8, 1e-12, 1e-12, 1e-12];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
        
        
    case "Test3"
        sigma0 = 1e-12;
        Mdegs = [0, 0, 3, 0];
        Mprops = [1e-8, 1e-12, sigma0, 1e-12];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3
                    alpha = 1.e-5;
                    beta = 1.e-5;
                    T0 = 20;
                    E0 = 1e6;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = Lg * ifield(:,2);   %modulo campo elettrico
                    Mprop = sigma0 *exp(alpha*(F1 - T0) + beta*(F2-E0));                   
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
    case "Customxy"
        sigma0 = 1;
        Mdegs = [0, 0, 1, 0];
        Mprops = [1e-8, 1e-12, sigma0, 1e-12];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3    
                    F1 = Lg * ifield(:,1);
                    Mprop = sigma0 *(1 + F1);%(xyg(:,1) + xyg(:,2));                   
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
        
        
    case "XLPE"
        sigma0 = 1e-16;
        Mdegs = [0, 0, 3, 0];
        Mprops = [2.3, 1, sigma0, 0.2857];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3
                    alpha = 0.084;
                    beta = 0.0645;
                    T0 = 273.15;
                    E0 = 0;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = Lg * ifield(:,2);   %modulo campo elettrico
                    Mprop = sigma0*exp(alpha*(F1 - T0) + beta*(F2-E0));                   
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
        
        
     case "EPDM"

         
        sigma0=1e-16; 
        Mdegs = [0, 0, 3, 0];
        Mprops = [2.24, 1, sigma0, 0.3];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3
                    A = 97;
                    Ea = 0.44;
                    Ea = Ea*1.60218e-19; %conversione da eletrtonvolt a Joule
                    a = 4.8e-10;
                    b=-5.1e-8;
                    alpha = -1.42;
                    kb = 1.380649e-23;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = Lg * ifield(:,2);   %modulo campo elettrico
                    beta = a * F1 + b;
                    Mprop = A * exp(-Ea/(kb * F1)) * sinh(beta*F2) * F2^(alpha);                   
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
        
        
    case "EPR"
        sigma0 = 1e-12;
        Mdegs = [0, 0, 3, 0];
        Mprops = [2.24, 1, 1e-12, 0.3];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3
                    alpha = 1.e-5;
                    beta = 1.e-5;
                    T0 = 293.15;
                    E0 = 1.8e7;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = Lg * ifield(:,2);   %modulo campo elettrico
                    Mprop = sigma0 *exp(alpha*(F1 - T0) + beta*(F2-E0));                   
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
end
if invflag
    Mprop = 1/Mprop;
end
end

