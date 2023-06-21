function [Mprop, Mdeg, dprop_dfield] = MatLib(MatKind, Pdata, xyg, Lg, varfield, ifield, evalflg, derflag)

iprop = Pdata.iprop;
invflag = Pdata.invflag;

% Material library
% prop(1) = material relative permittivity [Adim]
% prop(2) = material relative permeability [Adim]
% prop(3) = material electric conductivity [S/m]
% prop(4) = material thermal conductivity [W/(m K)]
%
% On input
% MatKind = tipo di materiale
%     "air" 
%     "aluminium"
%     "semicond"
%     "iron"
%     "copper"
%     "porcelain"
%     "Teflon"
%     "MagSteel" 
%     "Test1"
%     "Test2"
%     "Test3"
%     "Customxy" 
%     "XLPE" 
%     "EPDM"
%     "EPR"
%     "IsolanteScarso"
% iprop = indice della proprietà richiesta
%  iprop = 1 --> material relative permittivity [Adim]
%  iprop = 2 --> material relative permeability [Adim]
%  iprop = 3 --> material electric conductivity [S/m]
%  iprop = 4 --> material thermal conductivity [W/(m K)]

%  xyg = coordinate x;y (r;z) del (dei) punto/i dove si vuole valutare la proprietà richiesta 
%  Lg = coordinate baricentriche del (dei) punto/i dove si vuole valutare la proprietà richiesta 
%  GrdN = matrice gradiente delle funzioni di forma dell'elemento nel quale 
%   è contenuto il punto sul quale si valuta la proprietà 
%  varfield = valore dell'incognita del problema ai nodi dell'elemento nel quale 
%   è contenuto il punto sul quale si valuta la proprietà 
%  ifield = valore dei campi aggiuntivi in input ai nodi dell'elemento nel quale 
%   è contenuto il punto sul quale si valuta la proprietà  
%  evalflg se  = true --> viene valutato il valore della proprietà
%  invflag se  = true --> viene restituito in output il reciproco del valore della proprietà 
%  derflag se  = true --> viene valutata la derivata della proprietà

%   In output
% Mdeg = grado della funzione che descrive la proprietà
% Mgeg = 0 -> funzione costante
% Mgeg > 0 -> grado della funzione polinomiale espressa in termini di variabili
% spaziali
% Mgeg = -1 -> grado della funzione pari al grado delle funzioni di forma
% utilizzate
% dprop_dfield = derivata della proprietà rispetto ai campi in input e alla
% variabile


dprop_dfield = 0;

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
                    F2 = varfield;           %modulo campo elettrico
                    Mprop = sigma0 *exp(alpha*(F1 - T0) + beta*(F2-E0));
                    if derflag
                        dprop_dfield = [alpha beta]*Mprop;
                    end
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
                    beta = 0.0645*1e-6;
                    T0 = 273.15;
                    E0 = 0;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = varfield;           %modulo campo elettrico
                    Mprop = sigma0 * exp(alpha*(F1 - T0) + beta*(F2-E0));
                    if derflag
                        dprop_dfield = [beta alpha]*Mprop;
                    end
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
                    Ea = Ea*1.60218e-19; %conversione da elettronvolt a Joule
                    a = 4.8e-10;
                    b=-5.1e-8;
                    alpha = -1;%-1.42;
                    kb = 1.380649e-23;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = varfield;           %modulo campo elettrico
                    beta = a * F1 + b;
                    if varfield ~= 0
                        Mprop = A * exp(-Ea/(kb * F1)) * sinh(beta*F2) * F2^(alpha);
                        if derflag
                            d1 = Ea/(kb * F1 * F1) * Mprop + a * A * exp(-Ea/(kb * F1)) * cosh(beta*F2) * F2^(alpha);
                            d2 = A * exp(-Ea/(kb * F1)) * (alpha * sinh(beta*F2) * F2^(alpha - 1) + beta * cosh(beta*F2) * F2^(alpha));
                            dprop_dfield = [d2 d1];
                        end
                    else
                        Mprop = A * exp(-Ea/(kb * F1)) * beta;
                        if derflag
                            d1 =0;
                            d2 = 0;
                            dprop_dfield = [d2 d1];
                        end
                     
                    end
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
        
        
    case "EPR"
        sigma0 = 1e-12;
        Mdegs = [0, 0, 3, 0];
        Mprops = [2.24, 1, sigma0, 0.3];
        if Mdegs(iprop)> 0 && evalflg
            switch iprop
                case 3
                    alpha = 1.e-5;
                    beta = 1.e-5;
                    T0 = 293.15;
                    E0 = 1.8e7;
                    F1 = Lg * ifield(:,1);   %temperatura
                    F2 = varfield;           %modulo campo elettrico
                    Mprop0 = sigma0 *exp(alpha*(F1 - T0) + beta*(F2-E0));
                    Mprop = Mprop0+1e-12;
                    if derflag
                        dprop_dfield = [beta alpha]*Mprop0;
                    end
            end
        else
            Mprop = Mprops(iprop);
        end
        Mdeg = Mdegs(iprop);
    case "IsolanteScarso"
        Mdegs = [0, 0, 0, 0];
        Mprops = [1, 1, 1.e-12, 0.5];
        Mprop = Mprops(iprop);
        Mdeg = Mdegs(iprop);
end
if invflag
    Mprop = 1/Mprop;
end
end

