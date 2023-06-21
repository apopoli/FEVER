function [KJel] = KJelG(Lg, xyg, MatKind, Pdata, gradL, gradL2, phiv, ifield, flagJac)
% Calcolo della matrice di elemento per una terna di coordinate d'area L

% varargin : MatKind, Pdata, Grad_Lel, gradl2, phiv,  ifield , flagJac
switch Pdata.order
    case 1
        GradN = gradL;
        gradN2 = gradL2;
    case 2        
        [gradNL] = GradN_L(Lg, 2);
        GradN = gradL * gradNL;
        gradN2 = gradNL' * gradL2 * gradNL;
end

evalflg = true;
derflag = flagJac;


varfield  = norm(GradN * phiv);

[Mprop, ~, dprop_dfield] = MatLib(MatKind, Pdata, xyg, Lg, varfield, ifield, evalflg, derflag);

if Pdata.RZflag
    Kel1 = xyg(2) * Mprop * gradN2;
else
    Kel1 = Mprop * gradN2;
end

if flagJac
    a = gradN2 * phiv;
    if varfield ~= 0
        if Pdata.RZflag
            Jel = xyg(2) * dprop_dfield(1) * (a * a') / varfield;
        else
            Jel = dprop_dfield(1) * (a * a') / varfield;
        end
    else
        Jel=zeros(Pdata.Ksize);
    end
    KJel = [Kel1; Jel];
else
%     Jel = zeros(Pdata.Ksize);
    KJel = Kel1;
end

end

