function [Kel] = KelG_ter(Lg, xyg, MatKind, Pdata, gradl2, ifield, order)
% Calcolo della matrice di elemento per una terna di coordinate d'area L



[Mprop, ~] = MatLib_ter(MatKind, Pdata.iprop, xyg, Lg, ifield, true, Pdata.invflag);
if Pdata.RZflag
    Kel1 = xyg(2) * Mprop * gradl2;
else
    Kel1 = Mprop * gradl2;
end

switch order
    case 1
        Kel = Kel1;
    case 2
        [gradNL] = GradN_L(Lg, 2);
        Kel = gradNL' * Kel1 * gradNL;
end

end

