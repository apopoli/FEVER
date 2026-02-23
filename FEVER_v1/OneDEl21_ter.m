function [RHS1D, Kel1D] = OneDEl21_ter(pv, MatKind, Pdata, order, BCflag, BCval, ifield)
% calcolo della matrice di elemento (condizioni di Robin) e del termine
% noto di elemento per lati con condizioni al contorno di Robin e Neumann
% 

Length = norm(pv(1,:) - pv(2,:));



if BCflag == 1
    [~, Mdeg] = MatLib_ter(MatKind, Pdata.iprop, [0 0], 0, 0, false, false);  %ask for the polynomial degree ofthe function
    if Mdeg < 0
        Mdeg = -Mdeg * order;
    end
    if Pdata.RZflag
        deg = Mdeg + order + 1;
    else
        deg = Mdeg + order;
    end
    ngauss = ceil(0.5*(deg+1));
    Kel1D = zeros(2);
    RHS1D = IntGaussL_ter(Length, @RHSelNeum_ter, ngauss, pv, MatKind, order, Pdata, ifield)* BCval(1);
else
    if Pdata.RZflag
        deg = 2*order + 1;
    else
        deg = 2*order;
    end
    ngauss = ceil(0.5*(deg+1));
    Kel1D = IntGaussL_ter(Length, @KelRob_ter, ngauss, pv, order, Pdata) * BCval(2);
    if Pdata.RZflag
        deg = order + 1;
    else
        deg = order;
    end
    ngauss = ceil(0.5*(deg+1));
    RHS1D = IntGaussL_ter(Length, @RHSelRob_ter, ngauss, pv, order, Pdata) * BCval(1) * BCval(2);
end


end

