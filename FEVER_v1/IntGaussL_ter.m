function [ integrale ] = IntGaussL_ter(Lgt, f, ngauss, varargin)
%INTGAUSS15 integrale con regola di quadratura du Gauss (funzione valutata
%in come funzione delle coordinate baricentriche Li e Lj)
% %
% in Input: 
% Lgt  = lunghezza del dominio di integrazione 
% f  = funzione della variabile x da integrare. Se la funzione f è
% "vettorializzata", è più efficiente scommentare l'istruzione alla riga 53
% e commentare il ciclo for seguente.
% ngauss = numero di punti di Gauss. L'integrale fornisce risultati esatti
% a meno degli errori di troncamento per polinomi fino al grado 2*n - 1
% varargin = lista di argomenti opzionali (usati da f) 
% %
% Output:
% Int = result of the integration
%
%   Author: A. Cristofolini 17/09/2021
%----------------------------------------------------------------------------

switch (ngauss)
    case(1)
        w = 1;
        csi = [0.5 0.5];
    case(2)
        w = [0.5 0.5];
        csi = [0.78867513459481287 0.21132486540518713;
               0.21132486540518713 0.78867513459481287];
    case(3)
        w = [0.27777777777777778  0.44444444444444444  0.27777777777777778];
        csi = [0.88729833462074168 0.11270166537925832 ;
               0.5 0.5;
               0.11270166537925832 0.88729833462074168];
    case(4)
        w1 = 0.326072577431273;
        w2 = 0.173927422568727;
        w = [w2 w1 w1 w2];
        csi = [0.93056815579702629 0.06943184420297371;
               0.66999052179242810 0.33000947820757190;
               0.33000947820757190 0.66999052179242810;
               0.06943184420297371 0.93056815579702629];
%     case(5)
%         cstx1 = 0.906179845938664;
%         cstx2 = 0.538469310105683;
%         csth1 = 0.236926885056189;
%         csth2 = 0.478628670499366;
%         csth0 = 0.568888888888889;
%         w = [csth1 csth2 csth0 csth2 csth1];
%         csi = [-cstx1 -cstx2 0 cstx2 cstx1];
end




% integrale = w*(f(xgauss))';
integrale = 0;
for k = 1: ngauss
    integrale = integrale + w(k)*f(csi(k,:), varargin{:});
end

integrale = integrale * Lgt;

end

