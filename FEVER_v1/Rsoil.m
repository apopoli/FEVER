function [RT] = Rsoil(bb,Rext,l)
%Rsoil calcola la resistenza termica del terreno per unità di lunghezza (bisogna dividere per la lunghezza del tratto per ottenre la resistenza termica effettiva) a seconda del raggio
%esterno del cavo Rext (in millimetri) e della profondità di interramento bb (in metri)

rho_t = 1.3; 

uu = bb * 1000/Rext;

RT = rho_t/(2*pi) * (log(uu + sqrt(uu*uu -1))); %resistenza termica per unità di lunghezza
%RT = RT/l; nel caso si voglia la resistenza totale
end

