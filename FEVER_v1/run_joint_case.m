clear
close all
giuntocoperto
p0 = msh.POS;
p = p0(:,1:2);
t0 = msh.TRIANGLES;
t = t0(:,1:3);
ireg = t0(:,4);
ed = msh.LINES;
bedges = ed(:,1:2);
ibe = ed(:,3);

bcflag_e = zeros(length(bedges),1);
bcval_e = zeros(length(bedges),2);


%trovo i lati sui vari contorni e poi li concateno USARE QUESTO CON MESH
%DEL GIUNTO
tic;
ibeW12 = find(ibe==12);
ibeW13 = find(ibe==13);
ibeW14 = find(ibe==14);
ibeW15 = find(ibe==15);
ibeW16 = find(ibe==16);
ibeW17 = find(ibe==17);
ibeW18 = find(ibe==18);




ibeS19 = find(ibe==19);
ibeS20 = find(ibe==20);

ibeE21 = find(ibe==21);
ibeE22 = find(ibe==22);
ibeE23 = find(ibe==23);
ibeE24 = find(ibe==24);
ibeE25 = find(ibe==25);
ibeE26 = find(ibe==26);

ibeN27 = find(ibe==27);
ibeN29 = find(ibe==29);
ibeN30 = find(ibe==30);
ibeN32 = find(ibe==32);


ibeObl28 = find(ibe==28);
ibeObl31 = find(ibe==31);


ibeS = [ibeS19; ibeS20];   %  Trova tutti i lati sul contorno sud. I tag sono diversi, perchè sono diversi i materiali
ibeE = [ibeE21 ;ibeE22; ibeE23; ibeE24; ibeE25; ibeE26];   %  Trova tutti i lati sul contorno con tag (physical group assegnato in gmsh) = 5 (macro-lato est)
ibeN = [ibeN27; ibeN29; ibeN30; ibeN32];   %  Trova tutti i lati sul contorno con tag (physical group assegnato in gmsh) = 6 (macro-lato nord)
ibeW = [ibeW12; ibeW13; ibeW14; ibeW15; ibeW16; ibeW17; ibeW18];

ibeObl = [ibeObl28; ibeObl31];



Keff = 0.955;
V = (0.02)*(0.02) * pi * 4 * Keff; % pi * 0.24
Pt =  0.3065*1e5 * V; %J^2/sigma è la sorgente

RTt3 = Rsoil(1.3, 64, 1); %l = 1 tratto del giunto a raggio minore
%RTt2 = Rsoil(1.3, 77.5, 1); non uso il raggio medio per il tratto obliquo
RTt1 = Rsoil(1.3, 95, 1); %tratto del giunto a raggio maggiore

RTt4 = Rsoil(1.3, 58, 3.77); %tratto del cavo

%si calcolano le conduttanze per unità di lunghezza
g1 = 1/RTt1; %raggio maggiore
g3 = 1/RTt3; %raggio minore
g4 = 1/RTt4; %cavo


h1 = g1/(2*pi*0.095);
h3 = g3/(2*pi*0.064);
h4 = g4/(2* pi * 0.058);

Text = 293; %Temperatura esterna in kelvin

%condizioni al contorno per la parte con raggio variabile

for i=1:length(ibeObl)
    bcflag_e(ibeObl(i)) = 2;
    bcval_e(ibeObl(i),1) = Text;
    raggio_m_el = (p(bedges(ibeObl(i), 1), 2) + p(bedges(ibeObl(i), 2), 2))/2;%raggio medio dell'elemento sul contorno considerato
    
    raggio_m_el_mm = raggio_m_el * 1000;
    
    bcval_e(ibeObl(i),2) = (1/Rsoil(1.3, raggio_m_el_mm, 1))/(2*pi*raggio_m_el);
end

% iregbe = ones(size(tagbe));

bcflag_e(ibeS) = 1 * ones(length(ibeS),1);  %assegna il flag per la condizione al controno sul macro-lato sud
%bcval_e(ibeS,:) = 0 * ones(length(ibeS),2);  %assegna il valore per la condizione al controno sul macro-lato sud
bcval_e(ibeS,1) =  0 * ones(length(ibeS),1);  %assegna il valore per la condizione al controno sul macro-lato nord
bcval_e(ibeS,2) =  0 * ones(length(ibeS),1);  %assegna il valore per la condizione al controno sul macro-lato nord

bcflag_e(ibeE) = 1 * ones(length(ibeE),1);  %assegna il flag per la condizione al controno sul macro-lato est
%bcval_e(ibeE,:) =  0 * ones(length(ibeE),2);  %assegna il valore per la condizione al controno sul macro-lato est
bcval_e(ibeE,1) =  0 * ones(length(ibeE),1);  %assegna il valore per la condizione al controno sul macro-lato nord
bcval_e(ibeE,2) =  0 * ones(length(ibeE),1);  %assegna il valore per la condizione al controno sul macro-lato nord


%condizioni al contorno per il lato nord del cavo
bcflag_e(ibeN27) = 2 * ones(length(ibeN27),1);  %assegna il flag per la condizione al controno sul macro-lato nord
bcval_e(ibeN27,1) =  Text * ones(length(ibeN27),1);  %assegna il valore per la condizione al controno sul macro-lato nord
bcval_e(ibeN27,2) =  h4 * ones(length(ibeN27),1);  %assegna il valore per la condizione al controno sul macro-lato nord


%condizioni al contorno per la parte a raggio minore del giunto (lato nord)
bcflag_e(ibeN29) = 2 * ones(length(ibeN29),1);  
bcval_e(ibeN29,1) =  Text * ones(length(ibeN29),1);  
bcval_e(ibeN29,2) =  h3 * ones(length(ibeN29),1);  

bcflag_e(ibeN30) = 2 * ones(length(ibeN30),1);  
bcval_e(ibeN30,1) =  Text * ones(length(ibeN30),1);  
bcval_e(ibeN30,2) =  h3 * ones(length(ibeN30),1);  

%condizioni al contorno per la parte a raggio maggiore del giunto (lato nord)
bcflag_e(ibeN32) = 2 * ones(length(ibeN32),1);  
bcval_e(ibeN32,1) =  Text * ones(length(ibeN32),1);  
bcval_e(ibeN32,2) =  h1 * ones(length(ibeN32),1);  



bcflag_e(ibeW) = 1 * ones(length(ibeW),1);  %assegna il flag per la condizione al controno sul macro-lato ovest
%bcval_e(ibeW,:) =  0 * ones(length(ibeW),2);  %assegna il valore per la condizione al controno sul macro-lato ovest
bcval_e(ibeW,1) =  0 * ones(length(ibeW),1);  %assegna il valore per la condizione al controno sul macro-lato nord
bcval_e(ibeW,2) =  0 * ones(length(ibeW),1);  %assegna il valore per la condizione al controno sul macro-lato nord


ifield = p(:,1);



TipoProblema = "thermalRZ"; %electrostaticRZ electrodynamicSSRZ
sorgente = @src_ter;
% freq = 0;
order = 1;
materiali = ["copper", "copper", "semicond", "semicond", "XLPE", "semicond", "EPDM", "semicond", "semicond",       "aluminium", "XLPE",     "XLPE",      "aluminium", "semicond", "EPDM", "semicond", "semicond", "copper", "copper", "copper", "copper", "semicond", "XLPE", "semicond",       "aluminium", "XLPE", "XLPE", "XLPE", "XLPE", "XLPE", "XLPE", "XLPE", "EPDM", "XLPE"];
% materiali = ["copper", "iron", "copper", "iron", "iron", "iron","copper", "copper"];
% materiali = ["Teflon", "air", "Teflon", "air", "air", "air","Teflon", "Teflon"];
[ T, field, Kg ] = FEM2D_EM2020_ter( TipoProblema, p, t, bedges, bcflag_e, bcval_e, ireg, ibe, materiali, sorgente, order, ifield);

trisurf(t,p(:,1),p(:,2),T(:,1),'edgecolor','k','facecolor','interp','EdgeColor','none');
colorbar
axis equal  
view(2)

toc;
tic;

diel1=find(ireg==5);
diel2=find(ireg==7);
diel=[diel1;diel2];

tdiel = t(diel,:); %triangoli che appartengono ai dielettrici
ireg_diel=msh.TRIANGLES(diel,4);
vert_diel=unique(tdiel); %vertici che appartengono ai dielettrici

pdiel=p(vert_diel,:);
Tdiel = T(vert_diel); %temperatura nel dielettrico

ifield_diel(:,1) = Tdiel; %temperatura nei dielettrici
Ein = 18 * ones(length(vert_diel),1); %valore iniziale di campo elettrico

ifield_diel(:,2) = Ein;

for i=1:length(vert_diel)
    for j=1:3
        for k=1:length(tdiel)
            if vert_diel(i) == tdiel(k,j)
                tdiel(k,j) = i;
            end
        end
    end
end


ibediel33=find(ibe==33); %contorno dei dielettrici ancora da aggiungere in gmsh
ibediel34=find(ibe==34);
ibediel15=find(ibe==15);
ibediel23=find(ibe==23);

ibediel=[ibediel15;ibediel23;ibediel33;ibediel34];

ed_diel=msh.LINES(ibediel,:);
bedgesdiel=ed_diel(:,1:2); %lati sul contorno della seconda mesh
ibe_diel = ed_diel(:,3);

for i=1:length(vert_diel)
    for j=1:2
        for k=1:length(bedgesdiel)
            if vert_diel(i) == bedgesdiel(k,j)
                bedgesdiel(k,j) = i;
            end
        end
    end
end


bcflag_e_diel = zeros(length(bedgesdiel),1);
bcval_e_diel = zeros(length(bedgesdiel),2);

ibe_diel33=find(ibe_diel==33);
ibe_diel34=find(ibe_diel==34);
ibe_diel_Neum15 = find(ibe_diel==15);
ibe_diel_Neum23 = find(ibe_diel==23);

ibe_diel_Neum = [ibe_diel_Neum15;ibe_diel_Neum23];

bcflag_e_diel(ibe_diel33) = 2 * ones(length(ibe_diel33),1); %Dirichlet è due nel codice del prof, nel mio è 3. Da qui si usa il non linerare
bcval_e_diel(ibe_diel33,1) = 0 * ones(length(ibe_diel33),1);
bcval_e_diel(ibe_diel33,2) = 0 * ones(length(ibe_diel33),1);


bcflag_e_diel(ibe_diel34) = 2 * ones(length(ibe_diel34),1);
bcval_e_diel(ibe_diel34,1) = 550000 * ones(length(ibe_diel34),1);
bcval_e_diel(ibe_diel34,2) = 550000 * ones(length(ibe_diel34),1);

bcflag_e_diel(ibe_diel_Neum) = 1 * ones(length(ibe_diel_Neum),1);
bcval_e_diel(ibe_diel_Neum,1) = 0 * ones(length(ibe_diel_Neum),1);
bcval_e_diel(ibe_diel_Neum,2) = 0 * ones(length(ibe_diel_Neum),1);

sorgente = @src;
order = 1;
materiali = ["1", "2", "3", "4", "XLPE", "6", "EPDM", "8", "9", "10", "11", "12", "13", "14", "EPDM", "16", "17", "18", "19", "20", "21", "22", "XLPE", "24", "25", "26","27","28","29","30","31", "32", "EPDM", "XLPE"];

TipoProblema = "electrodynamicSSRZ";
toc; 

freq = 0;

%prova. Nel non lineare, i domini vano messi in ordine. Ad esempio, il 5 e
%7 devono diventare 1 e 2

[phi, field, IntSource, iter, err] = FEM2DNLv1( TipoProblema, pdiel, tdiel, bedgesdiel, bcflag_e_diel, bcval_e_diel, ireg_diel, ibe_diel, materiali, sorgente, freq, ifield_diel, order);



%plot del potenziale
figure
trisurf(tdiel,pdiel(:,1),pdiel(:,2),phi(:),'edgecolor','k','facecolor','interp','EdgeColor','none');
colorbar
axis equal  
view(2)

%plot di Ex
figure
trisurf(tdiel,pdiel(:,1),pdiel(:,2),field(:,1,1),'edgecolor','k','facecolor','interp','EdgeColor','none');
colorbar
axis equal  
view(2)

%plot di Ey
figure
trisurf(tdiel,pdiel(:,1),pdiel(:,2),field(:,2,1),'edgecolor','k','facecolor','interp','EdgeColor','none');
colorbar
axis equal  
view(2)


Emod = sqrt( field(:,1,1).^2 + field(:,2,1).^2 );
figure
trisurf(tdiel,pdiel(:,1),pdiel(:,2),Emod(:),'edgecolor','k','facecolor','interp','EdgeColor','none');
colorbar
axis equal
view(2)
% hold on
% quiver(p(:,1), p(:,2), field(:,1,3), field(:,2,3) ,1,'r')
% hold off

