function [x, iter, err]= NR_finale1(NRfun, x0 , tol, itermax, varargin )
% dati input:
%   NRfun funzione definita dall'utente che restituisce in output
%       F= funzione vettoriale del sistema
%       Jac= matricejacobiana riferita a F
%   x0= vettore COLONNA approssimazione iniziale
%   tol= tolleranza
%   itermax= numero massimo di iterazioni
%   
% dati di output:
%   x= vettore che approssima la soluzione
%   iter= numero di iterazioni effettuate
%   err= vettore valore assoluto della differenza tra iter e iter-1


% impongo condizioni iniziali
x=x0;
flagJac = false;
[Fx0,~] =NRfun(x, flagJac, varargin{:});
nFx0 = norm(Fx0);
if nFx0 == 0
    return
end
nx0 = norm(x0);
if nx0 == 0
    nx0 = nFx0;
end


err = zeros(itermax,1);

%  Newton_Raphson

for iter = 1:itermax
    flagJac = true;
    [Fx, Jac] = NRfun(x, flagJac, varargin{:});
    
    delta = - Jac\Fx;
    x = x+delta;
    err(iter)=norm(Fx)/nFx0;

    Dnorm = norm(delta)/norm(x); 
    fprintf('iter %8d, F norm = %16.10E,  Delta norm = %16.10E\n',iter, err(iter), Dnorm);

    

      
    % verifica convergenza
    if (iter>=itermax)
        fprintf('NUMERO MASSIMO DI ITERAZIONI RAGGIUNTO');
    
    else
        if( (Dnorm<=tol)&& err(iter)<=tol)
                break;
        end
        
    end
end

  
end