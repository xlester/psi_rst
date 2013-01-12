% Esta funcion genera una gausiana de amplitud A.
% M: Numero de Renglones.
% N: Numero de Columnas.
% A: Amplitud de la gausiana.
% s: Sigma
% Autor: Orlando Medina.
% fecha: 11 Enero 2013.

function [Gausiana] = makeGausiana(M,N,A,s)
    [X Y]    = meshgrid(1:N,1:M); 
    b        = (X-M/2).^2 + (Y-N/2).^2;
    Gausiana = A*exp(-b / (1*s*s) );
end