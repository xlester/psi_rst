clear all
close all

load ExperimentalResultsSorted.mat
I = I2(:,:,1:2:7);
[M,N,k]=size(I);

A = makeParabola(M,N,1);
B = makeGausiana(M,N,1,200);

for n=1:k
    I(:,:,n) = A + B .* I(:,:,n);
end


%% Inicializando parametros del metodo RST.

Muestreo = 8; % Numero de pixeles a satar para el muestreo.
iters1   = 20; % Numero de iteraciones para el metodo completo.
iters2   = 60; % Numero de iteraciones para el calculo de los pasos.
lambdaDC = 450; % Parametro de regulacizacion para el DC
lambdaSC = 500; % Parametro de regulacizacion para Seno y Coseno.
%% Inicializando parametros del metodo AIA.

iters = 8;
v     = 1;
Sk    = sin( v* (0:1:k-1) );
Ck    = cos( v* (0:1:k-1) );
Show  = 1; % 1 si se decea mostrar resultados parciales.
%% Aplicando metodos RST y AIA

% Aplicando algoritmo RST.
tic
[pasosRST f_RST S C a] = RST(I,Sk,Ck,lambdaDC,lambdaSC,Muestreo,iters1,iters2,Show);
tRST = toc;
pasosRST=AntiAliasing(pasosRST);

% Aplicando algoritmo AIA.
tic
[pasosAIA f_AIA] = AIA(I,Sk,Ck,iters,Show);
tAIA = toc;
pasosAIA=AntiAliasing(pasosAIA);
%% Eliminando Piston de fase.
pasosRST = pasosRST-pasosRST(1);
Sk = sin(pasosRST);
Ck = cos(pasosRST);

[a1 f_RST] = MinCuaCpp(I,Sk,Ck);
f_RSTreg = f_RST;
pasosRST = atan2(Sk,Ck);

pasosAIA = pasosAIA-pasosAIA(1);
% Sk = sin(pasosAIA);
% Ck = cos(pasosAIA);
% [a1 f_AIA] = MinCuaCpp(I,Sk,Ck);
% pasosAIA = atan2(Sk,Ck);

%% Mostrando Resultados.

SP_RST    = angle(f_RST);
SP_AIA    = angle(f_AIA);


%figure,imshow(-SP_RSTreg,[]),title('fase Estimada RST regularizada');
figure,imshow(SP_RST,[]),title('fase Estimada RST');
figure,imshow(SP_AIA,[]),title('fase Estimada AIA');
figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');

disp('Estimados AIA');
disp(pasosAIA);

disp('Estimados RST');
disp(pasosRST);