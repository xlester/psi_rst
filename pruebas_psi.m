%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 256; % Number of rows of each interferogram.
N       = 256; % Number of columns of each interferogram.
k       = 5;   % Number of frames.
A       = 40;  % Amplitud para la fase tipo Peaks.

step    = pi/2; % Valor del paso.
nv      = 0.0; % Varianza del Ruido.

%DC      = makeParabola(M,N,55);
DC      = makeGausiana(M,N,5,60);
phase   = makeRampa(4*pi/M,4*pi/N,M,N);
%phase   = makePeaks(N,M,A)+rampa;
b       = makeGausiana(M,N,1,95);
step_noise = 0;

[I,steps]       = makeI(DC,b,phase,step,step_noise,k,nv);
steps = atan2(sin(steps),cos(steps));


%% Inicializando parametros del metodo RST.

Muestreo = 6; % Numero de pixeles a satar para el muestreo.
iters1   = 20; % Numero de iteraciones para el metodo completo.
iters2   = 50; % Numero de iteraciones para el calculo de los pasos.
lambdaDC = 00; % Parametro de regulacizacion para el DC
lambdaSC = 500; % Parametro de regulacizacion para Seno y Coseno.
lambda   = 00; % Parametro de regulacizacion.
%% Inicializando parametros del metodo AIA.

iters = 20;
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

%% Eliminando Piston de fase.
pasosRST = pasosRST-pasosRST(1);
Sk = sin(pasosRST);
Ck = cos(pasosRST);

[a1 f_RST] = MinCuaCpp(I,Sk,Ck);
% f_RSTreg = f_RST;
% for x=1:200
%     [a f_RSTreg] = MinCuaReg(I,f_RSTreg,a,Sk,Ck,1,15,3);
% end

pasosRST = atan2(Sk,Ck);

pasosAIA = pasosAIA-pasosAIA(1);
Sk = sin(pasosAIA);
Ck = cos(pasosAIA);
[a1 f_AIA] = MinCuaCpp(I,Sk,Ck);
pasosAIA = atan2(Sk,Ck);

%% Mostrando Resultados.

%SP_RSTreg = angle(f_RSTreg);
SP_RST    = angle(f_RST);
SP_AIA    = angle(f_AIA);
wfase     = angle(exp(-1i*phase));
errorReg = mean2(abs(wfase+SP_RST));

%figure,imshow(-SP_RSTreg,[]),title('fase Estimada RST regularizada');
figure,imshow(SP_RST,[]),title('fase Estimada RST');
figure,imshow(SP_AIA,[]),title('fase Estimada AIA');
figure,imshow(wfase,[]),title('fase Esperada');
figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');

disp('Estimados AIA');
disp(pasosAIA);

disp('Estimados RST');
disp(pasosRST);

disp('Esperados');
disp(steps);

disp('Error AIA');
disp(abs(steps - pasosAIA));

disp('Error RST');
disp(abs(steps - pasosRST));

disp('Error medio de fase');
disp(errorReg);