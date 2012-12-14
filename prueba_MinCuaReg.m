%% Generando Interferogramas Sinteticos.
close all;
clear all;
M       = 512; % Number of rows of each interferogram.
N       = 512; % Number of columns of each interferogram.
k       = 5;   % Number of frames.
A       = 25;  % Amplitud para la fase tipo Peaks.

step    = pi/3; % Valor del paso.
nv      = 1.0; % Varianza del Ruido.

DC      = makeParabola(M,N,0.10);
rampa   = makeRampa(0.051,0.051,M,N);
phase   = makePeaks(N,M,A)+rampa;
b       = 1;
step_noise = 0.0;

[I,steps] = makeI(DC,b,phase,step,step_noise,k,nv);

%% Inicializando parametros para Minimos Cuadrados Regularizado.

iters = 200;
lambdaDC = 2; % Parametro de regulacizacion para el DC.
lambdaf  = 18; % Parametro de regulacizacion para campo complejo.

Phi   = ones(M,N);
Psi   = ones(M,N);
a     = ones(M,N);
f     = complex(Phi,Psi);
s = step;
Sk    = sin(s*(0:1:k-1));
Ck    = cos(s*(0:1:k-1));

%% Aplicando metodo Minimos Cuadrados Regularizado y Minimos Cuadrados Normal

[a f] = MinCuaReg(I,f,a,Sk,Ck,lambdaf,lambdaDC,iters);
[a1 f_Min] = MinCuaCpp(I,Sk,Ck);

%%  Mostrando los resultados
fase = angle(exp(1j*phase));
SP_reg   = angle(f);
SP_Min   = angle(f_Min);

errorFase = mean2(abs(fase-SP_reg));
disp('Pasos');
disp(steps);
disp('Error de fase');
disp(errorFase);
figure,imshow(I(:,:,1),[]),title('Interferograma de Entrada');
figure,imshow(a,[]),title('DC Calculado');
figure,imshow(fase,[]),title('Fase Esperada');
figure,imshow(SP_reg,[]),title('Fase Min Cua Regularizada');
figure,imshow(-SP_Min,[]),title('Fase Min Cua Sin Regularizar');

