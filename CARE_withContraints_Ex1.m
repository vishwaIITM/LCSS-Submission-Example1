%%%% EXAMPLE 1 
%%%% 
%%%% Code to check if the solutions of the coupled Algebriac Riccati equations
%%%% satisfy the structural constraints (7)
%%%%  B_i'P_i(I-C_i'(C_iC_i')^{-1}C_i) = 0, i=1,2,3
%%%%
syms p11 p12 p13 p14 p15 p16 p21 p22 p23 p24 p25 p26  p31 p32 p33 p34 p35 p36 real 

%%%  Problem parameters
%%%
A=[1 0 0;0 1 0;0 0 1];
B1=[1;0;0];
B2=[0;1;0];
B3=[0;0;1];
C1=[1 0 0; 0 1 0]; C2=eye(3); C3=[0 1 0;0 0 1];
R1=1;
R2=1;
R3=1;
Q1=[1 -1 0;-1 1 0;0 0 0];
Q2=[1 -1 0;-1 2 -1;0 -1 1];
Q3=[0 0 0;0 1 -1;0 -1 1];

pvect1=[p11; p12; p13; p14; p15; p16];
pvect2=[p21; p22; p23; p24; p25; p26];
pvect3=[p31; p32; p33; p34; p35; p36];
%
%%% Solutions of the CARE
%%% 
P1=[p11 p12 p13;p12 p14 p15; p13 p15 p16];
P2=[p21 p22 p23;p22 p24 p25; p23 p25 p26];
P3=[p31 p32 p33;p32 p34 p35; p33 p35 p36];

%%%
%%%
S1=B1*inv(R1)*B1'; S2=B2*inv(R2)*B2'; S3=B3*inv(R3)*B3';
A1=A-S2*P2-S3*P3; 
A2=A-S1*P1-S3*P3; 
A3=A-S1*P1-S2*P2;
%
%%% coupled Algebriac Riccati equations (CARE)
%%%
are1=transpose(A1)*P1+P1*A1+Q1-P1*S1*P1;
are2=transpose(A2)*P2+P2*A2+Q2-P2*S2*P2;
are3=transpose(A3)*P3+P3*A3+Q3-P3*S3*P3;
%%%
%%% 
aretot=[vectgen(are1);vectgen(are2);vectgen(are3)];   %%% 18 polynomial equations of degree 2
pvect=[pvect1; pvect2; pvect3];                       %%% 18 variables
%
%
%
%%%% Imposing the structural conditions (7) on the CARE we get
%%%%  B_1'P_1(I-C_1'(C_1C_1')^{-1}C_1 =[0  0  p13]= 0 ==> p13 = 0 
%%%%  B_2'P_2(I-C_2'(C_2C_2')^{-1}C_2 =[0  0  0]
%%%%  B_3'P_3(I-C_3'(C_3C_3')^{-1}C_3 =[p33  0  0]= 0 ==> p33 = 0
aretot1=subs(aretot,{p13,p33},{0,0});                 %%% 18 polynomail equations of degree 2

%%%% The number of variables are 16 as p13=0 and p33=0

pvect1=[p11; p12; p14; p15; p16;...                   %%% 16 variables 
    p21; p22; p23; p24; p25; p26;...
    p31; p32; p34; p35; p36];

%%%% Solve the CARE setting p13=0 and p33=0

sol=solve(aretot1,pvect1);                           %%% solve equation (6) satisfying
                                                     %%% constraints (7)

sol                                                  %%% display solutions


%%% 
%%% script to vectorize the (symmetric) matrix equations
%%%
function vect=vectgen(mat)
%
N=length(mat); 
%
vect=[];
for k=1:N
    vect=[vect mat(k,k:N)];
end
vect=transpose(vect);
end    
