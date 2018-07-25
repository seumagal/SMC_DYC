	clear all
	load simin
	% ParÃ¢metros do modelo.

    a                   =  1.035;  %distância entre cg e eixo frontal
    b                   =  1.655;  %distância entre cg e eixo traseiro
    tf                  =  1.52;   %comprimento do eixo frontal
    tr                  =  1.54;   %comprimento do eixo traseiro
    hg                  =  0.542;  %altura do cg
    m                   =  1704.7; %massa
    hs                  =  0.451;  %altura do centro de rolagem
    ms                  =  1526.9; %massaa sobre centro de rolagem
    kf					=  47298/4;  %coeficiente de rigidez frontal a rolagem
    kr 					=  47311/4;  %coeficiente de rigidez traseira a rolagem
    cf					=  2823;   %coeficiente de amortecimento frontal a rolagem
    cr 					=  2652;   %coeficiente de amortecimental traseiro a rolagem
    roll_damping        =  cf+cr;   %coeficiente de rigidez a rolagem
    roll_stiffness      =  kf+kr;   %coeficiente de amortecimento a rolagem
    Izz                 =  3048.1;  %momento de inércia do eixo de guinada
    Ixx                 =  744;     %momento de inércia do eixo de rolagem
    Ixz                 =  0;      %produto de inércia dos eixos de rolagem e guinada
    g                   =  9.80665 ;    %aceleração gravitacional
	Reff                 = 0.307;
    MFA                 = [ 1.6  ,-34 , 1250 ,  2320  , 12.8 ,     0 , -0.0053 , 0.1925         ];
 	MFB                 = [ 1.55 ,0   , 1000 ,   60  , 300  , 0.17 , 0      , 0      ,0.2 ];
    Ku = 0.06;  %Coeficiente utilizado para calcular a taxa de guinada desejada
    uspeed = 80/3.6;               %velociade longitudinal de linearização

 %	MFA                 = [ 1.3  ,-49 , 1216 , 16320*5, 11    , 0.006  , -0.004, -0.4 ];
 %	MFB                 = [ 1.57 ,-48 , 1338 ,  5.8 , 30   , 0.17   , 0 , 0, 0.66 ];
    l = a+b;
	pneumatic_trail     = 0;
	u = 80/3.6;
	friction_coefficient = 0.9;
    J = 1;
    Ts = 8e-5;
    
	eta = m*Ixx*Izz - m*Ixz*Ixz - ms*ms*hs*hs*Izz;
    

	d(1) = ms*hs*Izz;
	d(2) = m*Ixz;
	d(3) = m*Izz*ms*hs*g;
	d(4) = m*Izz*roll_stiffness;
	d(5) = m*Izz*roll_damping;
	d(1:5) = d(1:5)/eta;
	d(6) = 1/Izz;
	d(7) = Ixz/Izz;
	d(8) = 1/m;
    d(9) = 1;
    d(10) = ms*hs/m;
	d(11) = 1/m;
	d(12) = ms*hs/m;
	d(13) = 1;

	zf0 = m*g*b/(2*l);
	zr0 = m*g*a/(2*l);
	z(1) = m*hg*b/(2*l);
	z(2) = m*hg/(l*tf);
	z(3) = kf/tf;
	z(4) = cf/tf;
	z(5) = kr/tr;
	z(6) = cr/tr;
	z(7) = m*hg*a/(2*l);
	tf2 = tf/2;
    tr2 = tr/2;
    C_1 = 4;
    M_1 = 0.5;
    M_2 = 0.5;
    GAMMA = 8;    
    
    %states_max = [14.55 23.80 0.35 0.01 0.02]
    %states_min = [-12.31 6.18 -0.26 -0.01 -0.02];
    dac_max = [3  ,17 , 0.5 , 0.2 , 0.2 ,1.5]
    dac_min = [-7 ,5  ,-1.5 ,-0.4 ,-0.4 ,-3];

	params = [MFA,MFB,z,zf0,zr0,a,b,tf2,tr2,d,Reff,J];
    paramssmc = [a,b,Izz,m,C_1,M_1 ,M_2,GAMMA];

    C_alpha = MFA(4)*sin(2*atan(([b,b,a,a]*m*g/(2000*l))/MFA(5)));
   
%% Modelo linear Mx'= A1x + B1u + E1steer
    M1 = [      m*uspeed,    0,  -ms*hs,  0;
                       0,  Izz,    -Ixz,  0;
           -ms*hs*uspeed, -Ixz,     Ixx,  roll_damping;
                       0,    0,       0,  1 ];

    A1 = [-sum(C_alpha)         , -(C_alpha*[a;a;-b;-b]/uspeed) - m*uspeed  , 0  , 0;
          -C_alpha*[a;a;-b;-b]  , -C_alpha*[a*a;a*a;b*b;b*b]              , 0  , 0;
            0                   , ms*hs*uspeed                            , 0  , (ms*hs*g-roll_stiffness);
            0                   , 0                                       , 1  , 0 ];

    B1 = [0; 1; 0; 0];  

    E1 = [ C_alpha*[1;1;0;0]; C_alpha*[a;a;0;0];0;0 ];
        
%% Modelo espaço de estado: x` = Ax + Bu + Esteer y'= Cx + Du
    A = M1\A1;
    B = M1\B1;
    C = eye(size(A));
    E = M1\E1;
    D = zeros(size([B,E]));

%% Modelo discreto
    %sys = ss( A,[B,E]  ,C,D);
    %d_sys = c2d(sys,8e-4);
	%Q = diag([4,1,0.001,0.001]);
	%R = 1e-8;
	%
    %LQRopt = lqrd(A,B,Q,R,8e-5)
    %LQRopt = lqr(A,B,Q,R)
    %clear sys
    %clear d_sys
    
	
    save modelParam
