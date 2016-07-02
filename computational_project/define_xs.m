function xs=define_xs
%
% UO2
%
% xs.uo2.t(1) = 2.12450E-01;
% xs.uo2.t(2) = 3.55470E-01;
% xs.uo2.t(3) = 4.85540E-01;
% xs.uo2.t(4) = 5.59400E-01;
% xs.uo2.t(5) = 3.18030E-01;
% xs.uo2.t(6) = 4.01460E-01;
% xs.uo2.t(7) = 5.70610E-01;

xs.uo2.t(1) = 1.77949E-01;
xs.uo2.t(2) = 3.29805E-01;
xs.uo2.t(3) = 4.80388E-01;
xs.uo2.t(4) = 5.54367E-01;
xs.uo2.t(5) = 3.11801E-01;
xs.uo2.t(6) = 3.95168E-01;
xs.uo2.t(7) = 5.64406E-01;

xs.uo2.tr = xs.uo2.t;
xs.uo2.D = 1 ./ (3 * xs.uo2.tr);

xs.uo2.a(1) = 8.02480E-03;
xs.uo2.a(2) = 3.71740E-03;
xs.uo2.a(3) = 2.67690E-02;
xs.uo2.a(4) = 9.62360E-02;
xs.uo2.a(5) = 3.00200E-02;
xs.uo2.a(6) = 1.11260E-01;
xs.uo2.a(7) = 2.82780E-01;

xs.uo2.c(1) = 8.12740E-04;
xs.uo2.c(2) = 2.89810E-03;
xs.uo2.c(3) = 2.03158E-02;
xs.uo2.c(4) = 7.76712E-02;
xs.uo2.c(5) = 1.22116E-02;
xs.uo2.c(6) = 2.82252E-02;
xs.uo2.c(7) = 6.67760E-02;

xs.uo2.f(1) = 7.21206E-03;
xs.uo2.f(2) = 8.19301E-04;
xs.uo2.f(3) = 6.45320E-03;
xs.uo2.f(4) = 1.85648E-02;
xs.uo2.f(5) = 1.78084E-02;
xs.uo2.f(6) = 8.30348E-02;
xs.uo2.f(7) = 2.16004E-01;

xs.uo2.nu(1) = 2.78145E+00;
xs.uo2.nu(2) = 2.47443E+00;
xs.uo2.nu(3) = 2.43383E+00;
xs.uo2.nu(4) = 2.43380E+00;
xs.uo2.nu(5) = 2.43380E+00;
xs.uo2.nu(6) = 2.43380E+00;
xs.uo2.nu(7) = 2.43380E+00;

xs.uo2.x(1) = 5.87910E-01;
xs.uo2.x(2) = 4.11760E-01;
xs.uo2.x(3) = 3.39060E-04;
xs.uo2.x(4) = 1.17610E-07;
xs.uo2.x(5) = 0.00000E+00;
xs.uo2.x(6) = 0.00000E+00;
xs.uo2.x(7) = 0.00000E+00;

xs.uo2.s(1,:) = [1.27537E-01 4.23780E-02 9.43740E-06 5.51630E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.uo2.s(2,:) = [0.00000E+00 3.24456E-01 1.63140E-03 3.14270E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.uo2.s(3,:) = [0.00000E+00 0.00000E+00 4.50940E-01 2.67920E-03 0.00000E+00 0.00000E+00 0.00000E+00];
xs.uo2.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 4.52565E-01 5.56640E-03 0.00000E+00 0.00000E+00];
xs.uo2.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 1.25250E-04 2.71401E-01 1.02550E-02 1.00210E-08];
xs.uo2.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.29680E-03 2.65802E-01 1.68090E-02];
xs.uo2.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 8.54580E-03 2.73080E-01];


%
% 4.3% MOX
%
% xs.m43.t(1) = 2.11920E-01;
% xs.m43.t(2) = 3.55810E-01;
% xs.m43.t(3) = 4.88900E-01;
% xs.m43.t(4) = 5.71940E-01;
% xs.m43.t(5) = 4.32390E-01;
% xs.m43.t(6) = 6.84950E-01;
% xs.m43.t(7) = 6.88910E-01;

xs.m43.t(1) = 1.78731E-01;
xs.m43.t(2) = 3.30849E-01;
xs.m43.t(3) = 4.83772E-01;
xs.m43.t(4) = 5.66922E-01;
xs.m43.t(5) = 4.26227E-01;
xs.m43.t(6) = 6.78997E-01;
xs.m43.t(7) = 6.82852E-01;

xs.m43.tr = xs.m43.t;
xs.m43.D = 1 ./ (3 * xs.m43.tr);

xs.m43.a(1) = 8.43390E-03;
xs.m43.a(2) = 3.75770E-03;
xs.m43.a(3) = 2.79700E-02;
xs.m43.a(4) = 1.04210E-01;
xs.m43.a(5) = 1.39940E-01;
xs.m43.a(6) = 4.09180E-01;
xs.m43.a(7) = 4.09350E-01;

xs.m43.c(1) = 8.06860E-04;
xs.m43.c(2) = 2.88080E-03;
xs.m43.c(3) = 2.22717E-02;
xs.m43.c(4) = 8.13228E-02;
xs.m43.c(5) = 1.29177E-01;
xs.m43.c(6) = 1.76423E-01;
xs.m43.c(7) = 1.60382E-01;

xs.m43.f(1) = 7.62704E-03;
xs.m43.f(2) = 8.76898E-04;
xs.m43.f(3) = 5.69835E-03;
xs.m43.f(4) = 2.28872E-02;
xs.m43.f(5) = 1.07635E-02;
xs.m43.f(6) = 2.32757E-01;
xs.m43.f(7) = 2.48968E-01;

xs.m43.nu(1) = 2.85209E+00;
xs.m43.nu(2) = 2.89099E+00;
xs.m43.nu(3) = 2.85486E+00;
xs.m43.nu(4) = 2.86073E+00;
xs.m43.nu(5) = 2.85447E+00;
xs.m43.nu(6) = 2.86415E+00;
xs.m43.nu(7) = 2.86780E+00;

xs.m43.x(1) = 5.87910E-01;
xs.m43.x(2) = 4.11760E-01;
xs.m43.x(3) = 3.39060E-04;
xs.m43.x(4) = 1.17610E-07;
xs.m43.x(5) = 0.00000E+00;
xs.m43.x(6) = 0.00000E+00;
xs.m43.x(7) = 0.00000E+00;

xs.m43.s(1,:) = [1.28876E-01 4.14130E-02 8.22900E-06 5.04050E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m43.s(2,:) = [0.00000E+00 3.25452E-01 1.63950E-03 1.59820E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m43.s(3,:) = [0.00000E+00 0.00000E+00 4.53188E-01 2.61420E-03 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m43.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 4.57173E-01 5.53940E-03 0.00000E+00 0.00000E+00];
xs.m43.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 1.60460E-04 2.76814E-01 9.31270E-03 9.16560E-09];
xs.m43.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.00510E-03 2.52962E-01 1.48500E-02];
xs.m43.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 8.49480E-03 2.65007E-01];


%
% 7.0% MOX
%
% xs.m70.t(1) = 2.14540E-01;
% xs.m70.t(2) = 3.59350E-01;
% xs.m70.t(3) = 4.98910E-01;
% xs.m70.t(4) = 5.96220E-01;
% xs.m70.t(5) = 4.80350E-01;
% xs.m70.t(6) = 8.39360E-01;
% xs.m70.t(7) = 8.59480E-01;

xs.m70.t(1) = 1.81323E-01;
xs.m70.t(2) = 3.34368E-01;
xs.m70.t(3) = 4.93785E-01;
xs.m70.t(4) = 5.91216E-01;
xs.m70.t(5) = 4.74198E-01;
xs.m70.t(6) = 8.33601E-01;
xs.m70.t(7) = 8.53603E-01;

xs.m70.tr = xs.m70.t;
xs.m70.D = 1 ./ (3 * xs.m70.tr);

xs.m70.a(1) = 9.06570E-03;
xs.m70.a(2) = 4.29670E-03;
xs.m70.a(3) = 3.28810E-02;
xs.m70.a(4) = 1.22030E-01;
xs.m70.a(5) = 1.82980E-01;
xs.m70.a(6) = 5.68460E-01;
xs.m70.a(7) = 5.85210E-01;

xs.m70.c(1) = 8.11240E-04;
xs.m70.c(2) = 2.97105E-03;
xs.m70.c(3) = 2.44594E-02;
xs.m70.c(4) = 8.91570E-02;
xs.m70.c(5) = 1.67016E-01;
xs.m70.c(6) = 2.44666E-01;
xs.m70.c(7) = 2.22407E-01;

xs.m70.f(1) = 8.25446E-03;
xs.m70.f(2) = 1.32565E-03;
xs.m70.f(3) = 8.42156E-03;
xs.m70.f(4) = 3.28730E-02;
xs.m70.f(5) = 1.59636E-02;
xs.m70.f(6) = 3.23794E-01;
xs.m70.f(7) = 3.62803E-01;

xs.m70.nu(1) = 2.88498E+00;
xs.m70.nu(2) = 2.91079E+00;
xs.m70.nu(3) = 2.86574E+00;
xs.m70.nu(4) = 2.87063E+00;
xs.m70.nu(5) = 2.86714E+00;
xs.m70.nu(6) = 2.86658E+00;
xs.m70.nu(7) = 2.87539E+00;

xs.m70.x(1) = 5.87910E-01;
xs.m70.x(2) = 4.11760E-01;
xs.m70.x(3) = 3.39060E-04;
xs.m70.x(4) = 1.17610E-07;
xs.m70.x(5) = 0.00000E+00;
xs.m70.x(6) = 0.00000E+00;
xs.m70.x(7) = 0.00000E+00;

xs.m70.s(1,:) = [1.30457E-01 4.17920E-02 8.51050E-06 5.13290E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m70.s(2,:) = [0.00000E+00 3.28428E-01 1.64360E-03 2.20170E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m70.s(3,:) = [0.00000E+00 0.00000E+00 4.58371E-01 2.53310E-03 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m70.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 4.63709E-01 5.47660E-03 0.00000E+00 0.00000E+00];
xs.m70.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 1.76190E-04 2.82313E-01 8.72890E-03 9.00160E-09];
xs.m70.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.27600E-03 2.49751E-01 1.31140E-02];
xs.m70.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 8.86450E-03 2.59529E-01];


%
% 8.7% MOX
%
% xs.m87.t(1) = 2.16280E-01;
% xs.m87.t(2) = 3.61700E-01;
% xs.m87.t(3) = 5.05630E-01;
% xs.m87.t(4) = 6.11170E-01;
% xs.m87.t(5) = 5.08900E-01;
% xs.m87.t(6) = 9.26670E-01;
% xs.m87.t(7) = 9.60990E-01;

xs.m87.t(1) = 1.83045E-01;
xs.m87.t(2) = 3.36705E-01;
xs.m87.t(3) = 5.00507E-01;
xs.m87.t(4) = 6.06174E-01;
xs.m87.t(5) = 5.02754E-01;
xs.m87.t(6) = 9.21028E-01;
xs.m87.t(7) = 9.55231E-01;

xs.m87.tr = xs.m87.t;
xs.m87.D = 1 ./ (3 * xs.m87.tr);

xs.m87.a(1) = 9.48620E-03;
xs.m87.a(2) = 4.65560E-03;
xs.m87.a(3) = 3.62400E-02;
xs.m87.a(4) = 1.32720E-01;
xs.m87.a(5) = 2.08400E-01;
xs.m87.a(6) = 6.58700E-01;
xs.m87.a(7) = 6.90170E-01;

xs.m87.c(1) = 8.14110E-04;
xs.m87.c(2) = 3.03134E-03;
xs.m87.c(3) = 2.59684E-02;
xs.m87.c(4) = 9.36753E-02;
xs.m87.c(5) = 1.89142E-01;
xs.m87.c(6) = 2.83812E-01;
xs.m87.c(7) = 2.59571E-01;

xs.m87.f(1) = 8.67209E-03;
xs.m87.f(2) = 1.62426E-03;
xs.m87.f(3) = 1.02716E-02;
xs.m87.f(4) = 3.90447E-02;
xs.m87.f(5) = 1.92576E-02;
xs.m87.f(6) = 3.74888E-01;
xs.m87.f(7) = 4.30599E-01;

xs.m87.nu(1) = 2.90426E+00;
xs.m87.nu(2) = 2.91795E+00;
xs.m87.nu(3) = 2.86986E+00;
xs.m87.nu(4) = 2.87491E+00;
xs.m87.nu(5) = 2.87175E+00;
xs.m87.nu(6) = 2.86752E+00;
xs.m87.nu(7) = 2.87808E+00;

xs.m87.x(1) = 5.87910E-01;
xs.m87.x(2) = 4.11760E-01;
xs.m87.x(3) = 3.39060E-04;
xs.m87.x(4) = 1.17610E-07;
xs.m87.x(5) = 0.00000E+00;
xs.m87.x(6) = 0.00000E+00;
xs.m87.x(7) = 0.00000E+00;

xs.m87.s(1,:) = [1.31504E-01 4.20460E-02 8.69720E-06 5.19380E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m87.s(2,:) = [0.00000E+00 3.30403E-01 1.64630E-03 2.60060E-09 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m87.s(3,:) = [0.00000E+00 0.00000E+00 4.61792E-01 2.47490E-03 0.00000E+00 0.00000E+00 0.00000E+00];
xs.m87.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 4.68021E-01 5.43300E-03 0.00000E+00 0.00000E+00];
xs.m87.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 1.85970E-04 2.85771E-01 8.39730E-03 8.92800E-09];
xs.m87.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.39160E-03 2.47614E-01 1.23220E-02];
xs.m87.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 8.96810E-03 2.56093E-01];

%
% Guide Tube
%
% xs.gt.t(1) = 1.90730E-01;
% xs.gt.t(2) = 4.56520E-01;
% xs.gt.t(3) = 6.40670E-01;
% xs.gt.t(4) = 6.49670E-01;
% xs.gt.t(5) = 6.70580E-01;
% xs.gt.t(6) = 8.75050E-01;
% xs.gt.t(7) = 1.43450E+00;

xs.gt.t(1) = 1.26032E-01;
xs.gt.t(2) = 2.93160E-01;
xs.gt.t(3) = 2.84240E-01;
xs.gt.t(4) = 2.80960E-01;
xs.gt.t(5) = 3.34440E-01;
xs.gt.t(6) = 5.65640E-01;
xs.gt.t(7) = 1.17215E+00;

xs.gt.tr = xs.gt.t;
xs.gt.D = 1 ./ (3 * xs.gt.tr);

xs.gt.a(1) = 5.11320E-04;
xs.gt.a(2) = 7.58010E-05;
xs.gt.a(3) = 3.15720E-04;
xs.gt.a(4) = 1.15820E-03;
xs.gt.a(5) = 3.39750E-03;
xs.gt.a(6) = 9.18780E-03;
xs.gt.a(7) = 2.32420E-02;

xs.gt.c(1) = 5.11320E-04;
xs.gt.c(2) = 7.58010E-05;
xs.gt.c(3) = 3.15720E-04;
xs.gt.c(4) = 1.15820E-03;
xs.gt.c(5) = 3.39750E-03;
xs.gt.c(6) = 9.18780E-03;
xs.gt.c(7) = 2.32420E-02;

xs.gt.s(1,:) = [6.61659E-02 5.90700E-02 2.83340E-04 1.46220E-06 2.06420E-08 0.00000E+00 0.00000E+00];
xs.gt.s(2,:) = [0.00000E+00 2.40377E-01 5.24350E-02 2.49900E-04 1.92390E-05 2.98750E-06 4.21400E-07];
xs.gt.s(3,:) = [0.00000E+00 0.00000E+00 1.83297E-01 9.23970E-02 6.94460E-03 1.08030E-03 2.05670E-04];
xs.gt.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 7.88511E-02 1.70140E-01 2.58810E-02 4.92970E-03];
xs.gt.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 3.73330E-05 9.97372E-02 2.06790E-01 2.44780E-02];
xs.gt.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 9.17260E-04 3.16765E-01 2.38770E-01];
xs.gt.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 4.97920E-02 1.09912E+00];

xs.gt.f(1:7) = 0;
xs.gt.nu(1:7) = 0;
xs.gt.x(1:7) = 0;


%
% Moderator
%
% xs.mod.t(1) = 2.30070E-01;
% xs.mod.t(2) = 7.76460E-01;
% xs.mod.t(3) = 1.48420E+00;
% xs.mod.t(4) = 1.50520E+00;
% xs.mod.t(5) = 1.55920E+00;
% xs.mod.t(6) = 2.02540E+00;
% xs.mod.t(7) = 3.30570E+00;

xs.mod.t(1) = 1.59206E-01;
xs.mod.t(2) = 4.12970E-01;
xs.mod.t(3) = 5.90310E-01;
xs.mod.t(4) = 5.84350E-01;
xs.mod.t(5) = 7.18000E-01;
xs.mod.t(6) = 1.25445E+00;
xs.mod.t(7) = 2.65038E+00;

xs.mod.tr = xs.mod.t;
xs.mod.D = 1 ./ (3 * xs.mod.tr);

xs.mod.a(1) = 6.01050E-04;
xs.mod.a(2) = 1.57930E-05;
xs.mod.a(3) = 3.37160E-04;
xs.mod.a(4) = 1.94060E-03;
xs.mod.a(5) = 5.74160E-03;
xs.mod.a(6) = 1.50010E-02;
xs.mod.a(7) = 3.72390E-02;

xs.mod.c(1) = 6.01050E-04;
xs.mod.c(2) = 1.57930E-05;
xs.mod.c(3) = 3.37160E-04;
xs.mod.c(4) = 1.94060E-03;
xs.mod.c(5) = 5.74160E-03;
xs.mod.c(6) = 1.50010E-02;
xs.mod.c(7) = 3.72390E-02;

xs.mod.s(1,:) = [4.44777E-02 1.13400E-01 7.23470E-04 3.74990E-06 5.31840E-08 0.00000E+00 0.00000E+00];
xs.mod.s(2,:) = [0.00000E+00 2.82334E-01 1.29940E-01 6.23400E-04 4.80020E-05 7.44860E-06 1.04550E-06];
xs.mod.s(3,:) = [0.00000E+00 0.00000E+00 3.45256E-01 2.24570E-01 1.69990E-02 2.64430E-03 5.03440E-04];
xs.mod.s(4,:) = [0.00000E+00 0.00000E+00 0.00000E+00 9.10284E-02 4.15510E-01 6.37320E-02 1.21390E-02];
xs.mod.s(5,:) = [0.00000E+00 0.00000E+00 0.00000E+00 7.14370E-05 1.39138E-01 5.11820E-01 6.12290E-02];
xs.mod.s(6,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.21570E-03 6.99913E-01 5.37320E-01];
xs.mod.s(7,:) = [0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.32440E-01 2.48070E+00];

xs.mod.f(1:7) = 0;
xs.mod.nu(1:7) = 0;
xs.mod.x(1:7) = 0;

N_g = length(xs.uo2.a);
for i = 1:N_g
	xs.uo2.r(i) = xs.uo2.t(i) - xs.uo2.s(i,i);
	xs.m43.r(i) = xs.m43.t(i) - xs.m43.s(i,i);
	xs.m70.r(i) = xs.m70.t(i) - xs.m70.s(i,i);
	xs.m87.r(i) = xs.m87.t(i) - xs.m87.s(i,i);
	xs.mod.r(i) = xs.mod.t(i) - xs.mod.s(i,i);
	xs.gt.r(i) = xs.gt.t(i) - xs.gt.s(i,i);
end

end