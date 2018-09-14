function [b, sinepsm] = magField(t,sc,mg)
% Updated 11/21/2017 by S. Hur-Diaz, Emergent
% Magnetic field in inertial frame
% Assumes time zero is when the spacecraft is at the equator?

Beta1 = mg.we*t;
cosepsm = cos(sc.In)*cos(mg.Gta)+sin(sc.In)*sin(mg.Gta)*cos(Beta1);
nim = atan2(-sin(mg.Gta)*sin(Beta1),...
            sin(sc.In)*cos(mg.Gta)-cos(sc.In)*sin(mg.Gta)*cos(Beta1));

if sin(nim)==0;
    sinepsm = sin(sc.In)*cos(mg.Gta)-cos(sc.In)*sin(mg.Gta)*cos(Beta1)/cos(nim);
else
    sinepsm = -sin(mg.Gta)*sin(Beta1)/sin(nim);
end

wot = sc.wo*t;
b = mg.kmg*[sinepsm*cos(wot-nim);
            -cosepsm;
            2*sinepsm*sin(wot-nim)];
% Tesla == kg * s^(-2) * A^(-1)