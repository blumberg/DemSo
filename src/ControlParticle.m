tick = tick + 1;

Rx = 4;
Ry = 4;
Ang = tick/500;
Center = [5,5];

pos(1) = cos(Ang)*Rx + Center(1);
pos(2) = sin(Ang)*Ry + Center(2);
