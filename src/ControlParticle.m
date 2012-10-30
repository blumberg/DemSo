tick = tick + 1;

Rx = 1.5;
Ry = 0;
Ang = tick/1000;
Center = [5,5.5-tick/20000];

pos(1) = cos(Ang)*Rx + Center(1);
pos(2) = sin(Ang)*Ry + Center(2);
