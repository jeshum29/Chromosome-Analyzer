function R = rotMatrix_Y(th)

R = [...
    [cos(th),      0,      sin(th)];...
    [0,              1,          0];...
    [-sin(th),     0,      cos(th)];...
    ];

R(4,4) = 1;

