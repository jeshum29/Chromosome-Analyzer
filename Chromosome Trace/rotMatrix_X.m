
function R = rotMatrix_X(th)

R = [...
    [1,        0,                 0];...
    [0,       cos(th),        -sin(th)];...
    [0,       sin(th),          cos(th)];...
    ];

R(4,4) = 1;






