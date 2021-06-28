
Point(1) = {0,0,0};
Point(2) = {1,0,0};
Line(1) = {1, 2};
out[] = Extrude{0,1,0}{ Line{1}; };

//+
Extrude {0, 0, 1} {
  Surface{5};
}
//+
Physical Surface(28) = {26, 18};
//+
Physical Surface(29) = {27, 22, 14, 5};
