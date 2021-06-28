
Point(1) = {0,0,0};
Point(2) = {1,0,0};
Line(1) = {1, 2};
out[] = Extrude{0,1,0}{ Line{1}; };

//+
Extrude {0, 0, 1} {
  Surface{5};
}

Periodic Surface {18} = {26} Translate {1, 0, 0};
Periodic Surface {22} = {14} Translate {0, 1, 0};
Periodic Surface {27} = {5} Translate {0, 0, 1};

//+
Physical Surface(28) = {26, 18};
//+
Physical Surface(29) = {27, 22, 14, 5};
