function M = Mrotx(angle);

sa = sin(angle); ca = cos(angle);
M = [1 0 0 0; 0 ca sa 0; 0 -sa ca 0; 0 0 0 1];
