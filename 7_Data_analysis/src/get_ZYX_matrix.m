function Matrix = get_ZYX_matrix(Angles)
% According to WIKI on Euler angels, this is actually should be named as
% XYZ, because the rotation happens in XYZ order.

phi = Angles(1);
theta = Angles(2);
psi = Angles(3);

R1 = [cosd(phi) -sind(phi) 0;
      sind(phi) cosd(phi) 0;
      0          0         1];    
    
R2 = [cosd(theta) 0  sind(theta);
        0         1      0;
      -sind(theta) 0 cosd(theta)];

R3 = [1   0   0;
       0  cosd(psi) -sind(psi);
      0   sind(psi) cosd(psi)];    

Matrix = R1*R2*R3;

end