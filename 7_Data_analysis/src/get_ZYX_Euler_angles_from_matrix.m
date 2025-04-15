% get_GENFIRE_Euler_angles_from_matrix
% Xuezeng: for small phi and psi angles, tested no problem. 

% output parameter: new_Angles: 3x1 array of phi, theta, psi angles
% input parameters: rotMAT: 3x3 rotation matrix   


function new_Angles = get_ZYX_Euler_angles_from_matrix(rotMAT)

zero_tolerance = 1e-5;

new_Angles = zeros(3,1);

new_Angles(2) = asind(-1*rotMAT(3,1));

if abs(cosd(new_Angles(2))) > zero_tolerance
    new_Angles(1) = asind(rotMAT(2,1)/sqrt(1-rotMAT(3,1).^2));
    new_Angles(3) = asind(rotMAT(3,2)/sqrt(1-rotMAT(3,1).^2));
elseif acosd(rotMAT(1,3)) >= 90
    new_Angles(1) = (asind(rotMAT(2,3))+acosd(rotMAT(1,3))-180)/2;
    new_Angles(3) = (-asind(rotMAT(2,3))+acosd(rotMAT(1,3))-180)/2;
elseif acosd(rotMAT(1,3)) < 90
    new_Angles(1) = (asind(rotMAT(2,3))+acosd(rotMAT(1,3)))/2;
    new_Angles(3) = (-asind(rotMAT(2,3))+acosd(rotMAT(1,3)))/2;
end


% 
% if abs(cosd(new_Angles(2))) > zero_tolerance
%     new_Angles(1) = atan2(rotMAT(2,1),rotMAT(1,1));
%     new_Angles(3) = atan2(rotMAT(3,2),rotMAT(3,3));
% else
% end



% if nargin > 1
%     ZeroCenterConv = varargin{1};
% else
%     ZeroCenterConv = 0;
% end

% 
% if abs(sind(new_Angles(2))) > zero_tolerance
%     new_Angles(1) = atan2(rotMAT(2,3),rotMAT(1,3))*180/pi;
%     new_Angles(3) = atan2(rotMAT(3,2),-1*rotMAT(3,1))*180/pi;
% else
%     new_Angles(1) = atan2(rotMAT(3,3)*rotMAT(2,1),rotMAT(3,3)*rotMAT(1,1))*180/pi;
%     new_Angles(3) = 0;
% end
% 
% 
% if (abs(new_Angles(1)) > 90 || abs(new_Angles(3)) > 90 ) && ZeroCenterConv
%     MultiFactor = -1;
%     
%     new_Angles(2) = acosd(rotMAT(3,3))*MultiFactor;
% 
%     if abs(sind(new_Angles(2))) > zero_tolerance
%         new_Angles(1) = atan2(rotMAT(2,3)*MultiFactor,rotMAT(1,3)*MultiFactor)*180/pi;
%         new_Angles(3) = atan2(rotMAT(3,2)*MultiFactor,-1*rotMAT(3,1)*MultiFactor)*180/pi;
%     else
%         new_Angles(1) = atan2(rotMAT(2,1)*rotMAT(3,3),rotMAT(1,1)*rotMAT(3,3))*180/pi;
%         new_Angles(3) = 0;
%     end
% 
% end

end
