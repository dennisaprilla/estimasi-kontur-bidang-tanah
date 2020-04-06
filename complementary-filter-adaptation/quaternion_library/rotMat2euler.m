function euler = rotMat2euler(R, varargin)
%ROTMAT2EULER Converts a rotation matrix orientation to ZYX Euler angles
%
%   euler = rotMat2euler(R)
%
%   Converts a rotation matrix orientation to ZYX Euler angles where phi is
%   a rotation around X, theta around Y and psi around Z.
%
%   For more information see:
%   http://www.x-io.co.uk/node/8#quaternions
%
%   Date          Author          Notes
%   27/09/2011    SOH Madgwick    Initial release

    if nargin == 1
        type = 'zyx';
    else
        if ischar( varargin{1} )
            type = varargin{1};
        else
            error(message('aero:rotMat2euler:notChar'));
        end
    end

    switch lower( type )
        case 'zyx'
            % ini yang kepake untuk algoritma madgwick
            phi = atan2(R(3,2,:), R(3,3,:) );
            theta = -atan(R(3,1,:) ./ sqrt(1-R(3,1,:).^2) );
            psi = atan2(R(2,1,:), R(1,1,:) );

        case 'zxy'
            phi = asin(R(3,2,:));
            psi = atan2(-R(1,2,:), R(2,2,:));
            theta = atan2(-R(3,1,:), R(3,3,:));

        case 'yxz'
            % ini yang kepake untuk algoritma dcm
            phi   = asin(-R(2,3,:));
            theta = atan2(R(1,3,:), R(3,3,:));
            psi   = atan2(R(2,1,:), R(2,2,:));

        otherwise
            error(message('aero:rotMat2euler:unknownRotation', type));

    end
    
    euler = [phi(1,:)' theta(1,:)' psi(1,:)'];
end

