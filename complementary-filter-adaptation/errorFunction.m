function [ error_mean ] = errorFunction( quaternion_sensor2earth, quaternionGT_sensor2earth )
%ERRORFUNCTION Summary of this function goes here
%   Detailed explanation goes here

    quaterionGT_sensor2earth_identity = repmat([1 0 0 0], length(quaternionGT_sensor2earth), 1);

    quaterionGT_sensor2earth_error = quaternProd(quaternion_sensor2earth, quatinv(quaternionGT_sensor2earth)) - quaterionGT_sensor2earth_identity;
    quaterionGT_sensor2earth_l2norm = sqrt(sum(quaterionGT_sensor2earth_error.^2, 2));
    quaterionGT_sensor2earth_l1norm = quatnorm(quaterionGT_sensor2earth_error);
    error_mean = mean(quaterionGT_sensor2earth_l2norm);
end

