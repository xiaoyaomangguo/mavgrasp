function [desire_state] = getDesiredStateFromTrajectory(coefficient, time_stamp)
%���ǽ�ͨ��differential flatness ���õ���������ֵ��
%rq������xq,zq��, drq���ɻ������ٶ�dxq,dzq��, beta����е����ˮƽ�ĽǶȣ�, 
%theta��pitch��, F��������, M(����theta������)
% 
% waypoints = [0 4 40 0; 3 2 90 0; 4 1.5 120 0; 6 4 150 0];
% timepoints = [0 2.0 4.0 6.0];
% number_order = 7;
% mu_r = 1;
% mu_yaw = 1;

% basic_coef = getTrajectoryCoef(waypoints, timepoints, ...
%     number_order, mu_r, mu_yaw)
% for time_stamp = 0 : 0.1 : 2.0;
% coefficient = reshape(basic_coef(: , 1), 4, 8);           %�����4��8�еľ��󣬵�һ����x��ϵ��cx0-
m_q = 2;      %2kg
m_g = 0.2;    %0.2kg
m_s = (m_g + m_q);
L   = 0.1;      %the length of grasper
G = 9.8;
J_g = 1;
J_q = 1;

coefficient = coefficient';
t = [1; time_stamp; time_stamp.^2;  time_stamp.^3;  time_stamp.^4;  time_stamp.^5;  time_stamp.^6;  time_stamp.^7]';%������
% x z beta 
desired_state(1) = sum(coefficient(1, :).*t);   %x
desired_state(2) = sum(coefficient(2, :).*t);   %z
desired_state(3) = sum(coefficient(3, :).*t);   %beta
%d_x d_z d_beta
new_coef = polyder(coefficient(1, :));
desired_state(4) = sum(new_coef .*t(end - size(new_coef, 2) + 1 : end));   %d_x
new_coef = polyder(coefficient(2, :));
desired_state(5) = sum(new_coef .*t(end - size(new_coef, 2) + 1 : end));   %d_z
new_coef = polyder(coefficient(3, :));
desired_state(6) = sum(new_coef .*t(end - size(new_coef, 2) + 1 : end));   %d_beta
%dd_x dd_z dd_beta
new_coef = polyder(polyder(coefficient(1, :)));
desire_state(7) = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(coefficient(2, :)));
desire_state(8) = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(coefficient(3, :)));
desire_state(9) = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
%ddd_x ddd_z ddd_beta
new_coef = polyder(polyder(polyder(coefficient(1, :))));
d3xq = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(polyder(coefficient(2, :))));
d3zq = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(polyder(coefficient(3, :))));
d3beta = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
%dddd_x dddd_z  dddd_beta
new_coef = polyder(polyder(polyder(polyder(coefficient(1, :)))));
d4xq = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(polyder(polyder(coefficient(2, :)))));
d4zq = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));
new_coef = polyder(polyder(polyder(polyder(coefficient(3, :)))));
d4beta = sum(new_coef .* t(end - size(new_coef, 2) + 1 : end));

%%differential flatness
r_q = [desired_state(1); 0 ; desired_state(2)];
r_g = r_q + L * [cos(desired_state(3)); 0; -sin(desired_state(3))];
dr_q = [desired_state(4); 0; desired_state(5)];
dr_g = dr_q + L * [-sin(desired_state(3)) * desired_state(6); 
						0;  
						-cos(desired_state(3)) * desired_state(6) ] ;
d2r_q = [desired_state(7); 0; desired_state(8)];
d2r_g = d2r_q + L * [-cos(desired_state(3)) * desired_state(6)^2 - sin(desired_state(3)) * desire_state(9);  
							0; 
							sin(desired_state(3)) * desired_state(6)^2 + cos(desired_state(3)) * desire_state(9)];
d3r_q = [d3xq; d3zq];
d3r_g = d3r_q + L * [sin(desired_state(3)) * desired_state(6)^3 - 3 * cos(desired_state(3)) * desired_state(6) * desire_state(9) - sin(desired_state(3)) * d3beta;
                     0;
                     cos(desired_state(3)) * desired_state(6)^3 + 3 * sin(desired_state(3)) * desired_state(6) * desire_state(9) - cos(desired_state(3)) * d3beta];
d4r_q = [d4xq; d4zq];
d4r_g = d4r_q + L * [cos(desired_state(3)) * desired_state(6)^4 - 3 * cos(desired_state(3)) *  desire_state(9)^2 - 4 * cos(desired_state(3)) * desired_state(6) * d3beta - sin(desired_state(3)) * d4beta;
                                0;
                                -sin(desired_state(3)) * desired_state(6)^4 - 3 * sin(desired_state(3)) *  desire_state(9)^2 + 4 * sin(desired_state(3)) * desired_state(6) * d3beta - cos(desired_state(3)) * d4beta]
r_s = (m_q * r_q + m_g * r_g) / m_s;
dr_s = (m_q * dr_q + m_g * dr_g) / m_s;
d2r_s = (m_q * d2r_q + m_g * d2r_g) / m_s;
d3r_s = (m_q * d3r_q + m_g * d3r_g) / m_s;
d4r_s = (m_q * d4r_q + m_g * d4r_g) / m_s;


b2 = [0; 1; 0]; 
b3 = (d2r_s + G * [0; 0; 1] ) / (norm(d2r_s + G * [0; 0; 1]));
b1 = cross(b2, b3); 
F = m_s * norm(d2r_s + G * [0; 0; 1]);
dF = m_s * (b3 .* d3r_s);
theta = acos(dot(b3, [0; 0; 1]) / (norm(b3) * norm([0; 0; 1]))); % ����
dtheta = m_s * (b1 .* d3r_s) / F;
d2theta = (m_s * b1 .* d4r_s - 2 * dF * dtheta);

d2F = b3 .* (m_s * d4r_s) + dtheta^2 * F;
tao = J_g * desired_state(9) - L * m_g * (d2r_g(1,1) * sin(desired_state(3)) + (d2r_g(2,1) + G) * cos(desired_state(3)));
M = d2theta * J_q + tao;

end
