function [basis_coef] = getTrajectoryCoef(waypoints, timepoints, ...
    number_order, mu_r, mu_yaw)%, time_factor)
clc
%% Generate Coefficients of Polynomial Basis Function
% Generate Hessian Matrix
% P(t) = c0+c1*t+c2*t.^2+c3*t.^3+c4*t.^4+c5*t.^5+ .. + cn*t.^n
% nth derivative of P(t) = n!cn + (n+1)!/2!c(n+1)*t...
% = SIGMA( (i!/(i - n)! * ci * t.^ (i - n)) )
% square of nth derivative of P(t)
% coefficient of ci * cj
% i!j!/(i-n)!/(j-n)! * t.^(i + j - 2n) * ci * cj
% integral
% i!j!/(i-n)!/(j-n)! / (i + j - 2n + 1) * t.^(i + j - 2n + 1) * ci * cj

% problem: min( c'H c + f'c)
% c is 4 * n * m vector, [xTij yTij zTij yawTij]
% i is [1 .. n] j is [1 .. m]

%% make Hessian
% three waypoints only
waypoints = [0 4 30 0; 1.5 2.5 70 0; 3.0 2.0 90 0; 4.5 2.5 110 0; 6.0 4.0 140 0];
% waypoints = [0 4 30 0; 2.5 3.5 70 0; 3.0 2.0 90 0; 5.5 3.5 110 0; 6.0 4.0 140 0];
% waypoints = [0 4 30 0;  3.0 2.0 90 0;  6.0 4.0 140 0];

timepoints = [0 0.5 1.0 1.5 2.0];
number_order = 7;
mu_r = 1;
mu_yaw = 1;

m = size(waypoints, 1) - 1;
n = number_order; % if it is 6 than the biggest is .^6 with seven coefficients
H_all = zeros(4 * (n + 1) * m, 4 * (n + 1) * m);
A = [];
b = [];
num_segments = length(timepoints) - 1;

for id = 1 : num_segments
    % make Hessian for one segment
    t1 = timepoints(id);
    t2 = timepoints(id + 1);
    H = zeros(4 * (n + 1), 4 * (n + 1));
    H_xyz = zeros(n + 1, n + 1);
    H_yaw = zeros(n + 1, n + 1);
    order_diff = 4;
    for i = order_diff : n
        for j = order_diff : n
            H_xyz(i + 1, j + 1) = factorial(i) * factorial(j) / factorial(i - order_diff) / factorial(j - order_diff) ...
                / (i + j - 2 * order_diff + 1);
            H_xyz(i + 1, j + 1) = H_xyz(i + 1, j + 1) * (t2.^(i + j - 2 * order_diff + 1) - t1.^(i + j - 2 * order_diff + 1));
            H(i * 4 + 1, j * 4 + 1) = H_xyz(i + 1, j + 1);
            H(i * 4 + 2, j * 4 + 2) = H_xyz(i + 1, j + 1);
            H(i * 4 + 3, j * 4 + 3) = H_xyz(i + 1, j + 1);
        end
    end
    order_diff = 2;
    for i = order_diff : n
        for j = order_diff : n
            H_yaw(i + 1, j + 1) = factorial(i) * factorial(j) / factorial(i - order_diff) / factorial(j - order_diff) ...
                / (i + j - 2 * order_diff + 1);
            H_yaw(i + 1, j + 1) = H_yaw(i + 1, j + 1) * (t2.^(i + j - 2 * order_diff + 1) - t1.^(i + j - 2 * order_diff + 1));
            H(i * 4 + 4, j * 4 + 4) = H_yaw(i + 1, j + 1);
        end
    end
    % put into final Hession
    H_all((id - 1) * 4 * (n + 1) + 1 : id * 4 * (n + 1), (id - 1) * 4 * (n + 1) + 1 : id * 4 * (n + 1)) = H;
    
    % make constraints
    % every segment has 4 * 5 * 2 constraints , half of them are continuity
    % constraints
    % 5 means derivative of 0, 1, 2, 3, 4
    begin_base = []; end_base = [];
    for order = 0 : n
        begin_base= [t1.^order begin_base]; % descending
        end_base = [t2.^order end_base];
    end
    derivatives_begin = fliplr(begin_base); % ascending
    derivatives_end = fliplr(end_base);
    for i = 0 : 4
        begin_base = polyder(begin_base);   % diff based on descending
        end_base = polyder(end_base);
        derivatives_begin = [derivatives_begin; zeros(1, n + 1 - length(begin_base)) fliplr(begin_base)];% ascending,升序，与前阶导数后面对齐，前面补零
        derivatives_end = [derivatives_end; zeros(1, n + 1 - length(end_base)) fliplr(end_base)];
    end
    % Begin Point
    constraint = [derivatives_begin(1, :); zeros(3, size(derivatives_begin, 2))];
    constraint = reshape(constraint, 1, prod(size(constraint)));
    constraints = [constraint; 0 constraint(1 : end - 1); 0 0 constraint(1 : end - 2); 0 0 0 constraint(1 : end - 3)];
    constraints = [zeros(size(constraints, 1), size(constraints, 2) * (id - 1)) ...
        constraints, zeros(size(constraints, 1), size(constraints, 2) *(num_segments - id))];
    A = [A; constraints];
    b = [b; waypoints(id, :)'];
    % End Point
    constraint = [derivatives_end(1, :); zeros(3, size(derivatives_end, 2))];
    constraint = reshape(constraint, 1, prod(size(constraint)));
    constraints = [constraint; 0 constraint(1 : end - 1); 0 0 constraint(1 : end - 2); 0 0 0 constraint(1 : end - 3)];
    constraints = [zeros(size(constraints, 1), size(constraints, 2) * (id - 1)) ...
        constraints, zeros(size(constraints, 1), size(constraints, 2) *(num_segments - id))];
    A = [A; constraints];
    b = [b; waypoints(id + 1, :)'];
    % derivative continuity
    for i = 1 : 4
        constraint = [derivatives_begin(i + 1, :); zeros(3, size(derivatives_begin, 2))];
        constraint = reshape(constraint, 1, prod(size(constraint)));
        constraints = [constraint; 0 constraint(1 : end - 1); 0 0 constraint(1 : end - 2); 0 0 0 constraint(1 : end - 3)];
        if id == 1
            constraints = [zeros(size(constraints, 1), size(constraints, 2) * (id - 1)) ...
            constraints, zeros(size(constraints, 1), size(constraints, 2) *(num_segments - id))];
            A = [A; constraints];
            b = [b; 0; 0; 0; 0];
        else
            constraints = [zeros(size(constraints, 1), size(constraints, 2) * (id - 2)) ...
                constraints, -constraints, zeros(size(constraints, 1), size(constraints, 2) *(num_segments - id))];
            A = [A; constraints];
            b = [b; 0; 0; 0; 0];
        end
    end
end

% Find Solution
options = optimset ('Algorithm', 'active-set', 'largescale', 'off');
[basis_coef,fval,exitflag,output,lambda] = quadprog(H_all, [], [], [], A, b, [], [], [], options);
basis_coef = reshape(basis_coef, length(basis_coef) / m, m);

%%绘制图形
coef_0=[reshape(basis_coef(:,1),4,8);   % the first segment, four lines for x,z,beta,yaw. for ascending
        reshape(basis_coef(:,1),4,8);   % the second segment, four lines for x,z,beta,yaw. for ascending
        reshape(basis_coef(:,1),4,8);   % the third segment, four lines for x,z,beta,yaw. for ascending
        reshape(basis_coef(:,1),4,8)];  % the forth segment, four lines for x,z,beta,yaw. for ascending
coef_temp=coef_0;

coef_d1=[];     % First Derivative
coef_d2=[];     % Second Derivative
coef_d3=[];     % Third Derivative
coef_d4=[];     % Forth Derivative
for i=1:4  % four order derivation for every line
    coef_d=[];
    size(coef_0,1)
    for j=1:size(coef_0,1)     % for each line 
        
        d=polyder(fliplr(coef_temp(j,:)));
        coef_d=[coef_d;fliplr(d),zeros(1, size(coef_temp,2)-1 - length(d));];    % ascending
    end
    if i==1
        coef_d1=coef_d
    elseif i==2
        coef_d2=coef_d
    elseif i==3
        coef_d3=coef_d
    elseif i==4
        coef_d4=coef_d
        size(coef_d4,2)
    end
    coef_temp = coef_d;
    
end
% coef_1 = reshape(basis_coef(:,1),4,8)   % 
% d1x=polyder(fliplr(coef_1(1,:)));       % descending
% d1z=polyder(fliplr(coef_1(2,:)));
% d1b=polyder(fliplr(coef_1(3,:)));
% d1y=polyder(fliplr(coef_1(4,:)));
% coef_1d=[fliplr(d1x),zeros(1, n - length(d1x));fliplr(d1z),zeros(1, n - length(d1z));fliplr(d1b),zeros(1, n - length(d1b));fliplr(d1y),zeros(1, n - length(d1y))];% ascending
% coef_2 = reshape(basis_coef(:,2),4,8)
% d1x=polyder(fliplr(coef_2(1,:)));
% d1z=polyder(fliplr(coef_2(2,:)))
% d1b=polyder(fliplr(coef_2(3,:)))
% d1y=polyder(fliplr(coef_2(4,:)))
% coef_2d=[fliplr(d1x),zeros(1, n - length(d1x));fliplr(d1z),zeros(1, n - length(d1z));fliplr(d1b),zeros(1, n - length(d1b));fliplr(d1y),zeros(1, n - length(d1y))]
% coef_3 = reshape(basis_coef(:,3),4,8)
% d1x=polyder(fliplr(coef_3(1,:)));
% d1z=polyder(fliplr(coef_3(2,:)))
% d1b=polyder(fliplr(coef_3(3,:)))
% d1y=polyder(fliplr(coef_3(4,:)))
% coef_3d=[fliplr(d1x),zeros(1, n - length(d1x));fliplr(d1z),zeros(1, n - length(d1z));fliplr(d1b),zeros(1, n - length(d1b));fliplr(d1y),zeros(1, n - length(d1y))]
% coef_4 = reshape(basis_coef(:,4),4,8)
% d1x=polyder(fliplr(coef_4(1,:)));
% d1z=polyder(fliplr(coef_4(2,:)))
% d1b=polyder(fliplr(coef_4(3,:)))
% d1y=polyder(fliplr(coef_4(4,:)))
% coef_4d=[fliplr(d1x),zeros(1, n - length(d1x));fliplr(d1z),zeros(1, n - length(d1z));fliplr(d1b),zeros(1, n - length(d1b));fliplr(d1y),zeros(1, n - length(d1y))]
% 
% for id = 1 : num_segments
%     for  time_stamp = 0:0.01:timepoints(5);
%         t = [1; time_stamp; time_stamp.^2; time_stamp.^3; time_stamp.^4; time_stamp.^5; time_stamp.^6; time_stamp.^7]';%行向量
%         if 0 <= time_stamp <= timepoints(2)
%             coef = coef_1;
%             coef_d=coef_1d;
%         elseif  timepoints(2) < time_stamp <= timepoints(3)
%             coef = coef_2;
%             coef_d=coef_2d;
%         elseif timepoints(3) < time_stamp <= timepoints(4)
%             coef = coef_3;
%             coef_d=coef_3d;
%         elseif timepoints(2) < time_stamp <= timepoints(3)
%             coef = coef_4;
%             coef_d=coef_4d;
%         end
%         x = sum(coef(1, :).*t);
%         dx=sum(coef_d(1,:).*t(1:size(coef_d,2)));
%         z = sum(coef(2, :).*t);
%         dz=sum(coef_d(2,:).*t(1:size(coef_d,2)));
%         beta = sum(coef(3, :).*t);
%         figure(1); hold on; grid on;
%         plot(time_stamp,z);
%         xlabel('time (s)')
%         ylabel('z (m)');
%         figure(2);hold on; grid on;
%         plot(time_stamp,dz);
%         xlabel('time (s)')
%         ylabel('dz (m/s)');
% %         figure(3);hold on; grid on;
% %         plot(time_stamp,beta);
% %         xlabel('time (s)')
% %         ylabel('beta (degree)');
%     end
% end
hold off;



