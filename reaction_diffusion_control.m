%% Solution to Problem Set 3
% Control of Reaction Diffusion Equation 

%Start Parameter
lambda = 25;
control_on = true;
plot_temp = true;
plot_delta_t = 0.01; %must be >=0.002
plot_gain_kernel = false;
%End Parameter


delta_x = 0.01;
delta_t = 0.002;
x = 0:delta_x:1;
t = 0:delta_t:1;
u0 = sin(pi*x) + (0.5 - rand(1, length(x)) / 5);
u0(1) = 0;
u0(end) = 0;
A = full(gallery('tridiag',length(x),1,-2,1)) / delta_x^2;
A = A + eye(length(x)) * lambda;
A(1,:) = 0;
A(end, :) = 0;
A(1,1) = 1;
A(end,end) = 1;
u_be = zeros(length(x), length(t));
u_be(:,1) = u0;

k = zeros(1, length(x));
for i=1:length(x)
       k(1, i) = -lambda*x(1,i) * besseli(1, sqrt(lambda*(1-x(1,i)^2))) / sqrt(lambda*(1-x(1,i)^2));
       if i == length(x)
           k(1, i) = -lambda/2;
       end
end

for i=2:length(t)+1
    if control_on
        u_be(end, i-1) = trapz(delta_x, k(:).*u_be(:,i-1));
    else
        u_be(end, i-1)= 0;
    end
    if i == length(t)+1
        break
    end
    u_be(1, i-1) = 0;
    xi = eye(length(x)) - delta_t * A;
    u_be(:,i)= xi\u_be(:,i-1);
end

if plot_temp
    pdt = ceil(plot_delta_t / delta_t);
    for i=1:pdt:length(t)
        p = plot(x, u_be(:,i));
        ylim([min(u_be(:)) max(u_be(:))])
        title("Lambda=" + lambda + " and time=" + i*delta_t)
        pause(0.5);
        if ~ishghandle(p)
            return
        end
    end
end

close all;


%% Gain Kernel

if plot_gain_kernel
    delta_y = delta_x;
    y = x;
    lambda_arr = [5 10 15 20 25];
    k = zeros(length(lambda_arr), length(y));

    for l=1:length(lambda_arr)
        lambda = lambda_arr(l);
        for i=1:length(y)
           k(l, i) = -lambda*y(1,i) * besseli(1, sqrt(lambda*(1-y(1,i)^2))) / sqrt(lambda*(1-y(1,i)^2));
           if i == length(y)
               k(l, i) = -lambda/2;
           end
        end
        hold on
        plot(y, k(l, :), 'DisplayName', string(lambda))
    end

    hold off
    lgd = legend;
    lgd.Title.String = 'Lambdas';
end




    






