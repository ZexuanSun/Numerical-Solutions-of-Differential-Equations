% -------------Burgers program-------------

clear all;
%Selection of method.
method = menu('Choose a numerical method:', ...
'Lax-Friedrichs','Lax-Wendroff','Godunov','Beam-Warming');
xstart =-2;
xend = 2; % x-axis size.
N = 100;
dx = xend/N;
dt = 0.01;
tend = input('Enter the end time tend(better less than 1s,since I set xend = 2):');
choice = menu('Plot error or not', ...
'error','not');
x = xstart:dx:xend;
nt = floor(tend/dt);
dt = tend / nt;
%Set up the initial solution values.
uL = 2;
uR = 1;
s = (uL + uR )/2;
xshift = 0;
u0 = uR + (uL-uR) * ((x-xshift) <= 0.0); %Call to the function "uinit".
u = u0;
unew = 0*u;
%Implementation of the numerical methods.
%method = 4;
for i = 1 : nt 
switch method

case 1 %Lax-Friedrichs
    unew(2:end-1) = 0.5*(u(3:end)+u(1:end-2)) - 0.5*dt/dx * ...
    (f(u(3:end)) - f(u(1:end-2)));
    unew(1) = u(1);
    unew(end) = u(end);
case 2 %Lax-Wendroff
    % 2 到end-1 位
    unew(2:end-1) = u(2:end-1) ...
        - 0.5*dt/dx * (f(u(3:end)) - f(u(1:end-2))) ...
        + 0.5*(dt/dx)^2 * ...
    ( df(0.5*(u(3:end) + u(2:end-1))) .* (f(u(3:end)) - f(u(2:end-1))) - ...
      df(0.5*(u(2:end-1) + u(1:end-2))) .* (f(u(2:end-1)) - f(u(1:end-2))) );
    unew(1) = u(1);
    unew(end) = u(end);
case 3 %Godunov
    unew(2:end-1) =u(2:end-1)- dt/dx*(nf(u(2:end-1),u(3:end)) - nf(u(1:end-2),u(2:end-1)));
    unew(1) = u(1);
    unew(end) = u(end);
case 4 %Beam-Warming
    J = length(u);
    for j = 3:J - 2
        if u(j) >=  0
          unew(j) = u(j) - dt/(2*dx)*(3*f(u(j)) - 4*f(u(j-1)) + f(u(j-2))) +... 
          dt^2/(2*dx^2)*((u(j-1) + u(j))/2*(f(u(j)) - f(u(j-1))) - (u(j-2) +...
          u(j-1))/2*(f(u(j-1)) - f(u(j-2))) );
        elseif u(j) == 0
            unew(j) = 0;
        else u(j) < 0
          unew(j) = u(j) + dt/(2*dx)*(3*f(u(j)) - 4*f(u(j+1)) + f(u(j+2)) ) +...
              dt^2/(2*dx^2)*((u(j+1) + u(j+2))/2*(f(u(j+2)) -f(u(j+1))) -...
              (u(j) + u(j+1))/2*(f(u(j+1)) - f(u(j)) ));
        end
    end
    unew(2) = uL;
    unew(1) = uL;
%     unew(2) = unew(3);
%     unew(1) = unew(2);
    unew(J-1) = unew(J-2);
    unew(J) = u(J-1);
end
u = unew;
U(i,:) = u(:);
end
U=[u0;U];
T=0:dt:tend;
%Plot of the  real solutions.
if choice == 2
    ul =uL;
    ur = uR;
    t = tend;
    s = (ur + ul)/2;
    xx1 = xstart:0.01:s*t;
    xx2 = s*t+0.01:0.01:x(end);
    xx = [xx1 xx2];
    ans1 = ul * ones(1,length(xx1));
    ans2 = ur * ones(1,length(xx2));
    ans = [ans1,ans2];
    plot(xx,ans,'-','LineWidth',2)

    hold on 
    plot(x,u,'*')
    %xlim = ([x(1) x(end)]);
    ylim([0.5 2.5])
    switch method
        case 1
            title('Lax-Friedrichs')
        case 2
            title('Lax-Wendroff')
        case 3
            title('Godunov')
        case 4
            title('Beam-Warming')
    end
     toy = ['End time' num2str(tend) 's'];
     xlabel(toy)
      xlim([-2 2])
      ylim([uR-0.5 uL+0.5 ])
      legend('Exact. Sol.','Num. Sol.');
else
    error = 0*x;
    shock = s* tend;
    for i = 1:length(x)
        if(x(i) < shock)
            error(i) = abs(u(i) - uL);
        elseif (x(i) > shock)
            error(i) = abs(u(i)-uR);
        end
    end
    plot(x,error,'-o','LineWidth',2)
    switch method
        case 1
            title('Lax-Friedrichs error')
        case 2
            title('Lax-Wendroff error')
        case 3
            title('Godunov error')
        case 4
            title('Beam-Warming error')
    end
    toy = ['End time' num2str(tend) 's'];
     xlabel(toy);
     [m,p] = max(error);
     disp('Maximum Error')
     m
     disp('Index of Maximum Error')
     x(p)
     disp('Average Absolute Error')
     sum(error/length(error))
end   

% toy = ['End time' num2str(tend)];
% xlabel(toy)


%necessary functions
function ret = f( u )
    ret = 0.5 * u.^2;
end

function ret = df( u )
    ret = u;
end
function ret = nf( u, v )
    for i = 1:length(u)
        if (u(i) >= v(i))
            if ((u(i)+v(i))/2 > 0)
                ustar(i)=u(i);
            else
            ustar(i)=v(i);
            end
        else
            if (u(i)>0)
                ustar(i)=u(i);
            elseif (v(i)<0)
                ustar(i)=v(i);
            else
                ustar(i)=0;
            end
        end
    end
    ret =f(ustar);
end
% 
% figure(1)
% surf(x,T,U)
% shading interp
% xlabel('x'), ylabel('t'), zlabel ('u(x,t)');
% grid on
% colormap('Gray');
