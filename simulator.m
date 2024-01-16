t = linspace(0,200000,20001); %Gives the delta t values

moon_pos = @(x) [222*cos(2.6615*10^(-6)*x),222*sin(2.6615*10^(-6)*x)];
% Define moon pos as a function  of x
'''moon_pos = @(x) [0,222];''' % Alternative moon pos for stationary moon

theta = 52.4 % Launch angle
rocket_pos = [0,3.7] % Initial rocket position

[tout, pos] = simulate_rocket(rocket_pos,[0.0066*cosd(theta),0.0066*sind(theta)],moon_pos,t); 
% Calls function that will simulate the trajectory

plot(pos(:,1),pos(:,2)) %Plot the trajectory
hold on
plot(cos(pi*t/100000),sin(pi*t/100000)+222)
plot(3.7*cos(pi*t/100000),3.7*sin(pi*t/100000)) 
% These plot the earth and moon to help with visualisation

function[a] = accel(posr,posm)
% This is just the acceleration term calculated using Newton's law of gravitation

% Input:
% * posr: The position of the rocket at that time
% * posm: The position of the moon at that time

% Output:
% * a: The accleration felt by the rocket at that time

    a = (9.63*10^(-7))*(-83.3*posr/norm(posr)^3-(posr-posm)/norm(posr-posm)^3);
end

function[rout,vout] = loop(r,v,posm)
% This is the inner workings of the improved Euler method. Takes the
% previous pos, calculates an intermediate pos, uses that to find the new
% velocity then uses that again to find a new pos

% Input:
% * r: This is the position from the previous iteration
% * v: This is the velocity from the previous iteration
% * posm: This is the position of the moon at that specific time

% Output:
% * rout: The new position for this time
% * vout: The new velocity for this time
    rprime = r + 10*v;
    vout = v + 5*(accel(r,posm)+accel(rprime,posm)); 
    % Uses the function accel just to make code more digestible
    rout = r+5*(v+vout);
end

function [tout, pos] = simulate_rocket(init_pos, init_vel, moon_pos, t)
% Author: Joel Adams , Date: 21/3/2022
% Simulates the movement of the rocket under the Earth's and Moon's gravitational field.

% Input:
% * init_pos: 2−elements vector (x, y) indicating the initial position of the rocket. 
% * init_vel: 2−elements vector (vx, vy) of the initial velocity of the rocket.
% * moon_pos: a function that receives time, t, and return a 2−elements vector (x, y) 
% * t: an N−elements vector of the time step where the position of the rocket will be returned.

% Output:
% * tout: an M−elements vector that has all the time intervals that the rocket is in trajectory.
% * pos: (Mx2) matrix that gives (x,y) coord of rocket for each time interval

    pos = [init_pos];
    vel = init_vel;
    for n = 1:length(t)
        if norm(moon_pos(t(n))-pos(n,:)) < 1 
            break % Break out clause for if the rocket reaches the moon
        end
        if norm(pos(n,:)) < 3.7
            break % Break out clause for if the rocket returns to Earth
        end
        [pos(n+1,:),vel] = loop(pos(n,:),vel,moon_pos(t(n))); 
        % For each time interval calls the function loop with the previous
        % variables and saves the next iteration's variables.
    end
    tout = t;
end