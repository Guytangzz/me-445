function  [test] = debug_function(dyc_front, dyc_back, Theta_p, A0, An) 


Theta_front = linspace(0,Theta_p,100);
Theta_back = linspace(Theta_p,pi,100);

figure;
plot(Theta_front,dyc_front(Theta_front))
hold on 
plot(Theta_back,dyc_back(Theta_back))
hold off
title('Dyc front and Dyc back')

sprintf('Theta_p = %d',Theta_p)


% from equation slide 20 of lecture notes 7.
% all values should be the same 
sprintf('dyc_front = %d',dyc_front(Theta_p));
sprintf('dyc_back = %d',dyc_back(Theta_p));
%dyc_calculated

test = 0 ;
end