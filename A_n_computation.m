function  [A0, An] = A_n_computation(m, p, n) 
% n corresponding to the number of term in the Fourier transform
% An is a vector that contains the An values : An = [A1, A2, ..., An]
% Only A0 is computed separatly

% Derivative function for a 4 digit NACA profile :
    dyc_front = @(Theta) (2*m/(p^2)) .* ( p -( (1-cos(Theta)) /2 ) );
    dyc_back = @(Theta) (2*m/((1-p)^2)) .* (p-( (1-cos(Theta)) /2 ));

% computation of the angle at which the fonction that describes the profile
% change :
    Theta_p = acos(1-2.*p); 

% computation of A0
    A0 = (1/pi)*(integral(dyc_front,0,Theta_p) + integral(dyc_back,Theta_p,pi));

% computation of An 
% initialisation
An = zeros(1,n) ;
dyc_calculated = A0;
    for i=1:n 
    % Derivative function for a 4 digit NACA profile :
        dyc_front_n = @(Theta) (2*m/(p^2)) .* ( p -( (1-cos(Theta)) /2 ) ) .* cos(i*Theta) ;
        dyc_back_n = @(Theta) (2*m/((1-p)^2)) .* (p-( (1-cos(Theta)) /2 )) .* cos(i*Theta) ;
    % computation of An
        An(i) = (2/pi)*(integral(dyc_front_n,0,Theta_p) + integral(dyc_back_n,Theta_p,pi)) ;
        dyc_calculated = dyc_calculated + An(i)*cos(i*(pi/2));
    end

% from equation slide 20 of lecture notes 7.
% all values should be the same 
dyc_front(pi/2)
dyc_back(pi/2)
dyc_calculated

end 