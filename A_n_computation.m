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
    % A0_cours = (m/(pi*p^2)) * ((2*p-1)*Theta_p + sin(Theta_p)) + (m/(pi*(1-p)^2))*((2*p-1)*(pi-Theta_p)-sin(Theta_p)) ;
    % fprintf("A0 = %d \n",A0)
    % fprintf("A0_cours = %d \n",A0_cours)
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
        % dyc_calculated = dyc_calculated + An(i)*cos(i*(pi/2));
    end

    % A1_cours = (2*m/(pi*p^2)) * ((2*p-1)*sin(Theta_p) + 1/4 *sin(2*Theta_p) + Theta_p/2) - (2*m/(pi*(1-p)^2)) * ((2*p-1)*sin(Theta_p) + 1/4 *sin(2*Theta_p) - 1/2 *(pi-Theta_p)) ;
    % fprintf("A1 = %d \n",An(1))
    % fprintf("A1_cours = %d \n",A1_cours)
    % 
    % A2_cours = (2*m/(pi*p^2)) * ((2*p-1)*(1/2)*sin(2*Theta_p) + sin(Theta_p) - (2/3)*sin(Theta_p)^3) - (2*m/(pi*(1-p)^2)) * ((2*p-1)*(1/2)*sin(2*Theta_p) + sin(Theta_p) - (2/3)*sin(Theta_p)^3) ;
    % fprintf("A2 = %d \n",An(2))
    % fprintf("A2_cours = %d \n",A2_cours)
% modify input and output 
%[test] = debug_function(dyc_front, dyc_back, Theta_p, A0, An)

end 