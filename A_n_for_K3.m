function  [A_0] = A_n_for_K3(m, p) % n corresponding to the number of term in the Fourier transform

    % Then the derivative of the profile :
    dyc_front = @(Theta) (2*m/(p^2)) .* ( p -( (1-cos(Theta)) /2 ) );
    dyc_back = @(Theta) (2*m/((1-p)^2)) .* (p-( (1-cos(Theta)) /2 ));

    Theta_p = acos(1-2.*p);

    A_0 = (1/pi)*(integral(dyc_front,0,Theta_p) + integral(dyc_back,Theta_p,pi));

end 