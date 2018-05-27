function planeRotation(iter,c,s,s_bar,alpha,beta,gamma)

    #plane rotation vars
    sigma = c*s_bar + s*alpha;
    r_bar = - s*s_bar + c*alpha;
    tow = s*gamma;
    s_bar = c*gamma;

    if(iter == 1 )
        r_bar = alpha;
        s_bar = gamma;
    end

    rho = sqrt(r_bar^2 + beta^2);
    c = r_bar/rho;
    s = beta/rho;

    t = rh_s1;
    rh_s1 =  c*t;
    rh_s2 = -s*t;

    return Dict('s'=>s,'t'=>t,'sigma'=>sigma,'tow' => tow,'rho'=>rho,'rh_s1'=>rh_s1,'rh_s2'=>rh_s2)

end
