function QRfact(iter,dictVars)
    # dictVars has all of the following vars for func inputs below
    #INITIALIZATION OUTSIDE FUNC-back in tridiag looop

    #unpack vars here
    alpha = dictVars['alpha'];
    beta = dictVars['beta'];
    gamma = dictVars['gamma'];
    gamma_prev = dictVars['gamma_prev'];
    q1 = dictVars['q1'];
    q2 = dictVars['q2'];
    AA_norm = dictVars['AA_norm'];
    rnorm = dictVars['rnorm'];
    atol = dictVars['atol'];

    #compute scalars
    A_rnorm = rnorm*norm([gamma_prev*q1 + alpha*q2; gamma*q2]);
    Anorm = sqrt(AA_norm);
    AA_norm = AA_norm + alpha^2 + beta^2 + gamma^2;

    if (Arnorm/(Anorm*rnorm) < atol)
        stop = 2;
    end

    iter =  dictVars['iter'];
    c = dictVars['c'];
    s = dictVars['s'];
    s_bar = dictVars['s_bar'];

    # plane rotation
    planeDict = planeRotation(iter,c,s,s_bar,alpha,beta,gamma);

    v_prev = dictVars['v_prev'];
    tow_prev = dictVars['tow_prev'];

    # Update solution
    updateDict = updateSoln(v_prev,tow_prev,planeDict['sigma'],planeDict['rho'],planeDict['w1'],planeDict['w2']);

    #CHECK CONVERGENCE OUTSIDE-back in tridiag loop
    return Dict('gamma' = > gamma,'rh_s2' = > planeDict['rh_s2'],'s' = > planeDict['s'],'c' = > planeDict['c'],'t' = > planeDict['t']);
end
