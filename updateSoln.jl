function updateSoln(v_prev,tow_prev,sigma,rho,w1,w2)

    # Update solution
    w3 = ( v_prev - sigma*w2 - tow_prev*w1 )/rho;

    x = x + rh_s1*w3;
    w1 = w2;
    w2 = w3;

    return Dict('x' => x, 'w1' => w1, 'w2' => w2, 'w3' => w3);
end
