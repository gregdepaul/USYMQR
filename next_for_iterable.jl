function next(m::USYMQRIterable, iteration::Int)

    m.tau_prev = m.tau;
    m.A
    m.p = p;
    m.q
    m.sigma
    m.tau
    m.sbar

    m.x
    m.resnorm

    m.q1
    m.q2

end


function updateSoln!(x, rh_s1, v_prev, tau_prev, sigma, rho, w1, w2)
    # Update solution
    w3 = ( v_prev - sigma*w2 - tau_prev*w1 )/rho;
    x = x + rh_s1*w3;
    w1 = w2;
    w2 = w3;
end
