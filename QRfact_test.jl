include("QRfact.jl")
m = 3;
iter = 5;
dictVars = Dict{String,Any}();
dictVars["alpha"] = 0.2;
dictVars["beta"] = 0.5;
dictVars["gamma"] = 0.3;
dictVars["gamma_prev"] = 0.1;
dictVars["q1"] = 2.2;
dictVars["q2"] = 3.1;
dictVars["AA_norm"] = .98;
dictVars["rnorm"] = 0.67;
dictVars["atol"] = 10e-4;
dictVars["c"] = 0.56;
dictVars["s"] = 0.78;
dictVars["s_bar"] = 0.89;
dictVars["rh_s1"] = 0.12;
dictVars["v_prev"] = randn(m ,1);
dictVars["tow_prev"] = 0.43;
dictVars["w1"] = randn(m,1);
dictVars["w2"] = randn(m,1);

x=randn(m,1);
QRdict = QRfact(x,iter,dictVars);
