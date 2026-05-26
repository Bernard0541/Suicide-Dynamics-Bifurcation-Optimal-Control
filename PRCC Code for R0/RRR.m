function R=RRR(LHSmatrix,x)

Parameter_settings_LHS;
beta=LHSmatrix(x,1);
mu=LHSmatrix(x,2);
gamma=LHSmatrix(x,3);
alpha=LHSmatrix(x,4);
psi=LHSmatrix(x,5);
rho=LHSmatrix(x,6);
phi=LHSmatrix(x,7);
theta=LHSmatrix(x,8);
sigma=LHSmatrix(x,9);
omega=LHSmatrix(x,10);

R1=(beta*gamma*((theta+mu) + (1-rho)*alpha)*(sigma+mu)*(omega+mu))/((alpha+psi+mu)*(theta+mu)*(gamma+mu)*(sigma+mu)*(omega+mu)...
    - gamma*(1-phi)*omega*sigma*(psi*(theta+mu) + theta*(1-rho)*alpha));


R=(R1);
end