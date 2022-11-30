fn main() {
    println!("Hello, world!");
    let p = Parameters {
        EIR: 33.,
        ft: 0.,
        eta: 0.0001305,
        rho: 0.85,
        a0: 2920.,
        s2: 1.67,
        rA: 0.00512821,
        rT: 0.2,
        rD: 0.2,
        rU: 0.00906627,
        rP: 0.2,
        dE: 12.,
        tl: 12.5,
        cD: 0.0676909,
        cT: 0.0034482,
        cU: 0.006203,
        g_inf:	1.82425,
        d1: 0.160527,
        dd: 3650.,
        ID0: 1.577533,
        kd: 0.476614,
        ud: 9.44512,
        ad0: 8001.99,
        fd0: 0.007055,
        gd: 4.8183,
        aA: 0.757,
        aU: 0.186,
        b0: 0.590076,
        b1: 0.5,
        db: 3650.,
        IB0: 43.8787,
        kb: 2.15506,
        ub: 7.19919,
        phi0: 0.791666,
        phi1: 0.000737,
        dc: 10950.,
        IC0: 18.02366,
        kc: 2.36949,
        uc: 6.06349,
        PM: 0.774368,
        dm: 67.6952,
        tau: 10.,
        mu: 0.132,
        f: 0.33333333,
        Q0: 0.92,
    };
    let ages = Vec::from_iter((1..=100).map(f64::from));
    let nodes :Vec<f64> = vec![-4.8594628, -3.5818235, -2.4843258, -1.4659891, -0.4849357, 0.4849357, 1.4659891, 2.4843258, 3.5818235, 4.8594628];
    let weights :Vec<f64> = vec![4.310653e-06, 7.580709e-04, 1.911158e-02, 1.354837e-01, 3.446423e-01, 3.446423e-01, 1.354837e-01, 1.911158e-02, 7.580709e-04, 4.310653e-06];
    let pop = Population::new(&ages, &nodes, &weights);
    let solution = equilibrium(&p, &pop);
    print!("{:.5}", pfpr210(&solution));
}

struct Parameters {
    EIR: f64,
    ft: f64,

    // age, heterogeneity in exposure
    eta: f64,
    rho: f64,
    a0: f64,
    s2: f64,

    // rate of leaving infection states
    rA: f64,
    rT: f64,
    rD: f64,
    rU: f64,
    rP: f64,

    // human latent period and time lag from asexual parasites to infectiousness
    dE: f64,
    tl: f64,

    // infectiousness to mosquitoes
    cD: f64,
    cT: f64,
    cU: f64,
    g_inf: f64,

    // anti-parasite immunity
    d1: f64,
    dd: f64,
    ID0: f64,
    kd: f64,
    ud: f64,
    ad0: f64,
    fd0: f64,
    gd: f64,
    aA: f64,
    aU: f64,

    // anti-infection immunity
    b0: f64,
    b1: f64,
    db: f64,
    IB0: f64,
    kb: f64,
    ub: f64,

    // clinical immunity
    phi0: f64,
    phi1: f64,
    dc: f64,
    IC0: f64,
    kc: f64,
    uc: f64,
    PM: f64,
    dm: f64,

    // mosquito parameters
    tau: f64,
    mu: f64,
    f: f64,
    Q0: f64,
}

struct Population {
    ghnodes: Vec<f64>,
    ghweights: Vec<f64>,
    na: usize,
    nh: usize,
    age20: usize,

    // ageing
    age_days_midpoint: Vec<f64>,
    r: Vec<f64>,
}

impl Population {
    pub fn new(age: &Vec<f64>, nodes: &Vec<f64>, weights: &Vec<f64>) -> Self {

        //initialise ages and midpoints and 20 age group
        let na = age.len();
        let mut age_days :Vec<f64> = vec![0.;na];
        for i in 0..na {
            age_days[i] = age[i] * 365.;
        }
        let mut age_days_midpoint :Vec<f64> = vec![0.;na];
        for i in 0..na {
            if i == (na - 1) {
                age_days_midpoint[i] = age_days[i];
            } else {
                age_days_midpoint[i] = (age_days[i] + age_days[i+1]) * 0.5;
            }
        }

        // calculate age group of 20 year old
        let mut age20 = 0;
        for i in 0..na {
            if i < (na-1) && age_days_midpoint[i]<=20.*365. && age_days_midpoint[i+1]>20.*365. {
                age20 = i;
                break;
            }
        }

        // get rate of ageing in each group
        let mut r :Vec<f64> = vec![0.;na];
        for i in 0..na {
            r[i] = if i==(na-1) { 0. } else { 1./(age_days[i+1]-age_days[i]) };
        }

        let nh = nodes.len();

        Self {
            ghnodes: nodes.clone(),
            ghweights: weights.clone(),
            na: na,
            nh: nh,
            age20: age20,

            // ageing
            age_days_midpoint: age_days_midpoint,
            r: r
        }
    }
}

struct Solution {
    S: Vec<f64>,
    T: Vec<f64>,
    D: Vec<f64>,
    A: Vec<f64>,
    U: Vec<f64>,
    P: Vec<f64>,

    prop: Vec<f64>,
    pos_M: Vec<f64>,
    inc: Vec<f64>,

    ICA: Vec<f64>,
    ICM: Vec<f64>,

    FOIM: f64,
}

impl Solution {
    pub fn new(na: usize) -> Self {
        Self {
            S: vec![0.; na],
            T: vec![0.; na],
            D: vec![0.; na],
            A: vec![0.; na],
            U: vec![0.; na],
            P: vec![0.; na],
            prop: vec![0.; na],
            pos_M: vec![0.; na],
            inc: vec![0.; na],
            ICA: vec![0.; na],
            ICM: vec![0.; na],
            FOIM: 0.,
        }
    }
}

fn pfpr210(solution: &Solution) -> f64 {
    solution.pos_M[2..=10].iter().sum::<f64>() / solution.prop[2..=10].iter().sum::<f64>()
}

fn equilibrium(p: &Parameters, pop: &Population) -> Solution {
    //initialise some vectors
    let mut prop: Vec<f64> = vec![0.;pop.na];
    let mut psi: Vec<f64> = vec![0.;pop.na];
    let mut zeta: Vec<f64> = vec![0.;pop.nh];

    // calculate proportion and relative biting rate in each age group
    for i in 0..pop.na {
        if i == 0 {
            prop[i] = p.eta/(pop.r[i]+p.eta)
        } else {
            prop[i] = pop.r[i-1]*prop[i-1]/(pop.r[i]+p.eta);
        }
        psi[i] = 1. - p.rho*(-pop.age_days_midpoint[i]/p.a0).exp();
    }

    // calculate EIR scaling factor over Gaussian quadrature nodes
    for i in 0..pop.nh {
        zeta[i] = (-p.s2*0.5 + (p.s2).sqrt()*pop.ghnodes[i]).exp();
    }

    // loop through Gaussian quadrature nodes
    let mut solution = Solution::new(pop.na);
    for i in 0..pop.nh {

        // calculate immunity functions and onward infectiousness at equilibrium for 
        let mut IB = 0.;
        let mut IC = 0.;
        let mut ID = 0.;
        let mut FOI :Vec<f64> = vec![0.;pop.na];
        let mut ICA :Vec<f64> = vec![0.;pop.na];
        let mut ICM :Vec<f64> = vec![0.;pop.na];
        let mut phi:Vec<f64> = vec![0.;pop.na];
        let mut q :Vec<f64> = vec![0.;pop.na];
        let mut cA :Vec<f64> = vec![0.;pop.na];
        let mut S :Vec<f64> = vec![0.;pop.na];
        let mut T :Vec<f64> = vec![0.;pop.na];
        let mut D :Vec<f64> = vec![0.;pop.na];
        let mut A :Vec<f64> = vec![0.;pop.na];
        let mut U :Vec<f64> = vec![0.;pop.na];
        let mut P :Vec<f64> = vec![0.;pop.na];
        let mut pos_M :Vec<f64> = vec![0.;pop.na];
        let mut inc:Vec<f64> = vec![0.;pop.na];
        let mut inf:Vec<f64> = vec![0.;pop.na];

        // all age groups
        for j in 0..pop.na {

            // rate of ageing plus death
            let re = pop.r[j] + p.eta;

            // update pre-erythrocytic immunity IB
            let eps = zeta[i] * p.EIR/365. * psi[j];
            IB = (eps/(eps*p.ub+1.) + re*IB)/(1./p.db+re);

            // calculate probability of infection from pre-erythrocytic immunity IB via
            // Hill function
            let b = p.b0*(p.b1 + (1.-p.b1)/(1.+(IB/p.IB0).powf(p.kb)));

            // calculate force of infection (lambda)
            FOI[j] = b*eps;

            // update clinical immunity IC
            IC = (FOI[j]/(FOI[j]*p.uc+1.) + re*IC)/(1./p.dc + re);
            ICA[j] = IC;

            // update detection immunity ID
            ID = (FOI[j]/(FOI[j]*p.ud+1.) + re*ID)/(1./p.dd + re);

            // calculate probability that an asymptomatic infection (state A) will be
            // detected by microscopy
            let fd = 1. - (1.-p.fd0)/(1.+(pop.age_days_midpoint[j]/p.ad0).powf(p.gd));
            q[j] = p.d1 + (1.-p.d1)/(1.+(ID/p.ID0).powf(p.kd)*fd);

            // calculate onward infectiousness to mosquitoes
            cA[j] = p.cU + (p.cD-p.cU)*q[j].powf(p.g_inf);
        }

        // calculate maternal clinical immunity, assumed to be at birth a proportion
        // of the acquired immunity of a 20 year old
        let IM0 = ICA[pop.age20]*p.PM;
        for j in 0..pop.na {

            // rate of ageing plus death
            let re = pop.r[j] + p.eta;

            // maternal clinical immunity decays from birth
            let ICM_prev = if j==0 { IM0 } else { ICM[j-1] };
            ICM[j] = ICM_prev*re/(1./p.dm + re);
        }

        // calculate probability of acquiring clinical disease as a function of 
        // different immunity types
        for j in 0..pop.na {
            phi[j] = p.phi0*(p.phi1 + (1.-p.phi1)/(1.+((ICA[j]+ICM[j])/p.IC0).powf(p.kc)));
        }

        // calculate equilibrium solution of all model states. Again, see references
        // above for details
        for j in 0..pop.na {

            // rate of ageing plus death
            let re = pop.r[j] + p.eta;

            // calculate beta values
            let betaT = p.rT + re;
            let betaD = p.rD + re;
            let betaA = FOI[j]*phi[j] + p.rA + re;
            let betaU = FOI[j] + p.rU + re;
            let betaP = p.rP + re;

            // calculate a and b values
            let aT = p.ft*phi[j]*FOI[j]/betaT;
            let bT = if j == 0  { 0. } else { pop.r[j-1]*T[j-1]/betaT };
            let aD = (1.-p.ft)*phi[j]*FOI[j]/betaD;
            let bD = if j == 0 { 0. } else { pop.r[j-1]*D[j-1]/betaD };
            let aP = p.rT*aT/betaP;
            let bP = (p.rT*bT + (if j == 0 { 0. } else { pop.r[j-1]*P[j-1] }))/betaP;

            // calculate Y
            let Y = (prop[j] - (bT+bD+bP))/(1.+aT+aD+aP);

            // calculate final {T,D,P} solution
            T[j] = aT*Y+bT;
            D[j] = aD*Y+bD;
            P[j] = aP*Y+bP;

            // calculate final {A, U, S} solution
            let mut rA = 0.;
            let mut rU = 0.;
            if j > 0 {
                rA = pop.r[j-1]*A[j-1];
                rU = pop.r[j-1]*U[j-1];
            }
            A[j] = (rA + (1.-phi[j])*Y*FOI[j] + p.rD*D[j])/(betaA + (1.-phi[j])*FOI[j]);
            U[j] = (rU + p.rA*A[j])/betaU;
            S[j] = Y - A[j] - U[j];

            // calculate proportion detectable by mocroscopy and PCR
            pos_M[j] = D[j] + T[j] + A[j]*q[j];

            // calculate clinical incidence
            inc[j] = Y*FOI[j]*phi[j];

            // calculate incidence of infection
            inf[j] = p.cD*D[j] + p.cT*T[j] + cA[j]*A[j] + p.cU*U[j];
        }

        // add to sums over nodes
        let mut delta_foim = 0.;
        for j in 0..pop.na {
            solution.S[j] += pop.ghweights[i]*S[j];
            solution.T[j] += pop.ghweights[i]*T[j];
            solution.D[j] += pop.ghweights[i]*D[j];
            solution.A[j] += pop.ghweights[i]*A[j];
            solution.U[j] += pop.ghweights[i]*U[j];
            solution.P[j] += pop.ghweights[i]*P[j];

            solution.pos_M[j] += pop.ghweights[i]*pos_M[j];
            solution.inc[j] += pop.ghweights[i]*inc[j];

            delta_foim += inf[j]*psi[j];
        }
        solution.FOIM += delta_foim*pop.ghweights[i]*zeta[i];

    }  // end loop over nodes

    // complete overall force of infection on mosquitoes
    solution.FOIM *= p.f*p.Q0/(1. - p.rho*p.eta/(p.eta+1./p.a0));
    solution.prop = prop;

    solution
}
