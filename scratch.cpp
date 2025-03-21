
public:
    //! total Cauchy stress (do not overload!)
    mat3ds Stress(FEMaterialPoint& mp) override final;

    //! total spatial tangent (do not overload!)
    tens4ds Tangent(FEMaterialPoint& mp) override final;


//-----------------------------------------------------------------------------
//! The stress function calculates the total Cauchy stress as a sum of 
//! two terms, namely the deviatoric stress and the pressure. 
mat3ds FEFSG::Stress(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    GRMaterialPoint& gt = *mp.ExtractData<GRMaterialPoint>();

    // get current and end times
    const double t = GetFEModel()->GetTime().currentTime;
    //TODO: Get and set dt better time step?
    const double dt = m_dt;
    const int sn = int(t/dt);

    mat3ds stress =  mat3dd(pt.m_p) + DevStress(mp);

    // calculate the stress as a sum of deviatoric stress and pressure
    pt.m_p = UJ(pt.m_J, gt.m_J_s[sn]);

    gt.sigma_inv_curr = stress.tr();

    printf("sigma_inv_curr: %f\n", gt.sigma_inv_curr);
    printf("pt.m_J: %f\n", pt.m_J);
    printf("gt.m_J_s[sn]: %f\n", gt.m_J_s[sn]);
    printf("pt.m_p: %f\n", pt.m_p);

    return mat3dd(pt.m_p) + DevStress(mp);
}

//------------------------------------------------------------------------------
//! The tangent function calculates the total spatial tangent, that is it calculates
//! the push-forward of the derivative of the 2ndPK stress with respect to C. However,
//! for an uncoupled material, the 2ndPK stress decouples in a deviatoric and a 
//! dilatational component. The deviatoric tangent is provided by the particular
//! material and the dilatational component is added here.
//!
tens4ds FEFSG::Tangent(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    GRMaterialPoint& gt = *mp.ExtractData<GRMaterialPoint>();

    // get current and end times
    const double t = GetFEModel()->GetTime().currentTime;
    //TODO: Get and set dt better time step?
    const double dt = m_dt;
    const int sn = int(t/dt);

    // 2nd-order identity tensor
    mat3dd I(1);

    // 4th-order identity tensors
    tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    
    // pressure
    pt.m_p = UJ(pt.m_J, gt.m_J_s[sn]);
    
    // tangent is sum of three terms
    // C = c_tilde + c_pressure + c_k
    //
    // + c_tilde is the derivative of the deviatoric stress with respect to C
    // + c_pressure is p*d(JC)/dC
    // + c_k comes from the derivative of p with respect to C
    // 
    // Note that the c_k term is not necessary in the 3F formulation (since p is independant variable) 
    // but we do need to add it here.
    //
    //        c_tilde         c_pressure            c_k
    return DevTangent(mp) + (IxI - I4*2)*pt.m_p + IxI*(UJJ(pt.m_J, gt.m_J_s[sn])*pt.m_J);
}




public:
    // TODO: removing virtual from the following 3 functions causes changes
    // to the convergence criteria on macOS, despite these functions not
    // being overridden anywhere.
    //! strain energy density U(J)
    virtual double U(double J, double J_tar) {
        switch (m_npmodel) {
            case 0: return 0.5*m_K*pow(log(J/J_tar),2); break;    // FEBio default
            case 1: return 0.25*m_K*(J/J_tar*J/J_tar - 2.0*log(J/J_tar) - 1.0); break;    // NIKE3D's Ogden material
            case 2: return 0.5*m_K*(J/J_tar-1)*(J/J_tar-1); break;      // ABAQUS
            case 3: return 0.5*m_K*((J/J_tar*J/J_tar-1)/2-log(J/J_tar)); break;      // ABAQUS - GOH
            default: { assert(false); return 0; }
        }
    }

    //! pressure, i.e. first derivative of U(J)
    virtual double UJ(double J, double J_tar) {
        printf("UJ");
        switch (m_npmodel) {
            case 0: return m_K*log(J/J_tar)/J; break;
            case 1: return 0.5*m_K*(J/(J_tar*J_tar) - 1.0/J); break;
            case 2: return m_K*(J/J_tar-1); break;
            case 3: return 0.5*m_K*(J/(J_tar*J_tar)-1.0/J); break;
            default: { assert(false); return 0; }
        }
    }

    //! second derivative of U(J) 
    virtual double UJJ(double J, double J_tar) {
        switch (m_npmodel) {
            case 0: return m_K*(1-log(J/J_tar))/(J*J); break;
            case 1: return 0.5*m_K*(1/(J_tar*J_tar) + 1.0/(J*J)); break;
            case 2: return m_K/J_tar; break;
            case 3: return 0.5*m_K*(1/(J_tar*J_tar) + 1.0/(J*J)); break;
            default: { assert(false); return 0; }
        }
    }