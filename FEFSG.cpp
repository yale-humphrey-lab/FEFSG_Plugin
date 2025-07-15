//#include "stdafx.h"
#include "FEFSG.h"
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FEAnalysis.h>                  // to get end time
#include <FECore/FEModel.h>                     // to get current time
#include <iostream>                             // to use cin.get()
#include <sstream>
#include <limits>
#include <fstream>
#include <algorithm>  // Include for std::count
#include <iterator>   // Include for std::istreambuf_iterator

// define the material parameters
BEGIN_FECORE_CLASS(FEFSG, FEUncoupledMaterial)
    ADD_PARAMETER(m_elastin_injury_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "elastin_injury_val");
    ADD_PARAMETER(m_crosslinking_injury_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "crosslinking_injury_val");
    ADD_PARAMETER(m_mechanosensing_injury_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "mechanosensing_injury_val");
    ADD_PARAMETER(m_mechanosensitivity_injury_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "mechanosensitivity_injury_val");
    ADD_PARAMETER(m_mechanoregulation_injury_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "mechanoregulation_injury_val");
    ADD_PARAMETER(e_r , "e_r");
    ADD_PARAMETER(e_t , "e_t");
    ADD_PARAMETER(e_z , "e_z");
END_FECORE_CLASS();

// define the material parameters
FEFSG::FEFSG(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_K = 0;    // invalid value!
    m_npmodel = 0;
    //m_secant_tangent = true;
    m_elastin_injury_val = 0;
    m_crosslinking_injury_val = 0;
    m_mechanosensing_injury_val = 0;
    m_mechanosensitivity_injury_val = 0;
    m_mechanoregulation_injury_val = 0;
    e_r = vec3d(1,0,0);
    e_t = vec3d(0,1,0);
    e_z = vec3d(0,0,1);
}

FEMaterialPointData* GRMaterialPoint::Copy()
{
    GRMaterialPoint* pt = new GRMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void GRMaterialPoint::Init()
{
    FEMaterialPointData::Init();

    m_dt = 1.0;
    sn = 0;
    K_delta_tauw = 0.0;
    K_delta_sigma = 0.0;
    sigma_inv_h = 0.0;
    sigma_inv_curr = 0.0;
    rho_hat_h = 0.0;
    CB = 0.0;
    CS = 0.0;
    bar_tauw_curr = 0.0;
    bar_tauw_h = 0.0;
    lambda_m = 0.0;
    lambda_0 = 0.0;
    lambda_m_act = 1.0;
    T_act = 0.0;
    k_act = 0.0;
    m_W = 0.0;
    m_sigma = mat3ds(0.0);
    m_F_curr = mat3d(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    m_J_curr = 1.0;
    m_CC = tens4ds(0.0);
    A_max = 0.0;


    for (int i=0; i<MAX_TIMESTEPS; ++i)
    {
        m_F_s[i] = mat3d(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        m_J_s[i] = 1.0;
        rhoR[i] = 0.0;
        rho[i] = 0.0;
        ups_infl_d[i] = 0.0;
        ups_infl_p[i] = 0.0;
    }


    // Hardcoded filename as a string variable
    std::string filename = "configuration.txt";

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Read non-vector values from the first line
    inputFile >> m_dt >> A_max >> rho_hat_h >> bar_tauw_h >> sigma_inv_h >> K_delta_tauw >> K_delta_sigma >> k_act >> lambda_0 >> lambda_m >> CB >> CS >> T_act;

    sigma_inv_curr = sigma_inv_h;
    bar_tauw_curr = bar_tauw_h;
    rhoR[0] = rho_hat_h;
    rho[0] = rho_hat_h;

    // Determine the number of lines in the file (excluding the first line)
    m_nconstituents = std::count(std::istreambuf_iterator<char>(inputFile), std::istreambuf_iterator<char>(), '\n');
    if (m_nconstituents > MAX_CONSTITUENTS) {
        std::cerr << "Too many constituents for simulation. Edit plugin if you want more." << filename << std::endl;
        return;
    }

    // Move the file pointer back to the beginning
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);

    // Skip the first line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Read values for each constituent
    for (int i=0; i<m_nconstituents; ++i) {
        inputFile >> m_constituents[i].m_degradable >> m_constituents[i].m_inflammatory >> m_constituents[i].m_active >> m_constituents[i].m_polymer
                  >> m_constituents[i].c1_alpha_h >> m_constituents[i].c2_alpha_h >> m_constituents[i].eta_alpha_h >> m_constituents[i].g_alpha_h
                  >> m_constituents[i].g_alpha_r >> m_constituents[i].g_alpha_theta >> m_constituents[i].g_alpha_z >> m_constituents[i].phi_alpha
                  >> m_constituents[i].k_alpha_h >> m_constituents[i].K_tauw_p_alpha_h >> m_constituents[i].K_sigma_p_alpha_h
                  >> m_constituents[i].K_tauw_d_alpha_h >> m_constituents[i].K_sigma_d_alpha_h;

        m_constituents[i].eta_alpha_h = m_constituents[i].eta_alpha_h * PI / 180.0;
        m_constituents[i].rho_hat_alpha_h = rho_hat_h;
        m_constituents[i].rhoR_alpha_h = m_constituents[i].phi_alpha * rho_hat_h;

        m_constituents[i].c1_alpha = m_constituents[i].c1_alpha_h;
        m_constituents[i].c2_alpha = m_constituents[i].c2_alpha_h;
        m_constituents[i].g_alpha = m_constituents[i].g_alpha_h;


        for (int j=0; j<MAX_TIMESTEPS; ++j)
        {
            m_constituents[i].mR_alpha[j] = m_constituents[i].k_alpha_h * m_constituents[i].rhoR_alpha_h;
            m_constituents[i].rhoR_alpha[j] = m_constituents[i].rhoR_alpha_h;
            m_constituents[i].epsilonR_alpha[j] = m_constituents[i].phi_alpha;
            m_constituents[i].epsilon_alpha[j] = m_constituents[i].phi_alpha;
            m_constituents[i].k_alpha[j] = m_constituents[i].k_alpha_h;
        }

    }


    // Close the file
    inputFile.close();

}

void GRMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);

    ar & sn;
    ar & m_dt;
    ar & m_nconstituents;
    ar & K_delta_tauw;
    ar & K_delta_sigma;
    ar & sigma_inv_h;
    ar & sigma_inv_curr;
    ar & rho_hat_h;
    ar & m_constituents;
    ar & m_F_s;
    ar & m_J_s;
    ar & rhoR;
    ar & rho;
    ar & ups_infl_p;
    ar & ups_infl_d;
    ar & CB;
    ar & CS;
    ar & bar_tauw_curr;
    ar & bar_tauw_h;
    ar & lambda_m;
    ar & lambda_0;
    ar & T_act;
    ar & k_act;
    ar & m_sigma;
    ar & m_F_curr;
    ar & m_J_curr;
    // ar & m_CC;
}

FEMaterialPointData* FEFSG::CreateMaterialPointData() 
{ 
    return new GRMaterialPoint(new FEElasticMaterialPoint); 
}

//! Update material point data.
void GRMaterialPoint::Update(const FETimeInfo& timeInfo)
{

    FEElasticMaterialPoint& et = *(this->ExtractData<FEElasticMaterialPoint>());

    // get current and end times
    const double t = timeInfo.currentTime;
    const double dt = timeInfo.timeIncrement;
    const int n_iter = timeInfo.currentRestart;

    double omega = double(n_iter)/20.0;
    if (omega > 1.0) {
        omega = 1.0;
    }

    sn = int(t) - 1; //int(sn + dt);

    if (sn > 0){

        sigma_inv_curr = et.m_s.tr();

        update_kinetics(sn);

        // Ramp Jacobian to prevent numerical issues
        et.m_J_star = et.m_J_star * (1.0 - omega) +  omega * m_J_s[sn];
        m_F_s[sn] = m_F_curr;

    } else {

        sigma_inv_curr = et.m_s.tr();
        sigma_inv_h = et.m_s.tr();
    }

    // don't forget to call the base class
    FEMaterialPointData::Update(timeInfo);

}


void FEFSG::DevStressTangent(FEMaterialPoint& mp, mat3ds& devstress, tens4ds& devtangent)
{       

    // The FEMaterialPoint classes are stored in a linked list. The specific material
    // point data needed by this function can be accessed using the ExtractData member.
    // In this case, we want to FEElasticMaterialPoint data since it stores the deformation
    // information that is needed to evaluate the stress.
    FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
    GRMaterialPoint& pt = *mp.ExtractData<GRMaterialPoint>();

    // We'll need the deformation gradient and its determinant in this function.
    // Note that we don't take the determinant of F directly (using mat3d::det)
    // but instead use the m_J member variable of FEMaterialPoint.
    const mat3d &F = et.m_F;
    const double J = et.m_J;
    const int sn = pt.sn;

    const vec3d X = mp.m_r0;
    pt.bar_tauw_curr = 1/(pow(sqrt(pow(mp.m_r0(0),2) + pow(mp.m_r0(1),2)),3));

    // Evaluate the distortional deformation gradient
    mat3d F_bar = F*pow(J / et.m_J_star, -1. / 3.);

    // push deformation gradient to local coordinates
    mat3d Q = mat3d(e_r(mp), e_t(mp), e_z(mp));

    if (sn > 0){
        pt.K_delta_sigma = m_mechanosensitivity_injury_val(mp);
        pt.m_constituents[0].c1_alpha = pt.m_constituents[0].c1_alpha_h*(1.0 - m_elastin_injury_val(mp));

        //Loop through each constituent to update its mass density
        for  (int alpha=2; alpha<pt.m_nconstituents; ++alpha) {
            pt.m_constituents[alpha].c1_alpha = pt.m_constituents[alpha].c1_alpha_h*(1.0 - m_crosslinking_injury_val(mp));
            pt.m_constituents[alpha].g_alpha = pt.m_constituents[alpha].g_alpha_h*(1.0 - m_mechanoregulation_injury_val(mp));
        }

    }

    pt.m_F_curr = Q.transpose() * F_bar * Q;
    pt.m_J_curr = J;

    // calculate the stress as a sum of deviatoric stress and pressure
    pt.update_sigma(sn);

    mat3ds sbar = (Q * pt.m_sigma * Q.transpose()).sym();

    devstress = sbar.dev();

    const mat3dd  I(1.0);
    const tens4ds IxI = dyad1s(I);
    const tens4ds IoI = dyad4s(I);

    tens4ds cbar = pt.m_CC.pp(Q);

    devtangent = cbar - 1./3.*(ddots(cbar, IxI) - IxI*(cbar.tr()/3.))
    + 2./3.*((IoI-IxI/3.)*sbar.tr()-dyad1s(sbar.dev(),I));


    et.m_v.x = pt.m_J_s[sn];                // Target volume 
    et.m_v.y = m_elastin_injury_val(mp);    // Aneurysm injuryy
    et.m_v.z = pt.m_W;

    et.m_a.x = pt.m_CC(0,0,0,0) + 2.0 *  pt.m_sigma(0,0);     // Radial stiffness    
    et.m_a.y = pt.m_CC(1,1,1,1) + 2.0 *  pt.m_sigma(1,1);     // Circumfrential stiffness    
    et.m_a.z = pt.m_CC(2,2,2,2) + 2.0 *  pt.m_sigma(2,2);     // Axial stiffness

}

void GRMaterialPoint::update_kinetics(int sn) {

    //Differences in current mechanical state from the reference state
    //Stress invariant
    double delta_sigma = (sigma_inv_curr / sigma_inv_h) - 1;

    //Wall shear stress
    //TODO implement shear stress
    double delta_tauw = 1;

    //Initialize pars for looping later
    double K_sigma_p = 0, K_tauw_p = 0, K_sigma_d = 0, K_tauw_d = 0;
    double upsilon_p = 0, upsilon_d = 0;

    double k_alpha_s = 0;
    double mR_alpha_s = 0;
    double rhoR_alpha_calc = 0;

    double mq_0 = 0, mq_1 = 0, mq_2;
    double q_0 = 0, q_1 = 0, q_2;
    double k_0 = 0, k_1 = 0, k_2;
    double rhoR_s = 0, rhoR_alpha_s = 0;
    double J_s = 0;

    int taun_min = 0;
    double s = m_dt * double(sn);
    double tau_max = A_max; //max fiber age
    //Determine if beyond initial time
    if (s > tau_max) {
        taun_min = sn - int(tau_max / m_dt);
    }
    else {
        taun_min = 0;
    }

    //Loop through each constituent to update its mass density
    for  (int alpha=0; alpha<m_nconstituents; ++alpha) {

        if (sn > 0 && m_constituents[alpha].m_degradable) {

            if (m_constituents[alpha].m_inflammatory == 0) {

                //Get the gains for the current constituent
                K_sigma_p = m_constituents[alpha].K_sigma_p_alpha_h;
                K_tauw_p = m_constituents[alpha].K_tauw_p_alpha_h;

                K_sigma_d = m_constituents[alpha].K_sigma_d_alpha_h;
                K_tauw_d = m_constituents[alpha].K_tauw_d_alpha_h;

                //Update the stimulus functions for each constituent
                upsilon_p = 1 + (1 - K_delta_sigma) * K_sigma_p * delta_sigma - K_tauw_p * delta_tauw;
                //if (upsilon_p < 0.1) {
                //    upsilon_p = 0.1;
                //}

                upsilon_d = 1 + K_sigma_d * pow(delta_sigma, 2) + K_tauw_d * pow(delta_tauw, 2);

                rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha[sn];

                /* FOR TEVG
                if (m_constituents[alpha].rhoR_alpha[sn] <= m_constituents[alpha].rhoR_alpha_h) {
                    rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha_h;
                }
                else {
                    rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha[sn];
                }
                */
            }
            else {

                upsilon_p = ups_infl_p[sn];
                upsilon_d = 1 + ups_infl_d[sn];

                rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha_h;

            }

            //Make sure productions don't become negative
            upsilon_p = (upsilon_p > 0)* upsilon_p;
            upsilon_d = (upsilon_d > 0)* upsilon_d;

            //Reset these to zero for each constituent
            k_alpha_s = 0;
            mR_alpha_s = 0;
            rhoR_alpha_s = 0;

            //Update kinetic values for the current time
            k_alpha_s = m_constituents[alpha].k_alpha_h * upsilon_d;
            
            mR_alpha_s = m_constituents[alpha].k_alpha_h * rhoR_alpha_calc * upsilon_p; //rhoR_alpha_h[alpha];
            m_constituents[alpha].k_alpha[sn] = k_alpha_s;
            m_constituents[alpha].mR_alpha[sn] = mR_alpha_s;

            k_2 = m_constituents[alpha].k_alpha[sn];
            q_2 = 1.0;
            mq_2 = m_constituents[alpha].mR_alpha[sn] * q_2;

            //loop through and update constituent densities from previous time points
            //starting from the current time point and counting down is more efficient
            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                //Trapazoidal rule     
                k_1 = m_constituents[alpha].k_alpha[taun];
                q_1 = exp(-(k_2 + k_1) * m_dt / 2) * q_2;
                mq_1 = m_constituents[alpha].mR_alpha[taun] * q_1;

                rhoR_alpha_s += (mq_2 + mq_1) * m_dt / 2;

                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

            }

            //Account for the cohort of material present initially
            if (taun_min == 0) {
                rhoR_alpha_s += m_constituents[alpha].rhoR_alpha[0] * q_1;
            }

            //Update referential volume fraction
            m_constituents[alpha].epsilonR_alpha[sn] = rhoR_alpha_s / m_constituents[alpha].rho_hat_alpha_h;
            J_s += m_constituents[alpha].epsilonR_alpha[sn];
        }
        else {

            //Precalculate maintenance or loss of constituents not produced
            rhoR_alpha_s = m_constituents[alpha].rhoR_alpha[sn];
            m_constituents[alpha].epsilonR_alpha[sn] = rhoR_alpha_s / m_constituents[alpha].rho_hat_alpha_h;
            J_s += m_constituents[alpha].epsilonR_alpha[sn];
        }

        m_constituents[alpha].rhoR_alpha[sn] = rhoR_alpha_s;
        rhoR_s += rhoR_alpha_s;

    }
    //Update spatial volume fractions
    for  (int alpha=0; alpha<m_nconstituents; ++alpha) {
        m_constituents[alpha].epsilon_alpha[sn] = m_constituents[alpha].epsilonR_alpha[sn] / J_s;
    }

    rhoR[sn] = rhoR_s;
    rho[sn] = rhoR_s / J_s;
    m_J_s[sn] = J_s;

}

void GRMaterialPoint::update_sigma(int sn) {

    //Find the current deformation gradient
    mat3d F_s = m_F_curr;
    mat3ds C_s = (F_s.transpose()*F_s).sym();
    mat3d F_tau;
    mat3ds C_tau;

    //Local active variables
    double lambda_act = 0;
    double lambda_act_1 = 0;
    double lambda_act_2 = 0;

    //Stiffness tensor
    tens4ds CC(0.0);
    tens4ds hat_CC_1(0.0);
    tens4ds hat_CC_2(0.0);

    //Stress tensors
    mat3ds sigma(0);
    mat3ds hat_sigma_1(0);
    mat3ds hat_sigma_2(0);

    //Energy doubles
    double W = 0.0;
    double hat_W_1 = 0.0;
    double hat_W_2 = 0.0;

    //For active stress
    double q_act_1 = 0, q_act_2 = 0;

    //Integration variables
    //For mass
    double mq_1 = 0, mq_2 = 0;
    double q_1 = 1.0, q_2 = 1.0;
    double k_1 = 0, k_2 = 0;

    int taun_min = 0;
    double s = m_dt * double(sn);
    double tau_max = A_max; //max fiber age
    //Determine if beyond initial time
    if (s > tau_max) {
        taun_min = sn - int(tau_max / m_dt);
    }
    else {
        taun_min = 0;
    }

    for  (int alpha=0; alpha<m_nconstituents; ++alpha) {

        //Trapz rule allows for fast heredity integral evaluation
        k_2 = m_constituents[alpha].k_alpha[sn];
        q_2 = 1.0;
        mq_2 = m_constituents[alpha].mR_alpha[sn];

        F_tau = m_F_s[sn];
        C_tau = (F_tau.transpose()*F_tau).sym();

        //Find active radius from current cohort
        if (m_constituents[alpha].m_active == 1) {
            lambda_act_2 = m_F_s[sn](1, 1);
            q_act_2 = 1.0;
        }

        m_constituents[alpha].constitutive(F_s, F_tau, sn, hat_W_2, hat_sigma_2, hat_CC_2);
        hat_W_2     = (1. / m_J_curr) * hat_W_2;
        hat_sigma_2 = (1. / m_J_curr) * hat_sigma_2;
        hat_CC_2    = (1. / m_J_curr) * hat_CC_2;

        //Check if during G&R or at initial time point
        if (sn > 0 && m_constituents[alpha].m_degradable) {

            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                //Find 1st intermediate kinetics
                k_1 = m_constituents[alpha].k_alpha[taun];
                q_1 = exp(-(k_2 + k_1) * m_dt / 2) * q_2;
                mq_1 = m_constituents[alpha].mR_alpha[taun] * q_1;

                F_tau = m_F_s[taun];
                C_tau = (F_tau.transpose()*F_tau).sym();

                //Find intermediate active state
                if (m_constituents[alpha].m_active == 1) {
                    lambda_act_1 = m_F_s[taun](1, 1);
                    q_act_1 = exp(-k_act * m_dt) * q_act_2;
                }

                m_constituents[alpha].constitutive(F_s, F_tau, sn, hat_W_1, hat_sigma_1, hat_CC_1);
                hat_W_1     = (1. / m_J_curr) * hat_W_1;
                hat_sigma_1 = (1. / m_J_curr) * hat_sigma_1;
                hat_CC_1    = (1. / m_J_curr) * hat_CC_1;

                // Add integration contribution
                W     += (mq_2 * hat_W_2 + mq_1 * hat_W_1) / m_constituents[alpha].rho_hat_alpha_h * m_dt / 2;
                sigma += (mq_2 * hat_sigma_2 + mq_1 * hat_sigma_1) / m_constituents[alpha].rho_hat_alpha_h * m_dt / 2;
                CC    += (mq_2 * hat_CC_2    + mq_1 * hat_CC_1)    / m_constituents[alpha].rho_hat_alpha_h * m_dt / 2;

                //Store active vars for next iteration
                //Find intermediate active state
                if (m_constituents[alpha].m_active == 1) {
                    lambda_act += k_act * (q_act_2 * lambda_act_2 + q_act_1 * lambda_act_1) * m_dt / 2;
                    lambda_act_2 = lambda_act_1;
                    q_act_2 = q_act_1;
                }

                //Store intermediate kinetics for next iteration
                k_2 = k_1;
                q_2 = q_1;
                mq_2 = mq_1;

                //Store intermediate stress and stiffness for next iteration
                hat_W_2     = hat_W_1;
                hat_sigma_2 = hat_sigma_1;
                hat_CC_2    = hat_CC_1;

            }

            // Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                W += m_constituents[alpha].rhoR_alpha[0]
                    / m_constituents[alpha].rho_hat_alpha_h * q_1 * hat_W_1;

                sigma += m_constituents[alpha].rhoR_alpha[0]
                    / m_constituents[alpha].rho_hat_alpha_h * q_1 * hat_sigma_1;

                CC += m_constituents[alpha].rhoR_alpha[0]
                    / m_constituents[alpha].rho_hat_alpha_h * q_1 * hat_CC_1;

                if (m_constituents[alpha].m_active == 1) {
                    if (!m_constituents[alpha].m_degradable){
                        q_act_1 = exp(-k_act * m_dt) * q_act_2;
                    }
                    lambda_act += m_F_s[0](1, 1) * q_act_1;
                }
            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            F_tau  = m_F_s[0];

            m_constituents[alpha].constitutive(F_s, F_tau, sn, hat_W_2, hat_sigma_2, hat_CC_2);
            hat_W_2     = (1. / m_J_curr) * hat_W_2;
            hat_sigma_2 = (1. / m_J_curr) * hat_sigma_2;
            hat_CC_2    = (1. / m_J_curr) * hat_CC_2;

            // Add in stress and stiffness contributions
            W += m_constituents[alpha].rhoR_alpha[sn] /
                m_constituents[alpha].rho_hat_alpha_h * hat_W_2;

            // Add in stress and stiffness contributions
            sigma += m_constituents[alpha].rhoR_alpha[sn] /
                m_constituents[alpha].rho_hat_alpha_h * hat_sigma_2;                

            CC += m_constituents[alpha].rhoR_alpha[sn] /
                m_constituents[alpha].rho_hat_alpha_h * hat_CC_2;

            lambda_act = lambda_act_2;
        }


        if (m_constituents[alpha].m_active == 1) {

            lambda_m_act = m_F_s[sn](1, 1)/lambda_act;

            // Copy the local element basis directions to n
            vec3d n[2];
            n[0].x = 0; n[0].y = 1; n[0].z = 0;
            n[1].x = 0; n[1].y = 0; n[1].z = 1;
            // Evaluate the structural direction in the current configuration
            vec3d ar,a;
            double eta_alpha_curr = m_constituents[alpha].eta_alpha_h;
            double cg = cos(eta_alpha_curr); double sg = sin(eta_alpha_curr);
            ar = n[0]*sg + n[1]*cg;
            ar = ar/(m_F_s[sn]*ar).norm();
            a = m_F_s[sn]*ar;

            mat3ds h0 = dyad(a);

            //Find active stress contribtion
            //add in initial active stress radius contribution
            double phismc = m_constituents[alpha].rhoR_alpha[sn] / m_constituents[alpha].rho_hat_alpha_h / m_J_curr;
            double C = CB - CS * (bar_tauw_curr / bar_tauw_h - 1);

            mat3ds  sigma_act_mat = h0 * phismc * T_act * (1 - exp(-pow(C, 2))) * lambda_m_act * (1.0-pow((lambda_m-lambda_m_act)/(lambda_m-lambda_0),2));
            //tens4ds CC_act_mat = dyad1s(h0) * phismc * T_act * (1 - exp(-pow(C, 2))) * lambda_m_act * (1.0-pow((lambda_m-lambda_m_act)/(lambda_m-lambda_0),2));

            sigma += sigma_act_mat;
            //CC += CC_act_mat;

        }

    }

    m_W = W;
    m_sigma = sigma;
    m_CC = CC;

    return;
}

