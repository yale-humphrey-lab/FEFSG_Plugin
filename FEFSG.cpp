//#include "stdafx.h"
#include "FEFSG.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FECore/FEAnalysis.h"                  // to get end time
#include "FECore/FEModel.h"                     // to get current time
#include "FECore/log.h"                         // to print to log file and/or screen
#include <iostream>                             // to use cin.get()
#include <sstream>
#include <signal.h>
#define _USE_MATH_DEFINES                       // to introduce pi constant (1/2)
#include <math.h>                               // to introduce pi constant (2/2)
#include <limits>
#include <fstream>
#include <algorithm>  // Include for std::count
#include <iterator>   // Include for std::istreambuf_iterator

// define the material parameters
BEGIN_FECORE_CLASS(FEFSG, FEUncoupledMaterial)
    ADD_PARAMETER(m_a_val      , FE_RANGE_GREATER_OR_EQUAL(0.0), "a_val");
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
    m_a_val = 0;
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
    nts = 720;
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
    T_act = 0.0;
    k_act = 0.0;
    m_sigma = mat3ds(0.0);
    m_CC = tens4ds(0.0);


    for (int i=0; i<MAX_TIMESTEPS; ++i)
    {
        m_lambda_act[i] = 0.0;
        m_F_s[i] = mat3d(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        m_J_s[i] = 1.0;
        rhoR[i] = 0.0;
        rho[i] = 0.0;
        ups_infl_d[i] = 0.0;
        ups_infl_p[i] = 0.0;
    }


    // Hardcoded filename as a string variable
    std::string filename = "example.txt";

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Read non-vector values from the first line
    inputFile >> m_dt >> rho_hat_h >> bar_tauw_h >> sigma_inv_h >> K_delta_tauw >> K_delta_sigma >> k_act >> lambda_0 >> lambda_m >> CB >> CS >> T_act;

    sigma_inv_curr = sigma_inv_h;
    bar_tauw_curr = bar_tauw_h;

    // Determine the number of lines in the file (excluding the first line)
    m_nconstituents = std::count(std::istreambuf_iterator<char>(inputFile), std::istreambuf_iterator<char>(), '\n');
    if (m_nconstituents > MAX_CONSTITUENTS) {
        std::cerr << "Too many constituents for simulation, edit plugin." << filename << std::endl;
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

        m_constituents[i].eta_alpha_h = m_constituents[i].eta_alpha_h * M_PI / 180.0;
        m_constituents[i].rho_hat_alpha_h = rho_hat_h;
        m_constituents[i].rhoR_alpha_h = m_constituents[i].phi_alpha * rho_hat_h;

        if (m_constituents[i].eta_alpha_h >= 0) { //for anisotropic m_constituents[i]s
            m_constituents[i].G_alpha_h = m_constituents[i].g_alpha_h * mat3dd(0.0, sin(m_constituents[i].eta_alpha_h), cos(m_constituents[i].eta_alpha_h));
        }
        else { //for isotropic m_constituents[i]s (i.e. elastin)
            m_constituents[i].G_alpha_h = mat3dd(m_constituents[i].g_alpha_r, m_constituents[i].g_alpha_theta, m_constituents[i].g_alpha_z);
        }

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

    ar & nts;
    ar & sn;
    ar & m_dt;
    ar & m_nconstituents;
    ar & K_delta_tauw;
    ar & K_delta_sigma;
    ar & sigma_inv_h;
    ar & sigma_inv_curr;
    ar & rho_hat_h;
    ar & m_constituents;
    ar & m_lambda_act;
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
    //fflush(stdout);

    sn = int(sn + dt);

    sigma_inv_curr = et.m_s.tr();
    update_kinetics(sn);
    update_sigma(sn);
    et.m_J_star = m_J_s[sn];

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

    // push deformation gradient to local coordinates
    mat3d Q = mat3d(e_r(mp), e_t(mp), e_z(mp));
    pt.K_delta_sigma = 0.184*m_a_val(mp);
    pt.m_constituents[0].c1_alpha_h = 89.710*(1.0 - 0.595*m_a_val(mp));

    pt.m_F_s[sn] = Q.transpose() * et.m_F * Q;

    // calculate the stress as a sum of deviatoric stress and pressure
    pt.update_sigma(sn);

    mat3ds sbar = pow(J / et.m_J_star, -2.0 / 3.0) * (Q * pt.m_sigma * Q.transpose()).sym();

    devstress = sbar.dev();

    const mat3dd  I(1.0);
    const tens4ds IxI = dyad1s(I);
    const tens4ds IoI = dyad4s(I);

    tens4ds cbar = pow(J / et.m_J_star, -4.0 / 3.0) * pt.m_CC.pp(Q);

    devtangent = cbar - 1./3.*(ddots(cbar, IxI) - IxI*(cbar.tr()/3.))
    + 2./3.*((IoI-IxI/3.)*sbar.tr()-dyad1s(sbar.dev(),I));

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
    //tangent = dev_tangent + (IxI - IoI*2)*p_val + IxI*(UJJ(J, pt.m_J_curr)*J);

    // Write outputs to unused velocity and acceleration vectors for writeout
    // TODO: confirm what "stiffness" means in this context

    et.m_v.x = pt.m_J_s[sn];                             // Target volume
    et.m_v.y = m_a_val(mp);                             // Aneurysm injury
    et.m_v.z = pt.m_constituents[0].rhoR_alpha[sn];     // Generally, elasticity density
    et.m_a.x = pt.m_CC(0,0,0,0);     // Radial stiffness
    et.m_a.y = pt.m_CC(1,1,1,1);     // Circumfrential stiffness
    et.m_a.z = pt.m_CC(2,2,2,2);     // Axial stiffness


}

void GRMaterialPoint::update_kinetics(int sn) {


    //This function updates the kinetics for G&R.
    int taun_min = 0;

    //Differences in current mechanical state from the reference state
    //Stress invariant
    double delta_sigma = ((1 - K_delta_sigma) * (sigma_inv_curr / sigma_inv_h)) - 1;
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

    int n = 0; //number of points in integration interval

    n = (sn - taun_min) + 1; //find number of integration pts

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
                upsilon_p = 1 + K_sigma_p * delta_sigma - K_tauw_p * delta_tauw;
                if (upsilon_p < 0.1) {
                    upsilon_p = 0.1;
                }

                upsilon_d = 1 + K_sigma_d * pow(delta_sigma, 2) + K_tauw_d * pow(delta_tauw, 2);

                if (m_constituents[alpha].rhoR_alpha[sn] <= m_constituents[alpha].rhoR_alpha_h) {
                    rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha_h;
                }
                else {
                    rhoR_alpha_calc = m_constituents[alpha].rhoR_alpha[sn];
                }
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


    //Get current time index
    int taun_min = 0;

    mat3ds  full1(1.0,1.0,1.0,1.0,1.0,1.0);
    tens4ds full11 = dyad1s(full1);

    //Find the current deformation gradient
    mat3d F_s = m_F_s[sn];
    double J_s = F_s.det(); //m_J_s[sn];

    // compute U from polar decomposition of deformation gradient tensor
    mat3ds U; mat3d R; F_s.right_polar(R,U);
    mat3ds U_tau; mat3d R_tau;

    // TODO: Consider initializing from passed value? bc of special
    mat3ds C_s = (F_s.transpose()*F_s).sym();

    //Calculate constituent specific stretches for evolving constituents at the current time
    for  (int alpha=0; alpha<m_nconstituents; ++alpha) {

        //Check to see if constituent is isotropic
        if (m_constituents[alpha].eta_alpha_h >= 0) {

            //Stretch is equal to the sqrt of I4
            m_constituents[alpha].lambda_alpha_tau[sn] = sqrt(C_s(2, 2) * pow(cos(m_constituents[alpha].eta_alpha_h), 2)
                + C_s(1, 1) * pow(sin(m_constituents[alpha].eta_alpha_h), 2));

        }

    }

    //Find the mechanical contributions of each constituent for each direction
    //double a, h;
    double J_tau = 1;
    mat3d F_tau  = m_F_s[sn];
    mat3d G_alpha_h_N(0.0);
    mat3d G_alpha_h(0.0);
    double lambda_alpha_ntau_s = 0;
    mat3d F_alpha_ntau_s(1, 0, 0, 0, 1, 0, 0, 0, 1);
    double hat_S_alpha = 0;
    mat3ds sigma(0);
    double g_alpha = 1.0;

    //Local active variables
    double C = 0;
    double lambda_act = 0;
    double lambda_act_1 = 0;
    double lambda_act_2 = 0;
    double parab_act = 0;
    double hat_sigma_act = 0, sigma_act = 0, S_act = 0;
    mat3ds sigma_act_mat(0);
    double Cbar_act = 0;

    m_lambda_act[sn] = sqrt(C_s(1, 1));

    //Stiffness variables
    double hat_dS_dlambda2_alpha = 0;
    double dSdC_act = 0;
    vector<double> constitutive_return = { 0, 0 };

    //Full material stiffness tensor
    tens4ds CC(0.0);
    tens4ds CC_act_mat(0.0);
    tens4ds hat_CC_1(0.0);
    tens4ds hat_CC_2(0.0);

    //Integration variables
    //For mass
    double mq_1 = 0, mq_2 = 0;
    double q_1 = 1.0, q_2 = 1.0;
    double k_1 = 0, k_2 = 0;

    int n = 0; //number of pts in integration interval

    //For stress
    mat3ds hat_sigma_1(0);
    mat3ds hat_sigma_2(0);
    //For active stress
    double q_act_1 = 0, q_act_2 = 0;

    n = (sn - taun_min) + 1; //number of integration pts

    int alpha_count = 0;
        

    for  (int alpha=0; alpha<m_nconstituents; ++alpha) {

        alpha_count = alpha_count + 1;

        //Trapz rule allows for fast heredity integral evaluation
        k_2 = m_constituents[alpha].k_alpha[sn];
        q_2 = 1.0;
        mq_2 = m_constituents[alpha].mR_alpha[sn];

        //Find active radius from current cohort
        if (m_constituents[alpha].m_active == 1) {
            lambda_act_2 = sqrt(C_s(1, 1))/m_lambda_act[sn];
            q_act_2 = 1.0;
        }

        //Kinematics
        F_tau = F_s;
        J_tau = J_s;

        G_alpha_h = m_constituents[alpha].G_alpha_h;

        // compute U from polar decomposition of deformation gradient tensor
        F_tau.right_polar(R_tau,U_tau);
        G_alpha_h_N = G_alpha_h * R_tau;

        constitutive_return = m_constituents[alpha].constitutive(sn, sn);

        hat_S_alpha = constitutive_return[0];
        hat_dS_dlambda2_alpha = constitutive_return[1];
        F_alpha_ntau_s = F_s*F_tau.inverse()*G_alpha_h_N;

        hat_sigma_2 = 1.0/J_s * (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s.transpose()).sym();
        hat_CC_2 = hat_dS_dlambda2_alpha * 2.0 /J_s * full11.pp(F_alpha_ntau_s);


        //Check if during G&R or at initial time point
        if (sn > 0 && m_constituents[alpha].m_degradable) {

            for (int taun = sn - 1; taun >= taun_min; taun = taun - 1) {

                F_tau = m_F_s[taun];
                J_tau = F_tau.det(); //m_J_s[taun];
                F_tau.right_polar(R_tau,U_tau);
                G_alpha_h_N = G_alpha_h*R_tau;

                //Find 1st intermediate kinetics
                k_1 = m_constituents[alpha].k_alpha[taun];
                q_1 = exp(-(k_2 + k_1) * m_dt / 2) * q_2;
                mq_1 = m_constituents[alpha].mR_alpha[taun] * q_1;

                //Find intermediate active state
                if (m_constituents[alpha].m_active == 1) {
                    lambda_act_1 = sqrt(C_s(1, 1))/m_lambda_act[taun];
                    q_act_1 = exp(-k_act * m_dt) * q_act_2;
                }

                constitutive_return = m_constituents[alpha].constitutive(sn, taun);
                hat_S_alpha = constitutive_return[0];
                hat_dS_dlambda2_alpha = constitutive_return[1];
                F_alpha_ntau_s = F_s*F_tau.inverse()*G_alpha_h_N;

                hat_sigma_1 = 1.0 / J_s * (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s.transpose()).sym();
                hat_CC_1 = hat_dS_dlambda2_alpha * 2.0 /J_s * full11.pp(F_alpha_ntau_s);

                // Add integration contribution
                sigma += (mq_2 * hat_sigma_2 + mq_1 * hat_sigma_1) / m_constituents[alpha].rho_hat_alpha_h * m_dt / 2;
                CC += (mq_2 * hat_CC_2 + mq_1 * hat_CC_1) / m_constituents[alpha].rho_hat_alpha_h * m_dt / 2;

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
                hat_sigma_2 = hat_sigma_1;
                hat_CC_2 = hat_CC_1;


            }

            // Add in the stress and stiffness contributions of the initial material
            if (taun_min == 0) {
                sigma += m_constituents[alpha].rhoR_alpha[0]
                    / m_constituents[alpha].rho_hat_alpha_h * q_1 * hat_sigma_1;

                CC += m_constituents[alpha].rhoR_alpha[0]
                    / m_constituents[alpha].rho_hat_alpha_h * q_1 * hat_CC_1;
            }

        }
        //Initial time point and constituents with prescribed degradation profiles
        else {
            //Find stress from initial cohort
            constitutive_return = m_constituents[alpha].constitutive(sn, 0);
            hat_S_alpha = constitutive_return[0];
            hat_dS_dlambda2_alpha = constitutive_return[1];

            F_alpha_ntau_s = F_s*G_alpha_h;

            //Check if anisotropic
            hat_sigma_2 = 1.0/J_s * (F_alpha_ntau_s * hat_S_alpha * F_alpha_ntau_s.transpose()).sym();

            hat_CC_2 = hat_dS_dlambda2_alpha * 2.0 /J_s * full11.pp(F_alpha_ntau_s);

            // Add in stress and stiffness contributions
            sigma += m_constituents[alpha].rhoR_alpha[sn] /
                m_constituents[alpha].rho_hat_alpha_h * hat_sigma_2;                

            CC += m_constituents[alpha].rhoR_alpha[sn] /
                m_constituents[alpha].rho_hat_alpha_h * hat_CC_2;
        }


        if (taun_min == 0 && m_constituents[alpha].m_active == 1) {
            if (!m_constituents[alpha].m_degradable){
                q_act_1 = exp(-k_act * m_dt) * q_act_2;
            }
            lambda_act += sqrt(C_s(2,2))/m_lambda_act[0] * q_act_1;
        }


        if (m_constituents[alpha].m_active == 1) {
        
            //Find active stress contribtion
            //add in initial active stress radius contribution
            C = CB - CS * (bar_tauw_curr / bar_tauw_h - 1);

            if (sn == 0) {
                C = CB;
                lambda_act = 1.0;
            }

            parab_act = 1 - pow((lambda_m - lambda_act) / (lambda_m - lambda_0), 2);

            S_act = T_act * (1 - exp(-pow(C, 2))) * pow(lambda_act, -1) * parab_act *
                    m_constituents[alpha].rhoR_alpha[sn] / m_constituents[alpha].rho_hat_alpha_h;


            dSdC_act = T_act *  (1 - exp(-pow(C, 2))) * (1/lambda_act) *
                       ((-1/pow(lambda_act,2))*parab_act + 2*(1/lambda_act)*(lambda_m - lambda_act)/ pow(lambda_m - lambda_0, 2)) *
                       m_constituents[alpha].rhoR_alpha[sn] / m_constituents[alpha].rho_hat_alpha_h;


            sigma_act_mat = 1.0/J_s * (F_s * mat3dd(0,1,0) * F_s.transpose()).sym();
            //TODO : Check this

            CC_act_mat = tens4ds(0.0);
            CC_act_mat(1,1,1,1) = hat_dS_dlambda2_alpha;
            CC_act_mat = 1.0 / J_s * CC_act_mat.pp(F_s);

            sigma += sigma_act_mat;
            CC += CC_act_mat;

        }

    }

    m_sigma = sigma;
    m_CC = CC;

    return;
}

