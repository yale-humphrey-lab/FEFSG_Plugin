//#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material. 
// It requires FEBio 2.5 (or up)
//
// Author Steve Maas
// Copyright (c) 2015 - 2016
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial which is defined in this include files.
#include "FEBioMech/FEElasticMaterial.h"
#include <iostream>								// to use cin.get()

class GRConstituent {
public:
    double epsilon_pol_min;
    double eta_alpha_h;
    double c1_alpha_h;
    double c2_alpha_h;
    double g_alpha_h;
    mat3dd G_alpha_h;
    double k_alpha_h;
    double rho_hat_alpha_h;
    double rhoR_alpha_h;
    std::vector<double> mR_alpha;
    std::vector<double> k_alpha;
    std::vector<double> rhoR_alpha;
    std::vector<double> epsilonR_alpha;
    std::vector<double> lambda_alpha_tau;
    std::vector<double> epsilon_alpha;
    double K_sigma_p_alpha_h;
    double K_tauw_p_alpha_h;
    double K_sigma_d_alpha_h;
    double K_tauw_d_alpha_h;
    bool m_degradable;
    bool m_inflammatory;
    bool m_active;
    bool m_polymer;

    double g_alpha_r;
    double g_alpha_theta;
    double g_alpha_z;
    double phi_alpha;


    // Default constructor
    GRConstituent() {
        Init();
    }

    // Initialize function
    void Init() {
        epsilon_pol_min = 0.0;
        eta_alpha_h = 0.0;
        c1_alpha_h = 0.0;
        c2_alpha_h = 0.0;
        g_alpha_h = 0.0;
        G_alpha_h = mat3dd(0.0);
        k_alpha_h = 0.0;
        rho_hat_alpha_h = 0.0;
        rhoR_alpha_h = 0.0;
        mR_alpha = std::vector<double>(720, 0.0);  // Vector of zeros of length 720
        k_alpha = std::vector<double>(720, 0.0);   // Vector of zeros of length 720
        rhoR_alpha = std::vector<double>(720, 0.0); // Vector of zeros of length 720
        epsilonR_alpha = std::vector<double>(720, 0.0);  // Vector of zeros of length 720
        epsilon_alpha = std::vector<double>(720, 0.0);  // Vector of zeros of length 720
        lambda_alpha_tau = std::vector<double>(720, 1.0);  // Vector of zeros of length 720
        K_sigma_p_alpha_h = 0.0;
        K_tauw_p_alpha_h = 0.0;
        K_sigma_d_alpha_h = 0.0;
        K_tauw_d_alpha_h = 0.0;
        m_degradable = false;
        m_inflammatory = false;
        m_active = false;
        m_polymer = false;
    }


    // Serialization function for GRConstituent
    void Serialize(DumpStream& ar) {
	    ar & epsilon_pol_min;
	    ar & eta_alpha_h;
	    ar & c1_alpha_h;
	    ar & c2_alpha_h;
	    ar & g_alpha_h;
	    ar & G_alpha_h;
	    ar & k_alpha_h;
	    ar & rho_hat_alpha_h;
	    ar & rhoR_alpha_h;
	    ar & mR_alpha;
	    ar & k_alpha;
	    ar & rhoR_alpha;
	    ar & epsilonR_alpha;
	    ar & epsilon_alpha;
	    ar & lambda_alpha_tau;
	    ar & K_sigma_p_alpha_h;
	    ar & K_tauw_p_alpha_h;
	    ar & K_sigma_d_alpha_h;
	    ar & K_tauw_d_alpha_h;
	    ar & m_degradable;
	    ar & m_inflammatory;
	    ar & m_active;
	    ar & m_polymer;
    }


	vector<double> constitutive(int sn, int ts) {

	    double lambda_alpha_ntau_s = 0;
	    double Q1 = 0;
	    double Q2 = 0;
	    double hat_S_alpha = 0;
	    double hat_dS_dlambda2_alpha = 0;
	    double pol_mod = 0;
	    double epsilon_curr = 0;
	    vector<double> return_constitutive = { 0, 0 };

	    //Check if ansisotropic
	    if ( eta_alpha_h >= 0) {

	        lambda_alpha_ntau_s =  g_alpha_h * lambda_alpha_tau[sn] /  lambda_alpha_tau[ts];

	        if (lambda_alpha_ntau_s < 1) {
	            lambda_alpha_ntau_s = 1;
	        }

	        Q1 = (pow(lambda_alpha_ntau_s, 2) - 1);
	        Q2 =  c2_alpha_h * pow(Q1, 2);
	        hat_S_alpha =  c1_alpha_h * Q1 * exp(Q2);
	        hat_dS_dlambda2_alpha =  c1_alpha_h * exp(Q2) * (1 + 2 * Q2);

	        /*
			printf("sn: %d\n", sn);
			printf("ts: %d\n", ts);
			printf("g_alpha_h: %f\n", g_alpha_h);
			printf("lambda_alpha_tau[sn]: %f\n", lambda_alpha_tau[sn]);
			printf("lambda_alpha_tau[ts]: %f\n", lambda_alpha_tau[ts]);
			printf("lambda_alpha_ntau_s: %f\n", lambda_alpha_ntau_s);
			printf("Q1: %f\n", Q1);
			printf("Q2: %f\n", Q2);
			printf("hat_S_alpha: %f\n", hat_S_alpha);
			printf("hat_dS_dlambda2_alpha: %f\n", hat_dS_dlambda2_alpha);
			fflush(stdout);
			*/

	    }
	    else {

	        if ( m_polymer) {

	            if ( epsilon_alpha[sn] <  epsilon_pol_min) {
	                 epsilon_pol_min =  epsilon_alpha[sn];
	            }
	            //TODO Carefull with assignment (pointer?)
	            epsilon_curr =  epsilon_pol_min;
	            pol_mod = 0.03 * pow(epsilon_curr, 2);
	        }

	        else {
	            pol_mod = 1;
	        }

	        hat_S_alpha = pol_mod *  c1_alpha_h;
	    }

	    return_constitutive = { hat_S_alpha , hat_dS_dlambda2_alpha };
	    return return_constitutive;

	}
};

class FEBIOMECH_API GRMaterialPoint : public FEMaterialPointData
{
public:
	GRMaterialPoint(FEMaterialPointData *pt) : FEMaterialPointData(pt) {};

	FEMaterialPointData* Copy() override;

	void Init() override;
	void Serialize(DumpStream& ar) override;

	void update_sigma(double dt, int sn);
	void update_kinetics(double dt, int sn);


public:
    int nts;
    int sn;
    double K_delta_tauw;
    double K_delta_sigma;
    double sigma_inv_h;
    double sigma_inv_curr;
    double rho_hat_h;
    std::vector<GRConstituent> m_constituents;
    std::vector<double> m_lambda_act;
    std::vector<mat3d> m_F_s;  // Vector of mat3d variables
    std::vector<double> m_J_s;
    std::vector<double> rhoR;
    std::vector<double> rho;
    std::vector<double> ups_infl_p;
    std::vector<double> ups_infl_d;
    double CB;
    double CS;
    double bar_tauw_curr;
    double bar_tauw_h;
    double lambda_m;
    double lambda_0;
    double T_act;
    double k_act;
    mat3ds m_sigma;
    tens4dmm m_CC;
};


//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model. 
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial. 
class FEFSG : public FEElasticMaterial
{
public:
	FEFSG(FEModel* pfem);

	bool UseSecantTangent() override { return m_secant_tangent; }

	//! create material point data for this material
	FEMaterialPointData* CreateMaterialPointData() override;


public:
	// The constructor is called when an instance of this class is created.
	// All classes registered by the framework must take the FEModel* as the only
	// parameter in the constructor, even if the class does not need it (which most often
	// will be the case). For material classes, the FEModel parameter is passed to the 
	// base class in the initialization list.

	// setting m_secant_tangent = true so FESolidMaterial uses SecantTangent
	// (allows minor symmetry only tangents) instead of Tangent (minor and major symmetries)
	// 	 { return m_secant_tangent; }
    bool m_secant_tangent;   //!< flag for using secant tangent

public:
	// function to perform material evaluation. calculates stress and tangent to avoid code duplication
	void StressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4dmm& tangent);

	// This function calculates the spatial (i.e. Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint and returns a mat3ds object
	// which is a symmetric second-order tensor.
	virtual mat3ds Stress(FEMaterialPoint& pt) override {
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return stress;
	}

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
	// which is a fourth-order tensor with major and minor symmetries.
	virtual tens4ds Tangent(FEMaterialPoint& pt) override {
		tens4ds tangent;
		printf("!!!!!!!!!!!!!!!!!!!!Using TANGENT \n");
		fflush(stdout);
		return tangent;
	};

	// minor symmetries only
	virtual tens4dmm SecantTangent(FEMaterialPoint& pt, bool mat) override {
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return tangent;
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

public:
	double	m_K;			//!< bulk modulus
	int     m_npmodel;      //!< pressure model for U(J)
	double     m_dt;      //!< pressure model for U(J)

    DECLARE_FECORE_CLASS();


};
