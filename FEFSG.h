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
#include "FEBioMech/FEUncoupledMaterial.h"
#include <FECore/FEModelParam.h>
#include <iostream>								// to use cin.get()

class GRConstituent {
public:
	enum { MAX_TIMESTEPS = 720 };

    double epsilon_pol_min;
    double eta_alpha_h;
    double c1_alpha_h;
    double c2_alpha_h;
    double g_alpha_h;
    mat3dd G_alpha_h;
    double k_alpha_h;
    double rho_hat_alpha_h;
    double rhoR_alpha_h;
    double mR_alpha[MAX_TIMESTEPS];
    double k_alpha[MAX_TIMESTEPS];
    double rhoR_alpha[MAX_TIMESTEPS];
    double epsilonR_alpha[MAX_TIMESTEPS];
    double lambda_alpha_tau[MAX_TIMESTEPS];
    double epsilon_alpha[MAX_TIMESTEPS];
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


		for (int i=0; i<MAX_TIMESTEPS; ++i)
		{
	        mR_alpha[i] = 0;
	        k_alpha[i] = 0;
	        rhoR_alpha[i] = 0;
	        epsilonR_alpha[i] = 0;
	        epsilon_alpha[i] = 0;
	        lambda_alpha_tau[i] = 1;
		}

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
	    ar & g_alpha_r;
	    ar & g_alpha_theta;
	    ar & g_alpha_z;
	    ar & phi_alpha;
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
		        hat_S_alpha =  0;
		        hat_dS_dlambda2_alpha =  0;
	        } else {
		        Q1 = (pow(lambda_alpha_ntau_s, 2) - 1);
		        Q2 =  c2_alpha_h * pow(Q1, 2);
		        hat_S_alpha =  c1_alpha_h * Q1 * exp(Q2);
		        hat_dS_dlambda2_alpha =  c1_alpha_h * exp(Q2) * (1 + 2 * Q2);
	        }

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
	//! Update material point data
	void Update(const FETimeInfo& timeInfo);

	void update_sigma(int sn);
	void update_kinetics(int sn);

	enum { MAX_TIMESTEPS = 720 };
	enum { MAX_CONSTITUENTS = 6 };


public:
    int nts;
    int sn;
    double m_nconstituents;
    double m_dt;
    double K_delta_tauw;
    double K_delta_sigma;
    double sigma_inv_h;
    double sigma_inv_curr;
    double rho_hat_h;
    GRConstituent m_constituents[MAX_CONSTITUENTS];
    double m_lambda_act[MAX_TIMESTEPS];
    mat3d  m_F_s[MAX_TIMESTEPS];
    double m_J_s[MAX_TIMESTEPS];
    double rhoR[MAX_TIMESTEPS];
    double rho[MAX_TIMESTEPS];
    double ups_infl_p[MAX_TIMESTEPS];
    double ups_infl_d[MAX_TIMESTEPS];
    double CB;
    double CS;
    double bar_tauw_curr;
    double bar_tauw_h;
    double lambda_m;
    double lambda_0;
    double T_act;
    double k_act;
    mat3ds m_sigma;
    tens4ds m_CC;
};


//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model. 
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial. 
class FEFSG : public FEUncoupledMaterial
{
public:
	FEFSG(FEModel* pfem);

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
    //bool m_secant_tangent;   //!< flag for using secant tangent

public:
	// function to perform material evaluation. calculates stress and tangent to avoid code duplication
	void DevStressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4ds& tangent);

	// This function calculates the spatial (i.e. Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint and returns a mat3ds object
	// which is a symmetric second-order tensor.
	virtual mat3ds DevStress(FEMaterialPoint& pt) override {
		mat3ds devstress;
		tens4ds devtangent;
		DevStressTangent(pt, devstress, devtangent);
		return devstress;
	}

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
	// which is a fourth-order tensor with major and minor symmetries.
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override {
		mat3ds devstress;
		tens4ds devtangent;
		DevStressTangent(pt, devstress, devtangent);
		return devtangent;
	};

public:
	FEParamVec3     e_r;      //
	FEParamVec3     e_t;      //
	FEParamVec3     e_z;      //
	FEParamDouble     m_a_val;      //!< K_delta_sigma

    DECLARE_FECORE_CLASS();


};