#pragma once
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
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FEModelParam.h>
#include <iostream>								// to use cin.get()

class GRConstituent {
public:
	enum { MAX_TIMESTEPS = 2881 };

    double epsilon_pol_min;
    double eta_alpha_h;
    double c1_alpha_h;
    double c2_alpha_h;
    double c1_alpha;
    double c2_alpha;
    double g_alpha_h;
    double k_alpha_h;
    double rho_hat_alpha_h;
    double rhoR_alpha_h;
    double mR_alpha[MAX_TIMESTEPS];
    double k_alpha[MAX_TIMESTEPS];
    double rhoR_alpha[MAX_TIMESTEPS];
    double epsilonR_alpha[MAX_TIMESTEPS];
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
        c1_alpha = 0.0;
        c2_alpha = 0.0;
        g_alpha_h = 0.0;
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
	    ar & c1_alpha;
	    ar & c2_alpha;
	    ar & g_alpha_h;
	    ar & k_alpha_h;
	    ar & rho_hat_alpha_h;
	    ar & rhoR_alpha_h;
	    ar & mR_alpha;
	    ar & k_alpha;
	    ar & rhoR_alpha;
	    ar & epsilonR_alpha;
	    ar & epsilon_alpha;
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


	void constitutive(mat3d F_s, mat3d F_tau, int sn, mat3ds& stress, tens4ds& tangent) {

	    double pol_mod = 0;
	    double epsilon_curr = 0;
    
	    // Evaluate the elasticity tensor
	    mat3dd I(1);
	    tens4ds IxI = dyad1s(I);
	    tens4ds I4  = dyad4s(I);

	    stress = mat3ds(0.0);
	    tangent = tens4ds(0.0);

	    //Check if ansisotropic
	    if ( eta_alpha_h >= 0) {


	    	mat3dd G = mat3dd(g_alpha_h);
		    mat3d F = F_s*F_tau.inverse()*G;
		    mat3ds C = (F.transpose()*F).sym();
    		mat3ds U; mat3d R; F_tau.right_polar(R,U);

		    // Copy the local element basis directions to n
			vec3d n[2];
		    n[0].x = 0; n[0].y = 1; n[0].z = 0;
		    n[1].x = 0; n[1].y = 0; n[1].z = 1;
		    
		    // Evaluate the structural direction in the current configuration
		    vec3d ar,a;

			double lth = (F_tau * n[0]).norm();	// lth -> 1 for Fh -> Fo
			double lzh = (F_tau * n[1]).norm();	// lzh -> 1 for Fh -> Fo

			double eta_alpha_curr = eta_alpha_h;
			//double aexp = 1.0;
			//if ((eta_alpha_h != 0.0) and (eta_alpha_curr != PI/2)){
			//	eta_alpha_curr = atan(tan(eta_alpha_h)*pow(lth/lzh,aexp));	// Remodeled angle
			//}

		    double cg = cos(eta_alpha_curr); double sg = sin(eta_alpha_curr);
		    ar = n[0]*sg + n[1]*cg;
		    ar = R.transpose()*ar;
		    //ar = F_tau*ar;
		    //ar = ar/ar.norm();

		    a = F*ar;
		    // Evaluate the structural tensors in the current configuration
		    // and the fiber strains and stress contributions
		    double I40 = (ar*(C*ar));
		    double E0 = I40-1;

		    mat3ds h0;
		    //if (E0 >= 0) {
		        h0 = dyad(a);
		        stress += h0*(c1_alpha*E0*exp(c2_alpha*E0*E0));
		        tangent += dyad1s(h0)*(2.*c1_alpha*(1. + 2. * c2_alpha*E0*E0)*exp(c2_alpha*E0*E0));
    		//}


	    }
	    else {

	        if (m_polymer) {

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

	    	mat3dd G = mat3dd(g_alpha_r, g_alpha_theta, g_alpha_z);
		    mat3d  F = F_s*F_tau.inverse()*G;
		    mat3ds b = (F*F.transpose()).sym();

	        stress += pol_mod * c1_alpha * b;
	    }
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
	
	enum { MAX_TIMESTEPS = 2881 };
	enum { MAX_CONSTITUENTS = 6 };


public:
    int sn;
    int m_nconstituents;
    double m_dt;
    double K_delta_tauw;
    double K_delta_sigma;
    double sigma_inv_h;
    double sigma_inv_curr;
    double rho_hat_h;
    GRConstituent m_constituents[MAX_CONSTITUENTS];
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
    mat3d m_F_curr;
    double m_J_curr;
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
	FEParamDouble     m_e_injury_val;      //!< K_delta_sigma
	FEParamDouble     m_k_injury_val;      //!< K_delta_sigma

    DECLARE_FECORE_CLASS();


};