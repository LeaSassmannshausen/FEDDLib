#ifndef POWERLAW_DEF_hpp
#define POWERLAW_DEF_hpp

#include "PowerLaw_decl.hpp"

namespace FEDD {

template <class SC, class LO, class GO, class NO>
PowerLaw<SC,LO,GO,NO>::PowerLaw(ParameterListPtr_Type params):
DifferentiableFuncClass<SC,LO,GO,NO>(params)
{
    this->params_=params;
    // Reading through parameterlist
    shearThinningModel_ = this->params_->sublist("Material").get("ShearThinningModel","");
	powerlaw_constant_K = this->params_->sublist("Material").get("PowerLawParameter K",0.);            // corresponds to K in the formulas in the literature
    powerlaw_index_n    = this->params_->sublist("Material").get("PowerLaw index n",1.);                         // corresponds to n in the formulas being the power-law index     
    nu_zero             = this->params_->sublist("Material").get("Numerical_ZeroShearRateViscosity",1.);  
    nu_infty            =  this->params_->sublist("Material").get("Numerical_InftyShearRateViscosity",1.);  
    viscosity_ = 0.;
    //	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No discretisation Information for Velocity in Navier Stokes Element." );

	
}


template <class SC, class LO, class GO, class NO>
void PowerLaw<SC,LO,GO,NO>::evaluateFunction(ParameterListPtr_Type params, double shearRate, double &viscosity) {
	
    
    viscosity = this->powerlaw_constant_K*pow(shearRate, this->powerlaw_index_n-1.0);

    // Check bounds for viscosity
    if (viscosity < nu_infty)
    {
    viscosity = nu_infty;
    }
    else if (abs(viscosity) > nu_zero)
    {
    viscosity = nu_zero;
    }

    this-> viscosity_ = viscosity;
}


template <class SC, class LO, class GO, class NO>
void PowerLaw<SC,LO,GO,NO>::evaluateDerivative(ParameterListPtr_Type params, double shearRate, double &res) {
	
// The function is composed of d_eta/ d_GammaDot * d_GammaDot/ D_Tau while d_GammaDot * d_GammaDot/ D_Tau= - 2/GammaDot
// So a problematic case is if shear rate is in the denominator 0 because we may get inf values. 
// Therefore we have to check these cases and catch them
double epsilon = 10e-6;
if ( abs(shearRate) <= epsilon) //How to choose epsilon?
       {
            shearRate =  epsilon;
       }

res = (-2.0)*this->powerlaw_constant_K*(this->powerlaw_index_n-1.0)*pow(shearRate, this->powerlaw_index_n-3.0);


}


template <class SC, class LO, class GO, class NO>
void PowerLaw<SC,LO,GO,NO>::setParams(ParameterListPtr_Type params){
    this->params_=params;
    // Reading through parameterlist
    shearThinningModel_ = this->params_->sublist("Material").get("ShearThinningModel","");
	powerlaw_constant_K = this->params_->sublist("Material").get("PowerLawParameter K",0.);            // corresponds to K in the formulas in the literature
    powerlaw_index_n    = this->params_->sublist("Material").get("PowerLaw index n",1.);                         // corresponds to n in the formulas being the power-law index     
    nu_zero             = this->params_->sublist("Material").get("Numerical_ZeroShearRateViscosity",1.);  
    nu_infty            =  this->params_->sublist("Material").get("Numerical_InftyShearRateViscosity",1.);  
 }




template <class SC, class LO, class GO, class NO>
void PowerLaw<SC,LO,GO,NO>::echoParams(){
            std::cout << "************************************************************ "  <<std::endl;
            std::cout << "-- Chosen material model ..." << this->shearThinningModel_ << " --- "  <<std::endl;
            std::cout << "-- Specified material parameters:" <<std::endl;
            std::cout << "-- PowerLaw index n:" << this->powerlaw_index_n << std::endl;
            std::cout << "-- PowerLaw constant K:" << this->powerlaw_constant_K <<std::endl;
            std::cout << "************************************************************ "  <<std::endl;
  }




}
#endif

