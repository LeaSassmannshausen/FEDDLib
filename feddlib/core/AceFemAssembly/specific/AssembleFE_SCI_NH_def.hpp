#ifndef ASSEMBLEFE_SCI_NH_DEF_hpp
#define ASSEMBLEFE_SCI_NH_DEF_hpp

#include "AssembleFE_SCI_NH_decl.hpp"
#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.hpp"
#endif

#include <vector>
// #include <iostream>

namespace FEDD
{

	template <class SC, class LO, class GO, class NO>
	AssembleFE_SCI_NH<SC, LO, GO, NO>::AssembleFE_SCI_NH(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple)
	{
		// Extracting values from ParameterList
		int numMaterials =  this->params_->sublist("Parameter Solid").get("Number of Materials", 1);
		int materialID = 0;

		for(int i=1; i<= numMaterials; i++)
			if( this->params_->sublist("Parameter Solid").sublist(std::to_string(i)).get("Volume Flag", 15) == this->flag_)
				materialID = i;
		
		if(materialID == 0)
			cout << "!!! Warning: No corresponding parameterslist for the element flag="<< this->flag_ << ". Please Check volume flags of elements and Mesh Data !!! " << endl;


		E0_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("E", 0.38);
		E1_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("E1", 0.3);
		poissonRatio_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Poisson Ratio", 0.49e-0);
		c1_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("c1", 0.25e-0);
		D0_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("D0",6.e-05);
		m_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("m", 0.0);
		dofOrdering_ = this->params_->sublist("Parameter").sublist(std::to_string(materialID)).get("Ordering", 2);
		//cout << "--- Init AssembleFE_SCI_NH Element --- EMOD " << E0_  << " D0_ " << D0_ << endl;

		FEType_ = std::get<1>(this->diskTuple_->at(0));	   // FEType of Disk
		dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
		dofsChem_ = std::get<2>(this->diskTuple_->at(1));  // Degrees of freedom per node

		numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
		numNodesChem_ = std::get<3>(this->diskTuple_->at(1));  // Number of nodes of element

		dofsElement_ = dofsSolid_ * numNodesSolid_ + dofsChem_ * numNodesChem_; // "Dimension of return matrix"

		solution_n_.resize(60, 0.);
		solution_n1_.resize(60, 0.);

		this->postProcessingData_ = Teuchos::rcp( new SmallMatrix_Type(dofsElement_,0.));

		//cout << " Parameters:: E=" << E0_ << " E1=" << E1_ << " poissionRation=" << poissonRatio_ << endl; 
		/*timeParametersVec_.resize(0, vec_dbl_Type(2));
		numSegments_ = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

		for(int i=1; i <= numSegments_; i++){

			double startTime = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
			double dtTmp = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
			
			vec_dbl_Type segment = {startTime,dtTmp};
			timeParametersVec_.push_back(segment);
		}*/

	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_NH<SC, LO, GO, NO>::assembleJacobian()
	{

		SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(dofsElement_, 0.));

		assembleDeformationDiffusionNeoHook(elementMatrix);
		// elementMatrix->print();
		this->jacobian_ = elementMatrix;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_NH<SC, LO, GO, NO>::advanceInTime(double dt)
	{

		// If we have a time segment setting we switch to the demanded time increment
		/*for(int i=0; i<numSegments_ ; i++){
			if(this->timeStep_ +1.0e-12 > timeParametersVec_[i][0])
				this->timeIncrement_=timeParametersVec_[i][1];
		}*/
		if(this->timeStep_ -1.e-13 < 0) // only in this one instance T=0 we set the dt beforehand, as the initial dt is set through the paramterlist and this is error prone
			this->timeIncrement_=dt;

		this->timeStep_ = this->timeStep_ + this->timeIncrement_;

		this->timeIncrement_ = dt;

		if(this->globalElementID_==0){
			cout << " ---------------------------------------------- " << endl;
			cout << " AssembleFE_SCI_NH: Advancing time in elements" << endl;
			cout << " Timestep: " << this->timeStep_ << " \t timeincrement: "<< this->timeIncrement_ << endl;
			cout << " ---------------------------------------------- " << endl;

		}
        //cout << " Advance in time on element timestep: " << this->timeStep_ << " increment: " << this->timeIncrement_ << endl;

		for (int i = 0; i < 40; i++)
		{
			if (i < 30)
				solution_n_[i] = (*this->solution_)[i];
			else
				solution_n_[i] = (*this->solution_)[i];
		}
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_NH<SC, LO, GO, NO>::assembleRHS()
	{

		this->rhsVec_.reset(new vec_dbl_Type(dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

		std::vector<double> positions(30);
		int count = 0;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}

		std::vector<double> displacements(30);
		for (int i = 0; i < 30; i++){
			displacements[i] = (*this->solution_)[i];
		}
		std::vector<double> concentrations(10);
		for (int i = 0; i < 10; i++)
			concentrations[i] = (*this->solution_)[i + 30];

		std::vector<double> concentrationsLastConverged(10);
		for (int i = 0; i < 10; i++)
			concentrationsLastConverged[i] = (solution_n_)[i + 30];

		std::vector<double> domainData(6);
		domainData[0] = this->E0_;
		domainData[1] = this->E1_;
		domainData[2] = this->poissonRatio_;
		domainData[3] = this->c1_;
		domainData[4] = this->D0_;
		domainData[5] = this->m_;

		double timeIncrement = this->getTimeIncrement();

		int integrationCode = 19;

		AceGenInterface::DeformationDiffusionNeoHookTetrahedra3D10 neoHookeElement(&positions[0], &displacements[0], &concentrations[0], &concentrationsLastConverged[0], &domainData[0], timeIncrement, integrationCode);
		neoHookeElement.computeTangentResidual();
		
		double *residuum = neoHookeElement.getResiduum();

		for (int i = 0; i < 40; i++)
			(*this->rhsVec_)[i] = -residuum[i];


#endif
	
	}


	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_NH<SC, LO, GO, NO>::assembleDeformationDiffusionNeoHook(SmallMatrixPtr_Type &elementMatrix)
	{
		std::vector<double> positions(30);
#ifdef FEDD_HAVE_ACEGENINTERFACE

		int count = 0;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}

		std::vector<double> displacements(30);
		for (int i = 0; i < 30; i++)
			displacements[i] = (*this->solution_)[i];

		std::vector<double> concentrations(10);
		for (int i = 0; i < 10; i++)
			concentrations[i] = (*this->solution_)[i + 30];

		std::vector<double> concentrationsLastConverged(10);
		for (int i = 0; i < 10; i++)
			concentrationsLastConverged[i] = (solution_n_)[i + 30];

		std::vector<double> domainData(6);
		domainData[0] = this->E0_;
		domainData[1] = this->E1_;
		domainData[2] = this->poissonRatio_;
		domainData[3] = this->c1_;
		domainData[4] = this->D0_;
		domainData[5] = this->m_;

		double timeIncrement = this->getTimeIncrement();

		int integrationCode = 19;

		AceGenInterface::DeformationDiffusionNeoHookTetrahedra3D10 neoHookeElement(&positions[0], &displacements[0], &concentrations[0], &concentrationsLastConverged[0], &domainData[0], timeIncrement, integrationCode);
		neoHookeElement.computeTangentResidual();

		double **stiffnessMatrix = neoHookeElement.getStiffnessMatrix();
		

		for (UN i = 0; i < this->dofsElement_; i++)
		{
			for (UN j = 0; j < this->dofsElement_; j++)
			{
				(*elementMatrix)[i][j] = -stiffnessMatrix[i][j];
			}
		}

	
	#endif
	}
} // namespace FEDD
#endif // ASSEMBLEFE_SCI_NH_DEF_hpp
