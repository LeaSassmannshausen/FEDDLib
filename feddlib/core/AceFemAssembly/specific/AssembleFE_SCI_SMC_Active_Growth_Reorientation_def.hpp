#ifndef AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
#define AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp

#include "AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"

#include <vector>
// #include <iostream>

namespace FEDD
{

	template <class SC, class LO, class GO, class NO>
	AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::AssembleFE_SCI_SMC_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple)
	{

		int numMaterials = this->params_->sublist("Parameter Solid").get("Number of Materials", 1);
		int materialID = 0;
		for (int i = 1; i <= numMaterials; i++)
			if (this->params_->sublist("Parameter Solid").sublist(std::to_string(i)).get("Volume Flag", 15) == this->flag_)
				materialID = i;

		TEUCHOS_TEST_FOR_EXCEPTION(materialID == 0, std::logic_error, "!!! Warning: No corresponding parameterslist for the element flag = " << this->flag_ << ". Please Check volume flags of elements and Mesh Data !!!");

		this->iCode_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Integration Code", 18); // Only works for 18 currently!! Why? and is it still true?

#ifdef FEDD_HAVE_ACEGENINTERFACE


		this->element_ = Teuchos::rcp(new AceGenElement_Type(this->iCode_));

		this->historyLength_ = this->element_->getHistoryLength();
		this->numberOfIntegrationPoints_ = this->element_->getNumberOfGaussPoints();
		this->postDataLength_ = this->element_->getNumberOfPostData();
		this->domainDataLength_ = this->element_->getNumberOfDomainData();
		char** domainDataNames = this->element_->getDomainDataNames();
		char** postDataNames = this->element_->getPostDataNames();

		std::vector<std::string> subString(this->domainDataLength_);

		this->domainDataNames_.resize(this->domainDataLength_);
		this->postDataNames_.resize(this->postDataLength_);
		this->domainData_.resize(this->domainDataLength_, 0.0);

		for (int i = 0; i < this->domainDataLength_; i++)
		{
			this->domainDataNames_[i] = std::string(domainDataNames[i]);
			int pos1 = domainDataNames_[i].find("-");
			int pos2 = domainDataNames_[i].find("_");

			subString[i] = domainDataNames_[i].substr(pos1 + 1, pos2 - pos1 - 1);
			this->domainData_[i] = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get(subString[i], 0.0);
		}

		for (int i = 0; i < this->postDataLength_; i++)
			this->postDataNames_[i] = std::string(postDataNames[i]);

#endif

		this->subiterationTolerance_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Subiteration Tolerance", 1.e-7);
		this->FEType_ = std::get<1>(this->diskTuple_->at(0));	 // FEType of Disk
		this->dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
		this->dofsChem_ = std::get<2>(this->diskTuple_->at(1));	 // Degrees of freedom per node

		this->numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
		this->numNodesChem_ = std::get<3>(this->diskTuple_->at(1));	 // Number of nodes of element

		this->dofsElement_ = this->dofsSolid_ * this->numNodesSolid_ + this->dofsChem_ * this->numNodesChem_; // "Dimension of return matrix"

		// Einlesen durch Parameterdatei irgendwann cool

		// historyGP: Vector of history variables [Order: LambdaBarC1, LambdaBarC2, nA1, nA2, nB1, nB2, nC1, nC2, nD1, nD2, LambdaA1, LambdaA2, k251, k252, LambdaBarP1, LambdaBarP2, Theta1, Theta2, Theta3, Ag11, Ag12, Ag13, Ag21, Ag22, Ag23, Ag31, Ag32, Ag33, a11, a12, a13, a21, a22, a23] (The length must be equal to number of history variables per gauss point(34) * number of gauss points)
		std::vector<double> historyGP = {1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.82758, 1.82758, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

		this->history_ = historyGP;

		this->history_.reserve(this->historyLength_);
		for (int i = 0; i < this->numberOfIntegrationPoints_ - 1; i++)
			this->history_.insert(this->history_.end(), historyGP.begin(), historyGP.end());

		// Error out if history length is inconsistent
		TEUCHOS_TEST_FOR_EXCEPTION(this->history_.size() != this->historyLength_, std::logic_error, "History input length does not match history size of model! \n Hisory input length: " << this->history_.size() << "\n History size of model: " << this->historyLength_ << "\n");

		this->historyUpdated_.resize(this->historyLength_, 0.);

		this->solutionC_n_.resize(10, 0.);
		this->solutionC_n1_.resize(10, 0.);

		this->postProcessingData_ = Teuchos::rcp(new vec2D_dbl_Type(this->numNodesSolid_, vec_dbl_Type(this->postDataLength_)));
		this->solution_.reset(new vec_dbl_Type(this->dofsElement_, 0.));

#ifdef FEDD_HAVE_ACEGENINTERFACE

		// Element ID
		this->element_->setElementID(this->getGlobalElementID());

		// Nodal Positions in Reference Coordinates
		positions_.resize(30);
		int count = 0;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				positions_[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}
		this->element_->setPositions(this->positions_.data());

		// Domain Data
		this->element_->setDomainData(this->domainData_.data());

		// History Variables - initial values
		this->element_->setHistoryVector(this->history_.data());

		this->element_->setSubIterationTolerance(this->subiterationTolerance_);

		this->element_->setComputeCompleted(false);
#endif

	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleJacobian()
	{

		SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(this->dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

		assemble_SCI_SMC_Active_Growth_Reorientation(); // Use this if nothing works!

		double **stiffnessMatrixKuu = this->element_->getStiffnessMatrixKuu();
		double **stiffnessMatrixKuc = this->element_->getStiffnessMatrixKuc();
		double **stiffnessMatrixKcu = this->element_->getStiffnessMatrixKcu();
		double **stiffnessMatrixKcc = this->element_->getStiffnessMatrixKcc();
		double **massMatrixMc = this->element_->getMassMatrixMc();

		for (int i = 0; i < 30; i++)
			for (int j = 0; j < 30; j++)
				(*elementMatrix)[i][j] = stiffnessMatrixKuu[i][j];

		for (int i = 0; i < 30; i++)
			for (int j = 0; j < 10; j++)
				(*elementMatrix)[i][j + 30] = stiffnessMatrixKuc[i][j];

		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 30; j++)
				(*elementMatrix)[i + 30][j] = stiffnessMatrixKcu[i][j];

		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 10; j++)
				(*elementMatrix)[i + 30][j + 30] = stiffnessMatrixKcc[i][j] + (1. / this->getTimeIncrement()) * massMatrixMc[i][j];

#endif

		this->jacobian_ = elementMatrix;
		
	}
	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::advanceInTime(double dt)
	{

		// If we have a time segment setting we switch to the demanded time increment
		/*for(int i=0; i<numSegments_ ; i++){
			if(this->timeStep_+1.0e-12 > timeParametersVec_[i][0])
				this->timeIncrement_=timeParametersVec_[i][1];
		}*/
		if (this->timeStep_ - 1.e-13 < 0) // only in this one instance T=0 we set the dt beforehand, as the initial dt is set through the paramterlist and this is error prone
			this->timeIncrement_ = dt;

		this->timeStep_ = this->timeStep_ + this->timeIncrement_;

		this->timeIncrement_ = dt;

		if (this->globalElementID_ == 0)
		{
			cout << " ---------------------------------------------- " << endl;
			cout << " AssembleFE_SCI_SMC: Advancing time in elements" << endl;
			cout << " Timestep: " << this->timeStep_ << " \t timeincrement: " << this->timeIncrement_ << endl;
			cout << " ---------------------------------------------- " << endl;
		}
		// cout << " Update:: History " ;
		for (int i = 0; i < this->historyLength_; i++)
		{
			// if(this->timeStep_  > activeStartTime_ +dt )
			this->history_[i] = this->historyUpdated_[i];
			// cout << " | " << this->history_[i] ;
		}
		// cout << endl;
		for (int i = 0; i < 10; i++)
			this->solutionC_n_[i] = (*this->solution_)[i + 30]; // this is the LAST solution of newton iterations
#ifdef FEDD_HAVE_ACEGENINTERFACE
		this->element_->setComputeCompleted(false);
#endif

	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleRHS()
	{

		this->rhsVec_.reset(new vec_dbl_Type(this->dofsElement_, 0.));

#ifdef FEDD_HAVE_ACEGENINTERFACE

		assemble_SCI_SMC_Active_Growth_Reorientation();

		double *residuumRint = this->element_->getResiduumVectorRint();

		double *residuumRDyn = this->element_->getResiduumVectorRdyn();

		for (int i = 0; i < 30; i++)
			(*this->rhsVec_)[i] = residuumRint[i]; //+residuumRDyn[i];

		double *residuumRc = this->element_->getResiduumVectorRc();

		for (int i = 0; i < 10; i++)
			(*this->rhsVec_)[i + 30] = residuumRc[i];
#endif
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assemble_SCI_SMC_Active_Growth_Reorientation()
	{
#ifdef FEDD_HAVE_ACEGENINTERFACE

		double deltaT = this->getTimeIncrement();
		this->element_->setTimeIncrement(deltaT);

		double time = this->getTimeStep() + deltaT;
		this->element_->setTime(time);

		double displacements[30];
		for (int i = 0; i < 30; i++)
			displacements[i] = (*this->solution_)[i];
		this->element_->setDisplacements(displacements);

		double concentrations[10];
		for (int i = 0; i < 10; i++)
		{
			concentrations[i] = (*this->solution_)[i + 30];
			solutionC_n1_[i] = (*this->solution_)[i + 30]; // in each newtonstep solution for n+1 is updated.
		}
		this->element_->setConcentrations(concentrations);

		double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		this->element_->setAccelerations(accelerations);

		double rates[10];
		for (int i = 0; i < 10; i++)
			rates[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / deltaT;
		this->element_->setRates(rates);

		this->element_->setHistoryVector(this->history_.data());

		// AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(positions, displacements, concentrations, accelerations, rates, &this->domainData_[0], &this->history[0], subIterationTolerance, deltaT, time, this->iCode_, this->getGlobalElementID());
		int errorCode = this->element_->compute();

		this->element_->setComputeCompleted(true);

		double *historyUpdated = this->element_->getHistoryUpdated();
		for (int i = 0; i < this->historyLength_; i++)
			this->historyUpdated_[i] = historyUpdated[i];
#endif
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::postProcessing()
	{

#ifdef FEDD_HAVE_ACEGENINTERFACE

		double displacements[30];
		for (int i = 0; i < 30; i++)
			displacements[i] = (*this->solution_)[i];

		double concentrations[10];
		for (int i = 0; i < 10; i++)
		{
			concentrations[i] = (*this->solution_)[i + 30];
			solutionC_n1_[i] = (*this->solution_)[i + 30]; // in each newtonstep solution for n+1 is updated.
		}

		double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double rates[10];
		for (int i = 0; i < 10; i++)
			rates[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / this->getTimeIncrement();

		double **postProcessingResults = this->element_->postProcess(&displacements[0], &concentrations[0], this->historyUpdated_.data(), &rates[0], &accelerations[0]);

		for (int i = 0; i < 10; i++){
			cout << " Node " << i << " " ;
			for (int j = 0; j < this->postDataLength_; j++){
				(*this->postProcessingData_)[i][j] = postProcessingResults[i][j];
				cout << postProcessingResults[i][j] << " " ;
			}
			cout << endl;
		}
#endif
	}

	template <class SC, class LO, class GO, class NO>
	int AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::findPosition(string subString, std::vector<string> &stringArray)
	{
		for (int i = 0; i < stringArray.size(); i++)
		{
			if (stringArray[i] == subString)
				return i;
		}
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::initializeGrowth()
	{
#ifdef FEDD_HAVE_ACEGENINTERFACE

		std::vector<double> historyNew = this->element_->initializeGrowthOrientationVectors();
		for (int i = 0; i < this->historyLength_; i++)
			this->history_[i] = historyNew[i];
#endif

	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::initializeActiveResponse()
	{
#ifdef FEDD_HAVE_ACEGENINTERFACE

		std::vector<double> stretches = this->element_->getGaussPointStretches();
		int historyPerGP = (int)this->historyLength_ / this->numberOfIntegrationPoints_;
		for (int i = 0; i < this->numberOfIntegrationPoints_; i++)
		{
			history_[i * historyPerGP + 10] = stretches[i * 2];
			history_[i * historyPerGP + 11] = stretches[i * 2 + 1];
		}
#endif

	}
	
	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::updateDomainData(std::string dataName, double dataValue)
	{
		int position = findPosition(dataName, this->domainDataNames_);
		this->domainData_[position] = dataValue;
		
#ifdef FEDD_HAVE_ACEGENINTERFACE

		this->element_->setDomainData(this->domainData_.data());
#endif

	}

} // namespace FEDD
#endif // AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
