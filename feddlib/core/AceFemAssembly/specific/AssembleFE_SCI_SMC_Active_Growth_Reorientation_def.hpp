#ifndef AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
#define AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp

#include "AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"

#include <vector>
#include <string>
// #include <iostream>

namespace FEDD
{

	template <class SC, class LO, class GO, class NO>
	AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::AssembleFE_SCI_SMC_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple),
																																																						   timeParametersVecActive_(0),
																																																						   timeParametersVecGrowth_(0),
																																																						   timeParametersVecReorientation_(0)
	{
		activeAcceleratedEndTime_ = this->params_->sublist("Parameter").get("Accelerated Active Until",0.);
		TEUCHOS_TEST_FOR_EXCEPTION(activeAcceleratedEndTime_ < 0, std::logic_error, "!!! Warning: Accelerated Active Until not set correctly. Please Check Parameterlist !!!");

		activeAcceleratedMultiplier_ = this->params_->sublist("Parameter").get("Active Acceleration Multiplier",0.);
		TEUCHOS_TEST_FOR_EXCEPTION(activeAcceleratedMultiplier_ < 0, std::logic_error, "!!! Warning: Active Acceleration Multiplier not set correctly. Please Check Parameterlist !!!");
		
		int numMaterials = this->params_->sublist("Parameter Solid").get("Number of Materials", 1);
		int materialID = 0;
		for (int i = 1; i <= numMaterials; i++)
			if (this->params_->sublist("Parameter Solid").sublist(std::to_string(i)).get("Volume Flag", 15) == this->flag_)
				materialID = i;

		TEUCHOS_TEST_FOR_EXCEPTION(materialID == 0, std::logic_error, "!!! Warning: No corresponding parameterslist for the element flag = " << this->flag_ << ". Please Check volume flags of elements and Mesh Data !!!");

		this->iCode_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Integration Code", 18); // Only works for 18 currently!! Why? and is it still true?

#ifdef FEDD_HAVE_ACEGENINTERFACE

		// this->element_ = AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10(this->iCode_);
		AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 tempElem = AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10(this->iCode_);
		this->historyLength_ = tempElem.getHistoryLength();
		this->numberOfIntegrationPoints_ = tempElem.getNumberOfGaussPoints();
		this->postDataLength_ = tempElem.getNumberOfPostData();
		this->domainDataLength_ = tempElem.getNumberOfDomainData();
		char **domainDataNames = tempElem.getDomainDataNames();
		char **postDataNames = tempElem.getPostDataNames();
		

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
			this->domainDataNames_[i] = subString[i];
			this->domainData_[i] = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get(subString[i], 1.0);
			// std::cout << " DomainDataNames_ " << i << " "  << this->domainDataNames_[i] << " with value " << this->domainData_[i] << std::endl;

			TEUCHOS_TEST_FOR_EXCEPTION(this->domainData_[i] > 1.e12, std::logic_error, " Parameter not set correctly. Parameter " << this->domainDataNames_[i] << " received default value!!");

			// Pre-compute which parameters need acceleration/deceleration (optimization)
			if (domainDataNames_[i].find("LambdaBarCDotMax") != std::string::npos || 
			    domainDataNames_[i].find("LambdaBarCDotMin") != std::string::npos || 
			    domainDataNames_[i].find("Eta") != std::string::npos || 
			    domainDataNames_[i].find("K3") != std::string::npos || 
			    domainDataNames_[i].find("K4") != std::string::npos || 
			    domainDataNames_[i].find("K7") != std::string::npos || 
			    domainDataNames_[i].find("Beta1") != std::string::npos || 
			    domainDataNames_[i].find("Gamma6") != std::string::npos || 
			    domainDataNames_[i].find("KDotMin") != std::string::npos || 
			    domainDataNames_[i].find("KDotMax") != std::string::npos || 
			    domainDataNames_[i].find("LambdaBarDotPMin") != std::string::npos || 
			    domainDataNames_[i].find("LambdaBarDotPMax") != std::string::npos) {
				acceleratedParamIndices_.push_back(i);
			}
			else if (domainDataNames_[i].find("Gamma5") != std::string::npos || 
			         domainDataNames_[i].find("Gamma2") != std::string::npos) {
				deceleratedParamIndices_.push_back(i);
			}
		}

		for (int i = 0; i < this->postDataLength_; i++)
		{
			this->postDataNames_[i] = std::string(postDataNames[i]);
		}

		// Create map from field name to position in postDataNames_
    	for (int i = 0; i < postDataNames_.size(); i++) {
        	fieldNameToPosition_[postDataNames_[i]] = i;
    	}

		this->residuumRint_.resize(30, 0.0);
		this->residuumRc_.resize(10, 0.0);
		this->residuumRdyn_.resize(30, 0.0);

		this->stiffnessMatrixKuu_.resize(30, vec_dbl_Type(30, 0.0));
		this->stiffnessMatrixKuc_.resize(30, vec_dbl_Type(10, 0.0));
		this->stiffnessMatrixKcu_.resize(10, vec_dbl_Type(30, 0.0));
		this->stiffnessMatrixKcc_.resize(10, vec_dbl_Type(10, 0.0));
		this->massMatrixMc_.resize(10, vec_dbl_Type(10, 0.0));

#endif

		this->subiterationTolerance_ = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("Subiteration Tolerance", 1.e-7);
		this->FEType_ = std::get<1>(this->diskTuple_->at(0));	 // FEType of Disk
		this->dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
		this->dofsChem_ = std::get<2>(this->diskTuple_->at(1));	 // Degrees of freedom per node

		this->numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
		this->numNodesChem_ = std::get<3>(this->diskTuple_->at(1));	 // Number of nodes of element

		this->dofsElement_ = this->dofsSolid_ * this->numNodesSolid_ + this->dofsChem_ * this->numNodesChem_; // "Dimension of return matrix"

		this->positions_ = std::vector<double>(30, 0.0);
		this->displacements_ = std::vector<double>(30, 0.0);
		this->accelerations_ = std::vector<double>(30, 0.0);
		this->concentrations_ = std::vector<double>(10, 0.0);
		this->rates_ = std::vector<double>(10, 0.0);

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
		// this->element_.setElementID(this->getGlobalElementID());

		// Nodal Positions in Reference Coordinates
		int count = 0;
		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 3; j++)
			{
				this->positions_[count] = this->getNodesRefConfig()[i][j];
				count++;
			}

		// this->element_.setPositions(this->positions_.data());

		// Domain Data
		// this->element_.setDomainData(this->domainData_.data());

		// History Variables - initial values
		// this->element_.setHistoryVector(this->history_.data());

		// this->element_.setSubIterationTolerance(this->subiterationTolerance_);

		// this->element_.setComputeCompleted(false);

		// -----------
		// Active and Growth Time intervals
		int numSegmentsActive = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").get("Number of Segments", 0);
		int numSegmentsGrowth = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").get("Number of Segments", 0);
		int numSegmentsReorientation = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").get("Number of Segments", 0);

		// timeParametersVecActive[0]_.resize((0,vec_dbl_Type(2)));
		// timeParametersVecGrowth[0]_.resize((0,vec_dbl_Type(2)));

		for (int i = 1; i <= numSegmentsActive; i++)
		{

			double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").sublist(std::to_string(i)).get("Start Time", 0.);
			double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Active").sublist(std::to_string(i)).get("End Time", 0.);

			if (i == 1)
				TEUCHOS_TEST_FOR_EXCEPTION(startTime != this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).get("ActiveStartTime", 0.), std::logic_error, "!!!WARNING:: The ActiveStartTime and the start of Time stepping intervalls for active response do not match!!!");

			vec_dbl_Type segment = {startTime, endTime};
			timeParametersVecActive_.push_back(segment);

			// cout << " Active Segment " << i << ":[" << startTime << "," << endTime << "]" << endl;
		}

		for (int i = 1; i <= numSegmentsGrowth; i++)
		{

			double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").sublist(std::to_string(i)).get("Start Time", 0.);
			double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Growth").sublist(std::to_string(i)).get("End Time", 0.);

			vec_dbl_Type segment = {startTime, endTime};
			timeParametersVecGrowth_.push_back(segment);
			// cout << " Growth Segment " << i << ":[" << startTime << "," << endTime << "]" << endl;
		}

		for (int i = 1; i <= numSegmentsReorientation; i++)
		{

			double startTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").sublist(std::to_string(i)).get("Start Time", 0.);
			double endTime = this->params_->sublist("Parameter Solid").sublist(std::to_string(materialID)).sublist("Timestepping Intervalls Reorientation").sublist(std::to_string(i)).get("End Time", 0.);

			vec_dbl_Type segment = {startTime, endTime};
			timeParametersVecReorientation_.push_back(segment);
			// cout << " Reorientation Segment " << i << ":[" << startTime << "," << endTime << "]" << endl;
		}
		// if(this->globalElementID_<10)
		// 	cout << " Element initialized with timestep " << this->timeStep_ << endl;
		// We check this here also in case we do a restart.

	
#endif
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleJacobian()
	{
		SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(this->dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

		assemble_SCI_SMC_Active_Growth_Reorientation(true);

		for (int i = 0; i < 30; i++)
			for (int j = 0; j < 30; j++)
				(*elementMatrix)[i][j] = this->stiffnessMatrixKuu_[i][j];

		for (int i = 0; i < 30; i++)
			for (int j = 0; j < 10; j++)
				(*elementMatrix)[i][j + 30] = this->stiffnessMatrixKuc_[i][j];

		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 30; j++)
				(*elementMatrix)[i + 30][j] = this->stiffnessMatrixKcu_[i][j];

		for (int i = 0; i < 10; i++)
			for (int j = 0; j < 10; j++)
				(*elementMatrix)[i + 30][j + 30] = this->stiffnessMatrixKcc_[i][j] + (1. / this->getTimeIncrement()) * this->massMatrixMc_[i][j];

#endif

		this->jacobian_ = elementMatrix;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::checkingReorientationActiveGrowth()
	{
		bool restart = this->params_->sublist("Timestepping Parameter").get("Restart",false);
		double timeStepRestart = this->params_->sublist("Timestepping Parameter").get("Time step",0.0);
		bool found = false;
		for (int i = 0; i < timeParametersVecActive_.size(); i++)
		{
			if (this->timeStep_ + 1.0e-12 > timeParametersVecActive_[i][0] && this->timeStep_ - 1.0e-12 < timeParametersVecActive_[i][1])
			{
				this->activeBool_ = 1;
				if (!this->activeInitialized_)
				{
					// We need to ensure that in case of restart the is not initialized multiple times
					if (restart && timeStepRestart > timeParametersVecActive_[0][0] + 1.0e-10){  // only if the start of Active intervall coincides with the restart time
						this->activeInitialized_=true;
						this->growthInitialized_=true;
					}
					else{// Otherwise we know it was already previously initiated
						this->initializeActiveResponse();
						this->initializeGrowth();
					}

				}
				found = true;
			}
			else if (!found)
				this->activeBool_ = 0;
		}

		// Checking for growth
		found = false;
		for (int i = 0; i < timeParametersVecGrowth_.size(); i++)
		{
			if (this->timeStep_ + 1.0e-12 > timeParametersVecGrowth_[i][0] && this->timeStep_ - 1.0e-12 < timeParametersVecGrowth_[i][1])
			{
				this->growthBool_ = 1;
				if (!growthInitialized_){
					// We need to ensure that in case of restart the is not initialized multiple times
					if (restart && timeStepRestart > timeParametersVecActive_[0][0] + 1.0e-10){  // only if the start of Active intervall coincides with the restart time
						this->growthInitialized_=true;
					}
					else{ // Otherwise we know it was already previously initiated
						this->initializeGrowth();
					}
				}
				found = true;
			}
			else if (!found)
				growthBool_ = 0;
		}

		found = false;
		for (int i = 0; i < timeParametersVecReorientation_.size(); i++)
		{
			if (this->timeStep_ + 1.0e-12 > timeParametersVecReorientation_[i][0] && this->timeStep_ - 1.0e-12 < timeParametersVecReorientation_[i][1])
			{

				this->reorientationBool_ = 1;
				found = true;
			}
			else if (!found)
				this->reorientationBool_ = 0;
		}

		if (activeBool_ == 1 && growthBool_ == 1)
			std::cout << " WARNING: GROWTH AND ACTIVE RESPONSE HAPPENING SIMULTANEOUSLY!!" << std::endl;

		this->updateDomainData("growthBool", growthBool_);
		this->updateDomainData("ActiveBool", activeBool_);
		this->updateDomainData("reorientationBool", reorientationBool_);

		if (this->globalElementID_ == 0)
		{
			std::cout << " ---------------------------------------------- " << std::endl;
			std::cout << " CheckingReorientationActiveGrowth in elements" << std::endl;
			std::cout << " Growth " << growthBool_ << std::endl;
			std::cout << " Active " << activeBool_ << std::endl;
			std::cout << " Reorientation " << reorientationBool_ << std::endl;
			std::cout << " ---------------------------------------------- " << std::endl;
		}
	}
	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::advanceInTime(double dt)
	{
		// If we have a time segment setting we switch to the demanded time increment
		/*for(int i=0; i<numSegments_ ; i++){
			if(this->timeStep_+1.0e-12 > timeParametersVec_[i][0])
				this->timeIncrement_=timeParametersVec_[i][1];
		}*/
		// if (this->timeStep_ - 1.e-13 < 0) // only in this one instance T=0 we set the dt beforehand, as the initial dt is set through the paramterlist and this is error prone
		// this->timeIncrement_ = dt;

		if(this->timeStep_ - 1.e-13 < 0 && dt > 1.e-13) 
		{
			vec2D_dbl_Type gaussFibers = computeElementGaussFibersP2Tet();
			for(int i=0; i<gaussFibers.size(); i++)
			{
				this->historyUpdated_[39 + i*this->history_.size()/4] = gaussFibers[i][0]; // a11
				this->historyUpdated_[40 + i*this->history_.size()/4] = gaussFibers[i][1]; // a12
				this->historyUpdated_[41 + i*this->history_.size()/4] = gaussFibers[i][2]; // a13
				this->historyUpdated_[42 + i*this->history_.size()/4] = gaussFibers[i][3]; // a21
				this->historyUpdated_[43 + i*this->history_.size()/4] = gaussFibers[i][4]; // a22
				this->historyUpdated_[44 + i*this->history_.size()/4] = gaussFibers[i][5]; // a23
			}
			// 39 -- "a11"
			// 40 -- "a12"
			// 41 -- "a13"
			// 42 -- "a21"
			// 43 -- "a22"
			// 44 -- "a23"
			if(this->globalElementID_ <10)
			{
				std::cout << " Initial fiber direction in element " << this->globalElementID_ << ": " << std::endl;
				for(int i=0; i<gaussFibers.size(); i++)
				{
					std::cout << " Gauss Point " << i << ": a11: " << this->historyUpdated_[39 + i*this->history_.size()/4] << " a12: " << this->historyUpdated_[40 + i*this->history_.size()/4] << " a13: " << this->historyUpdated_[41 + i*this->history_.size()/4] << " a21: " << this->historyUpdated_[42 + i*this->history_.size()/4] << " a22: " << this->historyUpdated_[43 + i*this->history_.size()/4] << " a23: " << this->historyUpdated_[44 + i*this->history_.size()/4] << std::endl;

					std::cout << " Gauss Point before " << i << ": a11: " << this->history_[39 + i*this->history_.size()/4] << " a12: " << this->history_[40 + i*this->history_.size()/4] << " a13: " << this->history_[41 + i*this->history_.size()/4] << " a21: " << this->history_[42 + i*this->history_.size()/4] << " a22: " << this->history_[43 + i*this->history_.size()/4] << " a23: " << this->history_[44 + i*this->history_.size()/4] << std::endl;
				}
			}
		}
		this->timeStep_ = this->timeStep_ + this->timeIncrement_;

		this->timeIncrement_ = dt;

		// Checking for Active Response Reorientation and Growth
		checkingReorientationActiveGrowth();
		
		if (this->globalElementID_ == 0)
		{
			std::cout << " ---------------------------------------------- " << std::endl;
			std::cout << " AssembleFE_SCI_SMC: Advancing time in elements" << std::endl;
			std::cout << " Timestep: " << this->timeStep_ << " \t timeincrement: " << this->timeIncrement_ << " \t to be computed time: " << this->timeIncrement_ + this->timeStep_ << std::endl;
			std::cout << " ---------------------------------------------- " << std::endl;
		}
		// cout << " Update:: History " ;
		for (int i = 0; i < this->historyLength_; i++)
		{
			// if(this->timeStep_  > activeStartTime_ +dt )
			this->history_[i] = this->historyUpdated_[i];
			//  cout << " | " << this->history_[i] ;
		}
		// cout << endl;
		for (int i = 0; i < 10; i++)
			this->solutionC_n_[i] = (*this->solution_)[i + 30]; // this is the LAST solution of newton iterations
																// #ifdef FEDD_HAVE_ACEGENINTERFACE
																//		this->element_.setComputeCompleted(false);
																// #endif
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assembleRHS()
	{

		this->rhsVec_.reset(new vec_dbl_Type(this->dofsElement_, 0.));
#ifdef FEDD_HAVE_ACEGENINTERFACE

		assemble_SCI_SMC_Active_Growth_Reorientation(false);

		for (int i = 0; i < 30; i++)
			(*this->rhsVec_)[i] = this->residuumRint_[i]; //+residuumRDyn[i];

		for (int i = 0; i < 10; i++)
			(*this->rhsVec_)[i + 30] = this->residuumRc_[i];
#endif
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::assemble_SCI_SMC_Active_Growth_Reorientation(bool computeTangent)
	{
#ifdef FEDD_HAVE_ACEGENINTERFACE
	
		
		double deltaT = this->getTimeIncrement();
		// this->element_.setTimeIncrement(deltaT);

		double time = this->getTimeStep() + deltaT;
		// this->element_.setTime(time);

		for (int i = 0; i < 30; i++)
			this->displacements_[i] = (*this->solution_)[i];
		// this->element_.setDisplacements(this->displacements_.data());

		for (int i = 0; i < 10; i++)
		{
			this->concentrations_[i] = (*this->solution_)[i + 30];
			solutionC_n1_[i] = (*this->solution_)[i + 30]; // in each newtonstep solution for n+1 is updated.
		}
		// this->element_.setConcentrations(this->concentrations_.data());

		// this->element_.setAccelerations(this->accelerations_.data());

		for (int i = 0; i < 10; i++)
			this->rates_[i] = (this->solutionC_n1_[i] - this->solutionC_n_[i]) / deltaT;
		// this->element_.setRates(this->rates_.data());

		// this->element_.setHistoryVector(this->history_.data());
		std::vector<double> domainDataModified(this->domainDataLength_);

		// Copy base domain data
		std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

		// Apply modifications only if needed (using pre-computed indices - much faster!)
		if(time < this->activeAcceleratedEndTime_) {
			for(int idx : acceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
			}
			for(int idx : deceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
			}
		}

		AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem = AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10(this->positions_.data(), this->displacements_.data(), this->concentrations_.data(), this->accelerations_.data(), this->rates_.data(), domainDataModified.data(), this->history_.data(), this->subiterationTolerance_, deltaT, time, this->iCode_, this->getGlobalElementID());

		// std::cout << "elem.compute starts" << std::endl;

		int errorCode = elem.compute(computeTangent);

		// std::cout << "elem.compute ends" << std::endl;

		// vec_dbl_Type residuumRint_;
		// vec_dbl_Type residuumRdyn_;
		// vec_dbl_Type residuumRc_;

		// vec2D_dbl_Type stiffnessMatrixKuu_;
		// vec2D_dbl_Type stiffnessMatrixKuc_;
		// vec2D_dbl_Type stiffnessMatrixKcu_;
		// vec2D_dbl_Type stiffnessMatrixKcc_;
		// vec2D_dbl_Type massMatrixMc_;

		double *residuumRint = elem.getResiduumVectorRint();
		for (int i = 0; i < 30; i++)
			this->residuumRint_[i] = residuumRint[i];

		double *residuumRdyn = elem.getResiduumVectorRdyn();
		for (int i = 0; i < 30; i++)
			this->residuumRdyn_[i] = residuumRdyn[i];

		double *residuumRc = elem.getResiduumVectorRc();
		for (int i = 0; i < 10; i++)
			this->residuumRc_[i] = residuumRc[i];

		if(computeTangent){
			double **stiffnessMatrixKuu = elem.getStiffnessMatrixKuu();
			for (int i = 0; i < 30; i++)
				for (int j = 0; j < 30; j++)
				{
					this->stiffnessMatrixKuu_[i][j] = stiffnessMatrixKuu[i][j];
				}

			double **stiffnessMatrixKuc = elem.getStiffnessMatrixKuc();
			for (int i = 0; i < 30; i++)
				for (int j = 0; j < 10; j++)
					this->stiffnessMatrixKuc_[i][j] = stiffnessMatrixKuc[i][j];

			double **stiffnessMatrixKcu = elem.getStiffnessMatrixKcu();
			for (int i = 0; i < 10; i++)
				for (int j = 0; j < 30; j++)
					this->stiffnessMatrixKcu_[i][j] = stiffnessMatrixKcu[i][j];

			double **massMatrixMc = elem.getMassMatrixMc();
			for (int i = 0; i < 10; i++)
				for (int j = 0; j < 10; j++)
					this->massMatrixMc_[i][j] = massMatrixMc[i][j];

			double **stiffnessMatrixKcc = elem.getStiffnessMatrixKcc();
			for (int i = 0; i < 10; i++)
				for (int j = 0; j < 10; j++)
					this->stiffnessMatrixKcc_[i][j] = stiffnessMatrixKcc[i][j];

			double *historyUpdated = elem.getHistoryUpdated();
			for (int i = 0; i < this->historyLength_; i++)
				this->historyUpdated_[i] = historyUpdated[i];
		}

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

		double deltaT = this->getTimeIncrement();

		double time = this->getTimeStep() + deltaT;

		// if(this->growthInitialized_==true)
		// {
		// 	for(int i=0;i<this->historyLength_;i++)
		// 		std::cout << this->history_[i] << " ";
		// 	std::cout << std::endl;
		// } -- This seems to be okay

		std::vector<double> domainDataModified(this->domainDataLength_);
		
		// Copy base domain data
		std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

		// Apply modifications only if needed (using pre-computed indices - much faster!)
		if(time < this->activeAcceleratedEndTime_) {
			for(int idx : acceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
			}
			for(int idx : deceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
			}
		}

		AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(this->positions_.data(), &displacements[0], &concentrations[0], &accelerations[0], &rates[0], domainDataModified.data(), this->history_.data(), this->subiterationTolerance_, deltaT, time, this->iCode_, this->getGlobalElementID());

		// std::cout << "History values going into PP: " << std::endl;
		// std::cout << "Agn1: { " << this->history_[19] << ", " << this->history_[20] << ", " << this->history_[21] << " }" << std::endl;
		// std::cout << "Agn2: { " << this->history_[22] << ", " << this->history_[23] << ", " << this->history_[24] << " }" << std::endl;
		// std::cout << "Agn3: { " << this->history_[25] << ", " << this->history_[26] << ", " << this->history_[27] << " }" << std::endl;

		// elem.compute(); --THis shit right here could be the problem

		double **postProcessingResults = elem.postProcess(&displacements[0], &concentrations[0], this->history_.data(), &rates[0], &accelerations[0]); //Inside this the values are 0

		for (int i = 0; i < 10; i++)
		{
			//cout << " Node " << i << " ";
			for (int j = 0; j < this->postDataLength_; j++)
			{
				(*this->postProcessingData_)[i][j] = postProcessingResults[i][j];
				// cout << postProcessingResults[i][j] << " " ;
			}
			// if (this->growthInitialized_ == true && i == 0)
			// {
			// 	cout << "Mises Stress: " << postProcessingResults[i][10] << " ";
			// 	cout << "SCirc: " << postProcessingResults[i][11] << " ";
			// 	cout << "SAxial: " << postProcessingResults[i][12] << " ";
			// 	cout << "SRadial: " << postProcessingResults[i][13] << " ";
			// 	// 30 - 38      "Ag1n1","Ag1n2","Ag1n3","Ag2n1","Ag2n2","Ag2n3",
			// 	//            "Ag3n1","Ag3n2","Ag3n3"
			// 	cout << "Ag1n1: " << postProcessingResults[i][30] << " ";
			// 	cout << "Ag1n2: " << postProcessingResults[i][31] << " ";
			// 	cout << "Ag1n3: " << postProcessingResults[i][32] << " ";
			// 	cout << "Ag2n1: " << postProcessingResults[i][33] << " ";
			// 	cout << "Ag2n2: " << postProcessingResults[i][34] << " ";
			// 	cout << "Ag2n3: " << postProcessingResults[i][35] << " ";
			// 	cout << "Ag3n1: " << postProcessingResults[i][36] << " ";
			// 	cout << "Ag3n2: " << postProcessingResults[i][37] << " ";
			// 	cout << "Ag3n3: " << postProcessingResults[i][38] << " ";
			// 	cout << endl;
			// }
		}
#endif
	}

	template <class SC, class LO, class GO, class NO>
	int AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::findPosition(std::string subString, std::vector<std::string> &stringArray)
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

		double deltaT = this->getTimeIncrement();

		double time = this->getTimeStep() + deltaT;

		std::vector<double> domainDataModified(this->domainDataLength_);
		
		// Copy base domain data
		std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

		// Apply modifications only if needed (using pre-computed indices - much faster!)
		if(time < this->activeAcceleratedEndTime_) {
			for(int idx : acceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
			}
			for(int idx : deceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
			}
		}

		AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(this->positions_.data(), &displacements[0], &concentrations[0], &accelerations[0], &rates[0], domainDataModified.data(), this->history_.data(), this->subiterationTolerance_, deltaT, time, this->iCode_, this->getGlobalElementID());

		std::vector<double> historyNew = elem.initializeGrowthOrientationVectors();
		// std::cout << "Growth Orientation Vectors being set! \n HistoryOld: \n";
		// for (int i = 0; i < this->historyLength_; i++)
		// 	std::cout << this->history_[i] << " ";
		// std::cout << "\n HistoryNew: \n";
		// for (int i = 0; i < historyNew.size(); i++)
		// 	std::cout << historyNew[i] << " ";
		for (int i = 0; i < this->historyLength_; i++){
			this->history_[i] = historyNew[i];
			this->historyUpdated_[i] = historyNew[i];
		}
#endif
		growthInitialized_ = true;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::initializeActiveResponse()
	{
		double deltaT = this->getTimeIncrement();
		std::cout << " Initialize active Response " << std::endl;
		double time = this->getTimeStep() + deltaT;
#ifdef FEDD_HAVE_ACEGENINTERFACE
		std::vector<double> domainDataModified(this->domainDataLength_);
		
		// Copy base domain data
		std::copy(this->domainData_.begin(), this->domainData_.end(), domainDataModified.begin());

		// Apply modifications only if needed (using pre-computed indices - much faster!)
		if(time < this->activeAcceleratedEndTime_) {
			for(int idx : acceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] * this->activeAcceleratedMultiplier_;
			}
			for(int idx : deceleratedParamIndices_) {
				domainDataModified[idx] = this->domainData_[idx] / this->activeAcceleratedMultiplier_;
			}
		}
		
		AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(this->positions_.data(), this->displacements_.data(), this->concentrations_.data(), this->accelerations_.data(), this->rates_.data(), domainDataModified.data(), this->history_.data(), this->subiterationTolerance_, deltaT, time, this->iCode_, this->getGlobalElementID());

		std::vector<double> stretches = elem.getGaussPointStretches();
		std::cout << " Streches: ";
		for (int i = 0; i < stretches.size(); i++)
			std::cout << stretches[i] << " ";
		std::cout << std::endl;

		int historyPerGP = (int)this->historyLength_ / this->numberOfIntegrationPoints_;
		for (int i = 0; i < this->numberOfIntegrationPoints_; i++)
		{
			this->history_[i * historyPerGP + 10] = stretches[i * 2];
			this->history_[i * historyPerGP + 11] = stretches[i * 2 + 1];
			this->historyUpdated_[i * historyPerGP + 10] = stretches[i * 2];
			this->historyUpdated_[i * historyPerGP + 11] = stretches[i * 2 + 1];
		}
#endif
		activeInitialized_ = true;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::updateDomainData(std::string dataName, double dataValue)
	{
		int position = findPosition(dataName, this->domainDataNames_);
		// cout << " !!!! Position " << position << " of " << dataName << " found " << endl;
		this->domainData_[position] = dataValue;
	}

	template <class SC, class LO, class GO, class NO>
	vec2D_dbl_Type AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::computeElementGaussFibersP2Tet()
	{
		// Xe[a][i] = global coordinate i of node a
		// Xe size should be 10 x 3

		const double a = 0.5854101966249685;
		const double b = 0.1381966011250105;

		std::vector<std::vector<double>> gp = {
			{b, b, b},
			{a, b, b},
			{b, a, b},
			{b, b, a}
		};

		vec2D_dbl_Type fibers(4,vec_dbl_Type(6,0));

		for (int g = 0; g < 4; ++g)
		{
			double xi   = gp[g][0];
			double eta  = gp[g][1];
			double zeta = gp[g][2];

			double L1 = 1.0 - xi - eta - zeta;
			double L2 = xi;
			double L3 = eta;
			double L4 = zeta;

			std::vector<double> N(10);

			N[0] = L1 * (2.0*L1 - 1.0);
			N[1] = L2 * (2.0*L2 - 1.0);
			N[2] = L3 * (2.0*L3 - 1.0);
			N[3] = L4 * (2.0*L4 - 1.0);

			N[4] = 4.0 * L1 * L2;
			N[5] = 4.0 * L2 * L3;
			N[6] = 4.0 * L1 * L3;
			N[7] = 4.0 * L1 * L4;
			N[8] = 4.0 * L2 * L4;
			N[9] = 4.0 * L3 * L4;

			double x1 = 0.0;
			double x2 = 0.0;
			double x3 = 0.0;

			for (int A = 0; A < 10; ++A)
			{
				x1 += N[A] * this->getNodesRefConfig()[A][0];
				x2 += N[A] * this->getNodesRefConfig()[A][1];
				x3 += N[A] * this->getNodesRefConfig()[A][2];
			}

			fibers[g] = computeFiberDirections(x1, x2, x3);
		}

		return fibers;
	}

	template <class SC, class LO, class GO, class NO>
	vec_dbl_Type AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC, LO, GO, NO>::computeFiberDirections(double x1, double x2, double x3)
	{
		double pi = 3.14159265358979323846;
		double beta = fA_ * pi / 180.0;

		const double cb = std::cos(beta);
		const double sb = std::sin(beta);

		// ------------------------------------------------------------
		// Compute alpha (same logic as FEAP code)
		// ------------------------------------------------------------

		double alpha;

		if (x1 < 1.0) {
			alpha = -std::atan(std::abs(x3) / std::abs(1.0 - x1));
		}
		else {
			alpha = -0.5 * pi;
		}

		// ------------------------------------------------------------
		// Shift point by center (1,0,0)
		// ------------------------------------------------------------

		const double xs = x1 - 1.0;
		const double ys = x2;
		const double zs = x3;

		// ------------------------------------------------------------
		// rotate2643 convention
		// ROTATION_SIGN = -1
		// Rotate point by alpha around y-axis
		// ------------------------------------------------------------

		const double ca = std::cos(alpha);
		const double sa = std::sin(alpha);

		double xt = ca*xs - sa*zs;
		double yt = ys;
		double zt = sa*xs + ca*zs;

		// shift back
		xt += 1.0;

		// ------------------------------------------------------------
		// Radius in rotated frame
		// ------------------------------------------------------------

		double r = std::sqrt(xt*xt + yt*yt);

		if (r < 1.0e-12)
			r = 1.0e-12;

		// ------------------------------------------------------------
		// Local fibre directions
		// ------------------------------------------------------------

		// +43°
		double a1x_local = -yt/r * cb;
		double a1y_local =  xt/r * cb;
		double a1z_local = -sb;

		// -43°
		double a2x_local = -yt/r * cb;
		double a2y_local =  xt/r * cb;
		double a2z_local =  sb;

		// ------------------------------------------------------------
		// Rotate fibres back by -alpha
		// ------------------------------------------------------------

		const double ca2 = std::cos(-alpha);
		const double sa2 = std::sin(-alpha);

		// fibre family 1
		double a11 = ca2*a1x_local - sa2*a1z_local;
		double a12 = a1y_local;
		double a13 = sa2*a1x_local + ca2*a1z_local;

		// fibre family 2
		double a21 = ca2*a2x_local - sa2*a2z_local;
		double a22 = a2y_local;
		double a23 = sa2*a2x_local + ca2*a2z_local;

		// ------------------------------------------------------------
		// Normalize
		// ------------------------------------------------------------

		double n1 = std::sqrt(a11*a11 + a12*a12 + a13*a13);
		double n2 = std::sqrt(a21*a21 + a22*a22 + a23*a23);

		if (n1 > 1.0e-12) {
			a11 /= n1;
			a12 /= n1;
			a13 /= n1;
		}

		if (n2 > 1.0e-12) {
			a21 /= n2;
			a22 /= n2;
			a23 /= n2;
		}

		// ------------------------------------------------------------
		// Output:
		//
		// v[0] = {a11,a12,a13}
		// v[1] = {a21,a22,a23}
		// ------------------------------------------------------------

		vec_dbl_Type v(6,0);

		v[0] = a11;
		v[1] = a12;
		v[2] = a13;

		v[3] = a21;
		v[4] = a22;
		v[5] = a23;

		return v;
	}


} // namespace FEDD
#endif // AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp