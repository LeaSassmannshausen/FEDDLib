// ------------------------ Parameter ------------------------------
	// Stretch-dependent chemical kinetics model for smooth muscle cells
	// -----------------------------------------------------------------
	// Kinetics model
	k2_ = this->params_->sublist("Parameter Solid").get("K2",0.2e0); 
	k5_ = this->params_->sublist("Parameter Solid").get("K5",0.2e0);
	k3_ = this->params_->sublist("Parameter Solid").get("K3",0.134e0); // ??
	k4_ = this->params_->sublist("Parameter Solid").get("K4",0.166e-2); // ??
	k7_= this->params_->sublist("Parameter Solid").get("K7",0.66e-4); // ??
	// k_1 / k_6
	ca50_ = this->params_->sublist("Parameter Solid").get("Ca50",0.4e0); // ??
	// Ca^2+
	gamma1_ = this->params_->sublist("Parameter Solid").get("Gamma1",0.5131e0); 
	// \bar{lambda}_c determined by evolution equation with
	lambdaBarCDotMax_= this->params_->sublist("Parameter Solid").get("LambdaBarCDotMax",0.443e-1); // ??
	lambdaBarCDotMin_= this->params_->sublist("Parameter Solid").get("LambdaBarCDotMin",-0.443e-1); // ??
	gamma2_ = this->params_->sublist("Parameter Solid").get("Gamma2",50.0e0); // ??
	// and tau_c
	// the targe calcium concentration Ca^{2+}_tar is defined with
	gamma3_= this->params_->sublist("Parameter Solid").get("Gamma3",0.9e0);
	lambdaC50_ = this->params_->sublist("Parameter Solid").get("LambdaC50",0.12e1); // ??
	// The reduction of k_2 and k_5 due calcium sensitization is defined by two evolution equations
	zeta1_ = this->params_->sublist("Parameter Solid").get("Zeta1",100.e0); // penalty
	gamma4_ = this->params_->sublist("Parameter Solid").get("Gamma4",200.e0);
	gamma5_ = this->params_->sublist("Parameter Solid").get("Gamma5", 50.e0 );
	zeta2_ = this->params_->sublist("Parameter Solid").get("Zeta2", 1000.e0 );
	DeltaLambdaBarPMin_ =this->params_->sublist("Parameter Solid").get("DeltaLambdaBarPMin", -0.1e-4);
	kDotMin_ = this->params_->sublist("Parameter Solid").get("KDotMin",-0.10694e-1); // \dot{k}_2/5,min
	kDotMax_ = this->params_->sublist("Parameter Solid").get("KDotMax",0.9735e-3); // \dot{k}_2/5,max
	lambdaBarDotPMin_ = this->params_->sublist("Parameter Solid").get("LambdaBarDotMin",-0.2323e-3);
	lambdaBarDotPMax_ = this->params_->sublist("Parameter Solid").get("LambdaBarDotMax",0.699e-4);
	kMin_ = this->params_->sublist("Parameter Solid").get("KMin",0.15e0);
	// The target stretch dependent MLCP activity is defined with
	lambdaP50_ = this->params_->sublist("Parameter Solid").get("LambdaP50",1.0e0);
	gamma6_ = this->params_->sublist("Parameter Solid").get("Gamma6",0.15e1);
	// --------------
	// Effect of pharmacological agents
	c50_ = this->params_->sublist("Parameter Solid").get("C50",0.5e0);
	p1_ = this->params_->sublist("Parameter Solid").get("P1",0.6e0);
	p3_ = this->params_->sublist("Parameter Solid").get("P3",0.6e0);	
	// --------------
	// Smooth muscle cell activation
	muA_ = this->params_->sublist("Parameter Solid").get("MuA",0.11857e-1); 
	kappa_ = this->params_->sublist("Parameter Solid").get("Kappa",0.148262e0); 
	beta1_ = this->params_->sublist("Parameter Solid").get("Beta1",0.1006e-2); // ??
	beta2_ = this->params_->sublist("Parameter Solid").get("Beta2",0.2668e-1); 
	// At Starttime 1000 the diffused drug influences the material model. -> Active response at T=starttime	
	activeStartTime_ = this->params_->sublist("Parameter Solid").get("ActiveStartTime",1.0); 
	// ----------------------------------------------------------------------------
	// Diffusion Parameter
	d0_ = this->params_->sublist("Parameter Diffusion").get("D0",6.e-05);
	m_ = this->params_->sublist("Parameter Solid").get("m",0.e0);

	// -----------------------------------------------------------------
	// Growth and reorientation parameters (residual stresses and fiber reorientation)
	// -----------------------------------------------------------------
	growthStartTime_ = this->params_->sublist("Parameter Solid").get("GrowthStartTime",0.e0);
	reorientationStartTime_ = this->params_->sublist("Parameter Solid").get("ReorientationStartTime",0.e0);
	growthEndTime_ = this->params_->sublist("Parameter Solid").get("GrowthEndTime",0.e0);
	reorientationEndTime_ = this->params_->sublist("Parameter Solid").get("ReorientationEndTime",0.e0);
	// Evolution of the growth factors
	thetaPlus1_ = this->params_->sublist("Parameter Solid").get("ThetaPlus1",1.005063);
	thetaPlus2_ = this->params_->sublist("Parameter Solid").get("ThetaPlus2",1.020595);
	thetaPlus3_ = this->params_->sublist("Parameter Solid").get("ThetaPlus3",1.119918);
	thetaMinus1_ = this->params_->sublist("Parameter Solid").get("ThetaMinus1",0.98e0);
	thetaMinus2_ = this->params_->sublist("Parameter Solid").get("ThetaMinus2",0.98e0);
	thetaMinus3_ = this->params_->sublist("Parameter Solid").get("ThetaMinus3",0.98e0);
	kThetaPlus_ = this->params_->sublist("Parameter Solid").get("KThetaPlus",1.e-4);
	kThetaMinus_ = this->params_->sublist("Parameter Solid").get("KThetaMinus",1.e-4);
	mThetaPlus_ = this->params_->sublist("Parameter Solid").get("MThetaPlus",3.e0);
	mThetaMinus_ = this->params_->sublist("Parameter Solid").get("MThetaMinus",3.e0);
	// Evolution equation for the angle between fibre direction and the target fibre direction
	kEtaPlus_ = this->params_->sublist("Parameter Solid").get("KEtaPlus",0.1e-3);
	mEtaPlus_ = this->params_->sublist("Parameter Solid").get("MEtaPlus",5.0e0);
	eta_ = this->params_->sublist("Parameter Solid").get("Eta",0.18745e0); // ??
	// Fibre angle
	fA_= this->params_->sublist("Parameter Solid").get("FA",30.e0); // ??
	 
	// Strain energy density function ('The strain energy density stored within the arterial wall tissue is defined to be additively decoupled into a passive and an active part. wherein the passive response, is governed by an isotropic elastin–based ground substance and two transversely isotropic collagen fiber families.)
	alpha2_ = this->params_->sublist("Parameter Solid").get("Alpha2",0.15173775e0); 
	alpha3_ = this->params_->sublist("Parameter Solid").get("Alpha3",0.275662e1); 
	alpha1_ = this->params_->sublist("Parameter Solid").get("Alpha1",11.52507e-3);
	alpha4_ = this->params_->sublist("Parameter Solid").get("Alpha4",1.27631e-3);
	alpha5_ = this->params_->sublist("Parameter Solid").get("Alpha5",0.308798e1); 
	// Density
	rho_ = this->params_->sublist("Parameter Solid").get("Rho",1.e0);