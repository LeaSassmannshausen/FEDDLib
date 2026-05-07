#ifndef PrecOpFaCSI_DEF_hpp
#define PrecOpFaCSI_DEF_hpp
#include "PrecOpFaCSI_decl.hpp"
#include <Thyra_TpetraMultiVector_decl.hpp>
#include <Teuchos_VerboseObject.hpp>
/*!
 Definition of PrecOpFaCSI

 @brief  PrecOpFaCSI
 @author Christian Hochmuth
 @version 1.0
 @copyright CH

 The operator implemented in applyImpl() applies the FaCSI block
 preconditioner to a product vector. The standard block ordering is

   0: fluid velocity       (fv)
   1: fluid pressure       (fp)
   2: structure            (s)
   3: interface multiplier (l)

 If a geometry block is present, block 4 stores the geometry unknowns (g).
 If an SCI block is present without geometry, block 4 stores the chemical
 unknowns. If geometry and SCI are both present, block 5 stores the chemical
 unknowns. The algorithm is applied block-wise: structure/SCI solve, optional
 geometry correction, interface correction, fluid condensation and fluid solve,
 followed by the final coupling residual update.
 */

namespace FEDD {
using namespace Thyra;
        
// Constructors

template<class SC, class LO, class GO, class NO>
PrecOpFaCSI<SC,LO,GO,NO>::PrecOpFaCSI()
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(false),
useSolidPreconditioner_(true)
{
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}

template<class SC, class LO, class GO, class NO>
PrecOpFaCSI<SC,LO,GO,NO>::PrecOpFaCSI(CommConstPtr_Type comm, bool fluidPrecMonolithic, bool useFluidPreconditioner, bool useSolidPreconditioner, bool onlyDiagonal)
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(fluidPrecMonolithic),
useFluidPreconditioner_(useFluidPreconditioner),
useSolidPreconditioner_(useSolidPreconditioner),
onlyDiagonal_(onlyDiagonal)
{
    comm_=comm;
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGE(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type sInv,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT){

    // GE setup: FaCSI without a geometry inverse. The structure block is
    // handled by sInv_, the fluid velocity-pressure block by fInv_, and
    // C1_/C1T_/C2_ provide the interface coupling used in applyImpl().
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);

    initialize();
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGE(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type sciInv,
                                     ThyraLinOpPtr_Type sciS,
                                     ThyraLinOpPtr_Type sciC,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT){

    // GE setup with SCI: use the structure-chemical inverse instead of a
    // standalone structure inverse.
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setSCIC(sciC);
    setSCIS(sciS);
    setSCIInv(sciInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);

    initializeWithSCI();
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGI(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type C4,
                                     ThyraLinOpPtr_Type sInv,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT,
                                     ThyraLinOpPtr_Type gInv){
    
    // GI setup: extend GE by adding the geometry coupling C4_ and geometry
    // inverse gInv_.
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setC4(C4);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);
    setGeoInv(gInv);

    initialize();
}


template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGI(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type C4,
                                     ThyraLinOpPtr_Type sciInv,
                                     ThyraLinOpPtr_Type sciS,
                                     ThyraLinOpPtr_Type sciC,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT,
                                     ThyraLinOpPtr_Type gInv){
    
    // GI setup with SCI: geometry and structure-chemical blocks are both
    // active, so applyImpl() expects the chemical block after the geometry
    // block in the product vector.
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setC4(C4);
    setSCIC(sciC);
    setSCIS(sciS);
    setSCIInv(sciInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);
    setGeoInv(gInv);

    initializeWithSCI();
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGIShape(ThyraLinOpPtr_Type C1,
                                          ThyraLinOpPtr_Type C1T,
                                          ThyraLinOpPtr_Type C2,
                                          ThyraLinOpPtr_Type C4,
                                          ThyraLinOpPtr_Type sInv,
                                          ThyraLinOpPtr_Type fInv,
                                          ThyraLinOpPtr_Type fF,
                                          ThyraLinOpPtr_Type fBT,
                                          ThyraLinOpPtr_Type gInv,
                                          ThyraLinOpPtr_Type shape_v,
                                          ThyraLinOpPtr_Type shape_p){
    
    // Shape-derivative variant of GI: after the geometry solve, shape_v_ and
    // shape_p_ correct the fluid velocity and pressure residuals.
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setC4(C4);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);
    setGeoInv(gInv);
    setShapeDeriv(shape_v, shape_p);
 
    initialize();
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setCE(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type sciInv,
                                     ThyraLinOpPtr_Type sciS,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT){

    // CE setup: SCI preconditioning without an explicit chemical coupling
    // block in applyImpl().
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setSCIS(sciS);
    setSCIInv(sciInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);

    initializeWithSCI();
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC1(ThyraLinOpPtr_Type C1){
    C1_ = C1;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC1T(ThyraLinOpPtr_Type C1T){
    C1T_ = C1T;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC2(ThyraLinOpPtr_Type C2){
    C2_ = C2;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC4(ThyraLinOpPtr_Type C4){
    C4_ = C4;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setShapeDeriv(ThyraLinOpPtr_Type shape_v, ThyraLinOpPtr_Type shape_p){
    shape_v_ = shape_v;
    shape_p_ = shape_p;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setStructInv(ThyraLinOpPtr_Type sInv){
    sInv_ = sInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidInv(ThyraLinOpPtr_Type fInv){
    fInv_ = fInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGeoInv(ThyraLinOpPtr_Type gInv){
    gInv_ = gInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidF(ThyraLinOpPtr_Type fF){
    fF_ = fF;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidBT(ThyraLinOpPtr_Type fBT){
    fBT_ = fBT;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setSCIInv(ThyraLinOpPtr_Type sciInv){
    sciInv_ = sciInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setSCIC(ThyraLinOpPtr_Type sciC){
    sciC_ = sciC;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setSCIS(ThyraLinOpPtr_Type sciS){
    sciS_ = sciS;
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::initialize(){
    TEUCHOS_TEST_FOR_EXCEPTION(fInv_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: Fluid preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(sInv_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: Structure preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(C1_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: C1 not set.");

    // Product vector spaces mirror the block layout used in applyImpl():
    // fluid velocity, fluid pressure, structure and interface multiplier,
    // followed by an optional geometry block.
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( 4 );
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( 4 );
    vectorSpacesRange[0] = fF_->range();
    vectorSpacesRange[1] = fBT_->domain();
    vectorSpacesRange[2] = sInv_->range();
    vectorSpacesRange[3] = C1_->range();
    
    vectorSpacesDomain[0] = fF_->domain();
    vectorSpacesDomain[1] = fBT_->domain();
    vectorSpacesDomain[2] = sInv_->domain();
    vectorSpacesDomain[3] = C1T_->domain();
    if ( !gInv_.is_null() ) {
        vectorSpacesRange.push_back( gInv_->range() );
        vectorSpacesDomain.push_back( gInv_->domain() );
    }
    //     defaultProductRange_;
    //    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductDomain_;
    
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSR = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pR);
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSD = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pD);
//    
//    this->defaultProductRange_ = Thyra::multiVectorProductVectorSpace<SC>( pVSR , 1);
//    this->defaultProductDomain_ = Thyra::multiVectorProductVectorSpace<SC>( pVSD, 1);
    
    this->defaultProductRange_ = pR;
    this->defaultProductDomain_ = pD;

}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::initializeWithSCI(){


    TEUCHOS_TEST_FOR_EXCEPTION(fInv_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: Fluid preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(sciInv_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: Structure preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(C1_.is_null(), std::runtime_error,"Can not initialize FaCSCI preconditioner: C1 not set.");

    // SCI uses the same first four FaCSI blocks, but the structure space is
    // taken from sciS_. Optional geometry and chemical blocks are appended in
    // the same order expected by applyImpl().
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( 4 );
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( 4 );
    vectorSpacesRange[0] = fF_->range();
    vectorSpacesRange[1] = fBT_->domain();
    vectorSpacesRange[2] = sciS_->range();
    vectorSpacesRange[3] = C1_->range();
    
    vectorSpacesDomain[0] = fF_->domain();
    vectorSpacesDomain[1] = fBT_->domain();
    vectorSpacesDomain[2] = sciS_->domain();
    vectorSpacesDomain[3] = C1T_->domain();

    if ( !gInv_.is_null() ) {
        vectorSpacesRange.push_back( gInv_->range() );
        vectorSpacesDomain.push_back( gInv_->domain() );
    }
    if(!sciC_.is_null()){
        vectorSpacesRange.push_back( sciC_->range() );
        vectorSpacesDomain.push_back( sciC_->domain() );     
    }
    //     defaultProductRange_;
    //    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductDomain_;
    
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSR = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pR);
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSD = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pD);
//    
//    this->defaultProductRange_ = Thyra::multiVectorProductVectorSpace<SC>( pVSR , 1);
//    this->defaultProductDomain_ = Thyra::multiVectorProductVectorSpace<SC>( pVSD, 1);
    
    this->defaultProductRange_ = pR;
    this->defaultProductDomain_ = pD;

}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::applyIt(
                                         const EOpTransp M_trans,
                                         const MultiVectorBase<SC> &X_in,
                                         const Ptr<MultiVectorBase<SC> > &Y_inout,
                                         const SC alpha,
                                         const SC beta
                                         ) const
{
    applyImpl(M_trans, X_in, Y_inout, alpha, beta);

}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::applyImpl(
                                   const EOpTransp M_trans,
                                   const MultiVectorBase<SC> &X_in,
                                   const Ptr<MultiVectorBase<SC> > &Y_inout,
                                   const SC alpha,
                                   const SC beta
                                   ) const
{
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    using Teuchos::rcpFromRef;
    typedef Teuchos::ScalarTraits<SC> ST;
    typedef RCP<MultiVectorBase<SC> > MultiVectorPtr;
    typedef RCP<const MultiVectorBase<SC> > ConstMultiVectorPtr;
    typedef RCP<const LinearOpBase<SC> > ConstLinearOpPtr;

    int rank = comm_->getRank();

    // Traditional notation used below:
    //
    //   X_fv, X_fp, X_s, X_l, X_g
    //       = r_u_f, r_p, r_d_s, r_lambda, r_d_f
    //
    //   Y_fv, Y_fp, Y_s, Y_l, Y_g
    //       = u_f, p, d_s, lambda, d_f
    //
    // Here u_f, p, d_s, lambda and d_f denote the preconditioned corrections
    // produced by this operator. Hat variables denote intermediate corrected
    // residuals used before the corresponding block solve.

    // Diagnostic mode: skip the FaCSI coupling steps and copy the four primary
    // blocks. This is useful for comparing the fully coupled preconditioner
    // against the unpreconditioned block residual.
    //
    //   Y_fv = X_fv,  Y_fp = X_fp,  Y_s = X_s,  Y_l = X_l
    //   u_f = r_u_f,  p = r_p,  d_s = r_d_s,  lambda = r_lambda
    if (onlyDiagonal_) {
        
        Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
        = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( rcpFromRef(X_in) );
        
        Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
        = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );
        
        Y_inout->assign(0.);

        Teuchos::RCP< const MultiVectorBase< SC > > X_s = X->getMultiVectorBlock(2);
        Teuchos::RCP< MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);

        
        Teuchos::RCP< const MultiVectorBase< SC > > X_fv = X->getMultiVectorBlock(0);
        Teuchos::RCP< MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
        Teuchos::RCP< const MultiVectorBase< SC > > X_fp = X->getMultiVectorBlock(1);
        Teuchos::RCP< MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
        
        
        Teuchos::RCP< MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);
        Teuchos::RCP<const MultiVectorBase< SC > > X_l = X->getMultiVectorBlock(3);

        assign(Y_fv.ptr(), *X_fv);
        assign(Y_fp.ptr(), *X_fp);
        assign(Y_s.ptr(), *X_s);
        assign(Y_l.ptr(), *X_l);
                
    }
    else{
        Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
            = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( rcpFromRef(X_in) );

        Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
            = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );

        
        Y_inout->assign(0.);

        // Step 1: structure block.
        // Start from the structure right-hand side and apply either the
        // standalone structure inverse sInv_ or the SCI inverse. In the SCI
        // case the structure and chemical blocks are packed into the
        // monolithic SCI work vectors before applying sciInv_.
        //
        //   standard: Y_s = sInv_ * X_s
        //   SCI:      [Y_s, Y_chem]^T = sciInv_ * [X_s, X_chem]^T
        //
        //   standard: d_s = S^{-1} r_d_s
        //   SCI:      [d_s, c]^T = SCI^{-1} [r_d_s, r_c]^T
        //
        // If structure preconditioning is disabled, the code keeps
        //
        //   Y_s = X_s  and, for SCI,  Y_chem = X_chem.
        //   d_s = r_d_s  and, for SCI,  c = r_c.
        Teuchos::RCP< const MultiVectorBase< SC > > X_s = X->getMultiVectorBlock(2);
        Teuchos::RCP< MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);
        assign(Y_s.ptr(), *X_s);  
    
        // Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XsTpetra =
        // Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_s );
    
        
        if (useSolidPreconditioner_){
            if(!sInv_.is_null())
                sInv_->apply(NOTRANS, *X_s, Y_s.ptr(), 1., 0.);
            else{
                if(!sciC_.is_null()){
                    // std::cout << "FACSCI:: Implicit Case " << std::endl;
                    Teuchos::RCP< const MultiVectorBase< SC > > X_chem;
                    Teuchos::RCP< MultiVectorBase< SC > > Y_chem ;
            
                    if ( !gInv_.is_null() ) {
                        X_chem = X->getMultiVectorBlock(5);
                        Y_chem = Y->getNonconstMultiVectorBlock(5);
                    }
                    else{
                        X_chem = X->getMultiVectorBlock(4);
                        Y_chem = Y->getNonconstMultiVectorBlock(4);
                    }
                    
                    assign(Y_chem.ptr(), *X_chem);
                    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_sci( 2 );
                    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_sci( 2 );

                    X_sci[0] = Y_s;
                    X_sci[1] = Y_chem;

                    Y_sci[0] = Y_s;
                    Y_sci[1] = Y_chem;
                    Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > sciMonoVS = sciInv_->domain();
                    // std::cout << " Initi X_scimono_ " << std::endl;
                    if ( X_scimono_.is_null() ){
                        X_scimono_ = createMembers( sciMonoVS, X_sci[0]->domain()->dim() );
                        Y_scimono_ = createMembers( sciMonoVS, Y_sci[0]->domain()->dim() );
                    }
                    //We can/should speedup this process
                   // std::cout << " Copy to mono " << std::endl;

                    copyToMonoSCI(X_sci);
                    //std::cout << " Apply " << std::endl;

                    sciInv_->apply(NOTRANS, *X_scimono_, Y_scimono_.ptr(), 1., 0.);
                    copyFromMonoSCI(Y_sci); 
                }
                else{
                    sciInv_->apply(NOTRANS, *X_s, Y_s.ptr(), 1., 0.);        
                }  
            }

        }
        else{
            assign(Y_s.ptr(), *X_s);
            if(!sciC_.is_null()){
                Teuchos::RCP< const MultiVectorBase< SC > > X_chem = X->getMultiVectorBlock(4);
                Teuchos::RCP< MultiVectorBase< SC > > Y_chem = Y->getNonconstMultiVectorBlock(4);
                assign(Y_chem.ptr(), *X_chem);
            }

        }

        Teuchos::RCP< const MultiVectorBase< SC > > X_g;
        Teuchos::RCP< MultiVectorBase< SC > > Y_g;
        // Step 2: optional geometry block.
        // Remove the contribution of the already preconditioned structure
        // block through C4_, then apply the geometry inverse.
        //
        //   Y_g = gInv_ * (X_g - C4_ * Y_s)
        //   d_f = G^{-1} (r_d_f - C4 d_s)
        if ( !gInv_.is_null() ) {
            X_g = X->getMultiVectorBlock(4);
            Y_g = Y->getNonconstMultiVectorBlock(4);

            assign(Y_g.ptr(), *X_g);

            C4_->apply(NOTRANS, *Y_s, Y_g.ptr(), -1., 1.);
            
            gInv_->apply(NOTRANS, *Y_g, Y_g.ptr(), 1., 0.);
        }
        
        // Step 3: interface multiplier block.
        // Correct the multiplier right-hand side with the structure coupling
        // C2_ * Y_s before the fluid condensation step uses Y_l.
        //
        //   Y_l = X_l - C2_ * Y_s
        //   r_hat_lambda = r_lambda - C2 d_s
        Teuchos::RCP< const MultiVectorBase< SC > > X_l = X->getMultiVectorBlock(3);
        Teuchos::RCP< MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);

        assign(Y_l.ptr(), *X_l);
        
        C2_->apply( NOTRANS, *Y_s, Y_l.ptr(), -1., 1. );
          
        Teuchos::RCP< const MultiVectorBase< SC > > X_fv = X->getMultiVectorBlock(0);
        Teuchos::RCP< MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
        assign(Y_fv.ptr(), *X_fv);

        Teuchos::RCP< const MultiVectorBase< SC > > X_fp = X->getMultiVectorBlock(1);
        Teuchos::RCP< MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
        assign(Y_fp.ptr(), *X_fp);

        // If shape derivatives are available, subtract the geometry
        // contribution from the fluid velocity and pressure residuals.
        //
        //   Y_fv = X_fv - shape_v_ * Y_g
        //   Y_fp = X_fp - shape_p_ * Y_g
        //   r_hat_u_f = r_u_f - shape_v d_f
        //   r_hat_p   = r_p   - shape_p d_f
        //
        // Without shape derivatives, this remains
        //
        //   Y_fv = X_fv,  Y_fp = X_fp.
        //   r_hat_u_f = r_u_f,  r_hat_p = r_p.
        if (!shape_v_.is_null() && !shape_p_.is_null()) {
            shape_v_->apply(NOTRANS, *Y_g, Y_fv.ptr(), -1., 1.);
            shape_p_->apply(NOTRANS, *Y_g, Y_fp.ptr(), -1., 1.);
        }
        
        // Keep a copy of the original velocity residual. Step 6 uses it to
        // rebuild the momentum residual after the fluid block has been solved.
        //
        //   Z_fv_ = Y_fv
        //   z_u_f = r_hat_u_f
        if (Z_fv_.is_null())
            Z_fv_ = Y_fv->clone_mv();
        else
            assign(Z_fv_.ptr(), *Y_fv);
        
        // Step 4: fluid condensation.
        // Eliminate the interface contribution from the fluid velocity
        // equation:
        //
        //   Y_fv <- Y_fv - C1T_ * C1_ * Y_fv + C1T_ * Y_l
        //   r_tilde_u_f = r_hat_u_f - C1^T C1 r_hat_u_f
        //                 + C1^T r_hat_lambda
        //              = [r_tilde_u_f_I,r_lambda_Gamma]^T
        //
        // The modified (Y_fv, Y_fp) pair is the right-hand side for the fluid
        // preconditioner in the next step.
        // Replace Interface values of u with lamba values bc.: x_\Gamma = r_\lambda
        if (tmp_l_.is_null())
            tmp_l_ = Y_l->clone_mv();

        C1_->apply(NOTRANS, *Y_fv, tmp_l_.ptr(), 1., 0);
        C1T_->apply(NOTRANS, *tmp_l_, Y_fv.ptr(), -1., 1.);
        
        C1T_->apply(NOTRANS, *Y_l, Y_fv.ptr(), 1., 1.);
        
    //    std::cout << "W_fv" << std::endl;
    //    Y_fv->describe(*out,Teuchos::VERB_EXTREME);
    //    comm_->barrier();    comm_->barrier();    comm_->barrier();
        
        if (productRangeFluid_.is_null()) {
            Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRangeFluid( 2 );
            
            vectorSpacesRangeFluid[0] = X_fv->range();
            vectorSpacesRangeFluid[1] = X_fp->range();
            
            productRangeFluid_ = Thyra::productVectorSpace<SC>( vectorSpacesRangeFluid );
        }
        
        Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid( 2 );
        Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid( 2 );

        X_fluid[0] = Y_fv;
        X_fluid[1] = Y_fp;

        Y_fluid[0] = Y_fv;
        Y_fluid[1] = Y_fp;

        
    //    std::cout << "Y_fp before mono" << std::endl;
    //    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
    //    comm_->barrier();    comm_->barrier();    comm_->barrier();
        
        // Step 5: fluid block.
        // Apply the velocity-pressure fluid preconditioner either in the
        // monolithic storage used by fInv_ or as a Thyra product vector.
        //
        //   [Y_fv, Y_fp]^T = fInv_ * [Y_fv, Y_fp]^T
        //   [u_f, p]^T = F_f^{-1} [r_tilde_u_f, r_hat_p]^T
        //
        // If fluid preconditioning is disabled, the code restores
        //
        //   Y_fv = X_fv,  Y_fp = X_fp.
        //   u_f = r_u_f,  p = r_p.
        if (useFluidPreconditioner_){
        
            if (fluidPrecMonolithic_) {
                Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > fMonoVS = fInv_->domain();

                if ( X_fmono_.is_null() ){
                    X_fmono_ = createMembers( fMonoVS, X_fluid[0]->domain()->dim() );
                    Y_fmono_ = createMembers( fMonoVS, Y_fluid[0]->domain()->dim() );
                }
                //We can/should speedup this process
                copyToMono(X_fluid);
                fInv_->apply(NOTRANS, *X_fmono_, Y_fmono_.ptr(), 1., 0.);
                copyFromMono(Y_fluid);
            }
            else{
                Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodX_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, X_fluid );
                Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodY_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, Y_fluid );
                
                Teuchos::RCP< MultiVectorBase<SC> > X_f = Teuchos::rcp_dynamic_cast<MultiVectorBase< SC > >(prodX_f);
                Teuchos::RCP< MultiVectorBase<SC> > Y_f = Teuchos::rcp_dynamic_cast<MultiVectorBase< SC > >(prodY_f);
                
                fInv_->apply(NOTRANS, *X_f, Y_f.ptr(), 1., 0.);
            }
        }
        else{
            assign(Y_fv.ptr(), *X_fv);
            assign(Y_fp.ptr(), *X_fp);
        }

    
        // Step 6: final multiplier update.
        // Reconstruct the fluid momentum residual left after the fluid solve,
        // then map it back to the multiplier block with C1_.
        //
        //   Z_fv_ = Z_fv_ - fBT_ * Y_fp - fF_ * Y_fv
        //   Y_l   = C1_ * Z_fv_
        //
        //   z_u_f  = r_hat_u_f - B^T p - F u_f
        //   lambda = C1 z_u_f
        fBT_->apply(NOTRANS, *Y_fp, Z_fv_.ptr(), -1., 1.);
       
        fF_->apply(NOTRANS, *Y_fv, Z_fv_.ptr(), -1., 1.);

          
        C1_->apply(NOTRANS, *Z_fv_, Y_l.ptr(), 1., 0.);
           
            
    }
}
// private
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyToMono( Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid ) const{

    // Pack the product fluid vector (velocity, pressure) into the monolithic
    // layout expected by monolithic fluid preconditioners.
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fmono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );

    for (int j=0; j<X_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_fluid)(i,j) = (*thyData_v)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_fluid)( i + localSubDim_v, j ) = (*thyData_p)(i,j);

    }
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyFromMono(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid) const{

    // Scatter the monolithic fluid preconditioner result back into velocity
    // and pressure product-vector blocks.

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[0] ,Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[1] ,Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fmono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );
    
    for (int j=0; j<Y_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_v)(i,j) = (*thyData_fluid)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_p)(i,j) = (*thyData_fluid)( i + localSubDim_v, j );
        
    }    
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyToMonoSCI( Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_sci ) const{

    // Pack the SCI product vector (structure, chemical) into the monolithic
    // layout expected by sciInv_.
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_d = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_sci[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_c = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_sci[1]->range());
    
    const LO localOffset_d = mpiVS_d->localOffset();
    const LO localSubDim_d = mpiVS_d->localSubDim();
    
    const LO localOffset_c = mpiVS_c->localOffset();
    const LO localSubDim_c = mpiVS_c->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_d =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_sci[0] ,Range1D(localOffset_d,localOffset_d+localSubDim_d-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_c =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_sci[1] ,Range1D(localOffset_c,localOffset_c+localSubDim_c-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_scimono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_sci =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_scimono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );

    for (int j=0; j<X_sci[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_d; i++)
            (*thyData_sci)(i,j) = (*thyData_d)(i,j);
        
        for (LO i=0; i<localSubDim_c; i++)
            (*thyData_sci)( i + localSubDim_d, j ) = (*thyData_c)(i,j);

    }
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyFromMonoSCI(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_sci) const{

    // Scatter the SCI preconditioner result back into the structure and
    // chemical product-vector blocks.

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_d = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_sci[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_c = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_sci[1]->range());
    
    const LO localOffset_d = mpiVS_d->localOffset();
    const LO localSubDim_d = mpiVS_d->localSubDim();
    
    const LO localOffset_c = mpiVS_c->localOffset();
    const LO localSubDim_c = mpiVS_c->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_d =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_sci[0] ,Range1D(localOffset_d,localOffset_d+localSubDim_d-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_c =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_sci[1] ,Range1D(localOffset_c,localOffset_c+localSubDim_c-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_scimono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_sci =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_scimono_ ,Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );
    
    for (int j=0; j<Y_sci[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_d; i++)
            (*thyData_d)(i,j) = (*thyData_sci)(i,j);
        
        for (LO i=0; i<localSubDim_c; i++)
            (*thyData_c)(i,j) = (*thyData_sci)( i + localSubDim_d, j );
        
    }    
}


}

#endif
