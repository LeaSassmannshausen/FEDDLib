#ifndef NonLinLaplace_def_hpp
#define NonLinLaplace_def_hpp
#include "NonLinLaplace_decl.hpp"
/*!
 Definition of NonLinLaplace

 @brief NonLinLaplace
 @author Kyrill Ho
 @version 1.0
 @copyright KH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
NonLinLaplace<SC, LO, GO, NO>::NonLinLaplace(
    const DomainConstPtr_Type &domain, std::string FEType,
    ParameterListPtr_Type parameterList)
    : NonLinearProblem<SC, LO, GO, NO>(parameterList, domain->getComm()),
      u_rep_() {
    this->nonLinearTolerance_ =
        this->parameterList_->sublist("Parameter").get("relNonLinTol", 1.0e-6);
    this->initNOXParameters();
    this->addVariable(domain, FEType, "u", 1);
    this->dim_ = this->getDomain(0)->getDimension();
    u_rep_ = Teuchos::rcp(
        new MultiVector_Type(this->getDomain(0)->getMapRepeated()));
}

template <class SC, class LO, class GO, class NO>
NonLinLaplace<SC, LO, GO, NO>::~NonLinLaplace() {}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::info() {
    this->infoProblem();
    this->infoNonlinProblem();
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::assemble(std::string type) const {
    if (type == "") {
        if (this->verbose_) {
            std::cout << "-- Assembly nonlinear laplace ... " << std::flush;
        }
        this->initAssemble();
    } else {
        this->reAssemble(type);
    }
    if (this->verbose_) {
        std::cout << "done -- " << std::endl;
    }
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::initAssemble() const {

    MatrixPtr_Type A(
        new Matrix_Type(this->getDomain(0)->getMapUnique(),
                        this->getDomain(0)->getApproxEntriesPerRow()));
    if (this->verbose_) {
        std::cout << "-- Initial assembly " << std::flush;
    }

    if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(1));

    if (this->residualVec_.is_null())
        this->residualVec_.reset(new BlockMultiVector_Type(1));

    this->system_->addBlock(A, 0, 0);
    // TODO here the initialisation of the system may still be wrong.
    /* this->assembleSourceTerm(0.); */
    /* this->addToRhs(this->sourceTerm_); */
    this->feFactory_->assemblyNonlinearLaplace(
        this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
        this->system_, this->residualVec_, this->parameterList_, "Jacobian");
    /* this->feFactory_->assemblyNonlinearLaplace( */
    /*     this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_, */
    /*     this->system_, this->residualVec_, this->parameterList_, "Rhs"); */

    // Initialise solution to 1 everywhere
    this->solution_->putScalar(0.);

    if (this->verbose_) {
        std::cout << "done -- " << std::endl;
    }
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::reAssemble(std::string type) const {

    if (this->verbose_)
        std::cout << "-- Reassembly nonlinear laplace"
                  << " (" << type << ") ... " << std::flush;

    MatrixPtr_Type W = Teuchos::rcp(
        new Matrix_Type(this->getDomain(0)->getMapUnique(),
                        this->getDomain(0)->getApproxEntriesPerRow()));

    MultiVectorConstPtr_Type u = this->solution_->getBlock(0);

    this->u_rep_->importFromVector(u, true);

    this->system_->addBlock(W, 0, 0);

    // MultiVectorPtr_Type f = Teuchos::rcp(
    //     new MultiVector_Type(this->getDomain(0)->getMapRepeated(),
    //     1));

    if (type == "Rhs") {

        this->feFactory_->assemblyNonlinearLaplace(
            this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
            this->system_, this->residualVec_, this->parameterList_, "Rhs");

    } else if (type == "Newton") {

        this->feFactory_->assemblyNonlinearLaplace(
            this->dim_, this->getDomain(0)->getFEType(), 2, this->u_rep_,
            this->system_, this->residualVec_, this->parameterList_,
            "Jacobian");
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::reAssembleExtrapolation(
    BlockMultiVectorPtrArray_Type previousSolutions) {

    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Only Newton/NOX implemented for nonlinear material models!");
}

// TODO setup for nox solver?
template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    TEUCHOS_TEST_FOR_EXCEPTION(
        this->solution_->getBlock(0)->getMap()->getUnderlyingLib() != "Tpetra",
        std::runtime_error,
        "Use of NOX only supports Tpetra. Epetra support must be implemented.");
    RCP<Teuchos::FancyOStream> fancy =
        Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    TEUCHOS_TEST_FOR_EXCEPTION(inArgs.get_x().is_null(), std::logic_error,
                               "inArgs.get_x() is null.");

    RCP<const Thyra::VectorBase<SC>> vecThyra = inArgs.get_x();
    RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();

    RCP<Thyra::VectorBase<SC>> vecThyraNonConst =
        rcp_const_cast<Thyra::VectorBase<SC>>(vecThyra);

    this->solution_->fromThyraMultiVector(vecThyraNonConst);

    const RCP<Thyra::MultiVectorBase<SC>> f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC>> W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC>> W_prec_out = outArgs.get_W_prec();

    typedef Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, NO>
        tpetra_extract;
    typedef Xpetra::Matrix<SC, LO, GO, NO> XpetraMatrix_Type;
    typedef RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;
    typedef RCP<const XpetraMatrix_Type> XpetraMatrixConstPtr_Type;

    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);

    if (fill_f || fill_W || fill_W_prec) {

        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {

            this->calculateNonLinResidualVec("standard");

            RCP<Thyra::MultiVectorBase<SC>> f_thyra =
                this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }

        XpetraMatrixPtr_Type W;
        if (fill_W) {

            this->reAssemble("Newton");

            this->setBoundariesSystem();

            RCP<TpetraOp_Type> W_tpetra =
                tpetra_extract::getTpetraOperator(W_out);
            RCP<TpetraMatrix_Type> W_tpetraMat =
                rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

            XpetraMatrixConstPtr_Type W_systemXpetra =
                this->getSystem()->getBlock(0, 0)->getXpetraMatrix();

            XpetraMatrixPtr_Type W_systemXpetraNonConst =
                rcp_const_cast<XpetraMatrix_Type>(W_systemXpetra);
            Xpetra::CrsMatrixWrap<SC, LO, GO, NO> &crsOp =
                dynamic_cast<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> &>(
                    *W_systemXpetraNonConst);
            Xpetra::TpetraCrsMatrix<SC, LO, GO, NO> &xTpetraMat =
                dynamic_cast<Xpetra::TpetraCrsMatrix<SC, LO, GO, NO> &>(
                    *crsOp.getCrsMatrix());
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatXpetra =
                xTpetraMat.getTpetra_CrsMatrixNonConst();

            W_tpetraMat->resumeFill();

            for (auto i = 0;
                 i < tpetraMatXpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC, LO, GO,
                                           NO>::local_inds_host_view_type
                    indices; // ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC, LO, GO,
                                           NO>::values_host_view_type values;
                tpetraMatXpetra->getLocalRowView(i, indices, values);
                W_tpetraMat->replaceLocalValues(i, indices, values);
            }
            W_tpetraMat->fillComplete();
        }

        if (fill_W_prec) {
            this->setupPreconditioner("Monolithic");

            // ch 26.04.19: After each setup of the preconditioner we check if
            // we use a two-level precondtioner with multiplicative combination
            // between the levels. If this is the case, we need to pre apply the
            // coarse level to the residual(f_out).

            std::string levelCombination =
                this->parameterList_->sublist("ThyraPreconditioner")
                    .sublist("Preconditioner Types")
                    .sublist("FROSch")
                    .get("Level Combination", "Additive");
            if (!levelCombination.compare("Multiplicative")) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                           "Multiplicative Level Combination "
                                           "is not supported for NOX. In "
                                           "general we need to pre-apply the "
                                           "coarse problem. This must be "
                                           "implemented here.");
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC>>
NonLinLaplace<SC, LO, GO, NO>::create_W_op() const {

    Teuchos::RCP<const Thyra::LinearOpBase<SC>> W_opConst =
        this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC>> W_op =
        Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC>>(W_opConst);
    return W_op;
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC>>
NonLinLaplace<SC, LO, GO, NO>::create_W_prec() const {
    this->initializeSolverBuilder();
    this->initializePreconditioner();

    Teuchos::RCP<const Thyra::PreconditionerBase<SC>> thyraPrec =
        this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC>> thyraPrecNonConst =
        Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC>>(thyraPrec);

    return thyraPrecNonConst;
}

template <class SC, class LO, class GO, class NO>
void NonLinLaplace<SC, LO, GO, NO>::calculateNonLinResidualVec(
    std::string type, double time) const {

    this->reAssemble("Rhs");
    if (!type.compare("standard")) {
        this->residualVec_->update(-1., *this->rhs_, 1.);
        /* if (!this->sourceTerm_.is_null()) */
        /*   this->residualVec_->update(-1., *this->sourceTerm_, 1.); */
    } else if (!type.compare("reverse")) {
        this->residualVec_->update(1., *this->rhs_,
                                   -1.); // this = -1*this + 1*rhs
                                         /* if (!this->sourceTerm_.is_null()) */
        /*   this->residualVec_->update(1., *this->sourceTerm_, 1.); */
    } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                                   "Unknown type for residual computation.");
    }
    // this might be set again by the TimeProblem after adding of M*u
    this->bcFactory_->setBCMinusVector(this->residualVec_, this->solution_,
                                       time);
}
} // namespace FEDD
#endif
