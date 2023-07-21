#include "AssembleFEStokes_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFEStokes_def.hpp"
namespace FEDD {
    template class AssembleFEStokes<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

