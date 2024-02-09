#include "api/test.hpp"
#include "jlcxx/jlcxx.hpp"
#include "api/sirius_api.cpp"


std::string greet()
{
   return "blablabla";
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
   mod.method("sirius_test", &sirius_test);
   mod.method("greet", &greet);
   mod.method("from_api", &from_api);
}

// So it seems we can direclty import functions from the C-style API
// However, I dislike the idea of importing a cpp file. We probably need
// to provide headers for the whole API
//
// The most important question for now is: can we pass a MPI comm from julia to sirius ?!?
