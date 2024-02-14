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
   mod.method("sirius_initialize", &sirius_initialize);
   mod.method("sirius_finalize", &sirius_finalize);
   mod.method("sirius_create_context", &sirius_create_context);
   mod.method("sirius_free_object_handler", &sirius_free_object_handler);
   mod.method("sirius_context_initialized", &sirius_context_initialized);
   mod.method("sirius_initialize_context", &sirius_initialize_context);
}

// So it seems we can direclty import functions from the C-style API
// However, I dislike the idea of importing a cpp file. We probably need
// to provide headers for the whole API
//
// The most important question for now is: can we pass a MPI comm from julia to sirius ?!?
