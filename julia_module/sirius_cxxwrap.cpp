#include <sirius.hpp>
#include <jlcxx/jlcxx.hpp>
#include <mpi.h>
#include "core/json.hpp"


using json = nlohmann::json;
using namespace sirius;

//helpers to wrap the templated r3::vector class
struct WrapR3Vector
{
  template<typename TypeWrapperT>
  void operator()(TypeWrapperT&& wrapped)
  {
    typedef typename TypeWrapperT::type WrappedT;
    wrapped.template constructor<typename WrappedT::value_type, typename WrappedT::value_type, 
                                 typename WrappedT::value_type>();
    wrapped.method("length", &WrappedT::length);
    wrapped.method("length2", &WrappedT::length2);
  }
};

//2D mdarrays/matrices
struct WrapMdMatrix
{
  template<typename TypeWrapperT>
  void operator()(TypeWrapperT&& wrapped)
  {
    typedef typename TypeWrapperT::type WrappedT;
    wrapped.template constructor<>();
    //TODO: would be nice to define that function elsewhere for more readability
    //TODO: would be nice to return a full 2D array instead of element by element
    wrapped.module().method("get_element", [](const WrappedT& mat, size_t const idx) {return mat[idx];});
  }
};

namespace jlcxx
{
  template<typename T>
  struct BuildParameterList<matrix<T>>
  {
     typedef ParameterList<T> type;
  };

  template<typename T> struct IsMirroredType<matrix<T>> : std::false_type { };
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
   mod.method("initialize", &initialize);
   mod.method("finalize", &finalize);

   //Simulation context class
   mod.add_type<Simulation_context>("SimulationContext")
      .constructor<std::string const&>() //, mpi::Communicator const&>(); use COMM_WORLD for now
      .method("initialize", &Simulation_context::initialize);

   //r3::vector class
   mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("R3Vector") //assign Julia type of AbstractArray
      .apply<r3::vector<int>, r3::vector<double>>(WrapR3Vector()); //both int and double available
                                                                   
   //mdarray class
   mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("MdMatrix")
      .apply<matrix<double>, matrix<int>>(WrapMdMatrix());

   //Kpoint set class
   mod.add_type<K_point_set>("KPointSet")
      .constructor<Simulation_context& , r3::vector<int>, r3::vector<int>, int>();

   //json class, necessary for gs.find()
   //Note: the json class is actually quite opaque as it comes from an external library.
   //      I would rather querry the DFT_groud_state object for items of interest
   mod.add_type<json>("json");

   //Force class
   mod.add_type<Force>("Force")
      .method("forces_total", &Force::forces_total)
      .method("calc_forces_total", &Force::calc_forces_total);

   //Ground state class
   mod.add_type<DFT_ground_state>("GroundState")
      .constructor<K_point_set&>()
      .method("initial_state", &DFT_ground_state::initial_state)
      .method("total_energy", &DFT_ground_state::total_energy)
      .method("forces", &DFT_ground_state::forces)
      .method("find", &DFT_ground_state::find);
}
