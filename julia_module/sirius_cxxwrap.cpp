#include <sirius.hpp>
#include <jlcxx/jlcxx.hpp>
#include <mpi.h>
#include "core/json.hpp"
#include "core/mpi/communicator.hpp"


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
    wrapped.module().method("get_element", [](const WrappedT& vec, const int i) {return vec[i];});
  }
};

//helpers to wrap the templated r3::matrix class
struct WrapR3Matrix
{
  template<typename TypeWrapperT>
  void operator()(TypeWrapperT&& wrapped)
  {
    typedef typename TypeWrapperT::type WrappedT;
    wrapped.template constructor<>();
    wrapped.module().method("get_element", [](const WrappedT& mat, const int i, const int j) {return mat(i, j);});
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

//MPI communicator wrapper
class CommWrap
{
   private:
      mpi::Communicator mpi_comm_;

   public:
      CommWrap(int fcomm)
      {
         mpi_comm_ = mpi::Communicator::map_fcomm(fcomm);
      }

      inline mpi::Communicator& 
      get_comm()
      {
         return mpi_comm_;
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
   //r3::vector class
   mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("R3Vector")
      .apply<r3::vector<int>, r3::vector<double>>(WrapR3Vector()); //both int and double available
                                                                   
   //r3::matrix class
   mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("R3Matrix")
      .apply<r3::matrix<int>, r3::matrix<double>>(WrapR3Matrix());
                                                                   
   //mdarray class
   mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("MdMatrix")
      .apply<matrix<double>, matrix<int>>(WrapMdMatrix());

   //MPI comm
   mod.add_type<mpi::Communicator>("MPICommunicator");

   //MPI comm wrapped
   mod.add_type<CommWrap>("CommWrap")
      .constructor<int>()
      .method("get_comm", &CommWrap::get_comm);

   //General init/finalize
   mod.method("initialize", &initialize);
   mod.method("finalize", &finalize);

   //Simulation context class
   mod.add_type<Simulation_context>("SimulationContext")
      .constructor<std::string const&>() //, mpi::Communicator const&>(); use COMM_WORLD for now
      .constructor<std::string const&, mpi::Communicator const&>()
      .method("initialize", &Simulation_context::initialize);

   //useful getters from the Simulation_context
   mod.method("get_density_tol", [](const Simulation_context& ctx) 
              {return ctx.cfg().parameters().density_tol();});
   mod.method("get_energy_tol", [](const Simulation_context& ctx) 
              {return ctx.cfg().parameters().energy_tol();});
   mod.method("get_initial_tol", [](const Simulation_context& ctx) 
              {return ctx.cfg().iterative_solver().energy_tolerance();});
   mod.method("get_num_dft_iter", [](const Simulation_context& ctx) 
              {return ctx.cfg().parameters().num_dft_iter();});
   mod.method("get_ngridk", [](const Simulation_context& ctx)
              {auto gk = ctx.cfg().parameters().ngridk();
               return r3::vector<int>(gk[0], gk[1], gk[2]);}); 
   mod.method("get_shiftk", [](const Simulation_context& ctx)
              {auto sk = ctx.cfg().parameters().shiftk();
               return r3::vector<int>(sk[0], sk[1], sk[2]);});

   //Kpoint set class
   mod.add_type<K_point_set>("KPointSet")
      .constructor<Simulation_context& , r3::vector<int>, r3::vector<int>, int>();

   //json class, necessary for gs.find()
   mod.add_type<json>("json");

   //Force class
   mod.add_type<Force>("Force")
      .method("forces_total", &Force::forces_total)
      .method("calc_forces_total", &Force::calc_forces_total);

   //Stress class
   mod.add_type<Stress>("Stress")
      .method("stress_total", &Stress::stress_total)
      .method("calc_stress_total", &Stress::calc_stress_total);

   //Ground state class
   mod.add_type<DFT_ground_state>("GroundState")
      .constructor<K_point_set&>()
      .method("initial_state", &DFT_ground_state::initial_state)
      .method("total_energy", &DFT_ground_state::total_energy)
      .method("forces", &DFT_ground_state::forces)
      .method("stress", &DFT_ground_state::stress)
      .method("find", &DFT_ground_state::find);
}
