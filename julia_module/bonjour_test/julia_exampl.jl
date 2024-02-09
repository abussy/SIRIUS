# Load the module and generate the functions
module CppBonjour
  using CxxWrap
  @wrapmodule(() -> joinpath("/home/bussya/Documents/git/spack/opt/spack/linux-ubuntu23.10-skylake/gcc-13.2.0/sirius-develop-doszdkbjuev3izwyetkiflbdnn3iog4m/lib/julia_wrapper","libbonjour"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppBonjour.greet()
