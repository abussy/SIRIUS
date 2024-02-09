# Load the module and generate the functions
module SiriusJl
  using CxxWrap
  @wrapmodule(() -> joinpath("/home/bussya/Documents/git/spack/opt/spack/linux-ubuntu23.10-skylake/gcc-13.2.0/sirius-develop-t5xeyb3ko6za5zyyjyd44ht4e36mgzxk/lib/julia_wrapper","libsirius_jl"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show SiriusJl.sirius_test()
@show SiriusJl.greet()
@show SiriusJl.from_api()
