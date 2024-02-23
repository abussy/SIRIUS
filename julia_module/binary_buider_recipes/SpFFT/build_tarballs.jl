# Note that this script can accept some limited command-line arguments, run                          
# `julia build_tarballs.jl --help` to see a usage message.                                           
using BinaryBuilder, Pkg                                                                             
using Base.BinaryPlatforms                                                                           
const YGGDRASIL_DIR = "../.."                                                                        
include(joinpath(YGGDRASIL_DIR, "platforms", "mpi.jl"))

name = "SpFFT"
version = v"1.0.6"

sources = [
   GitSource("https://github.com/eth-cscs/SpFFT/", "aa6653f044dc8f6dbf5dc7befe45db7ce353938e")
]

script = raw"""
cd $WORKSPACE/srcdir

#sed -i 's/CMAKE_CXX_STANDARD 11/CMAKE_CXX_STANDARD 17/' ./CMakeLists.txt
#sed -i 's/CMAKE_CUDA_STANDARD 11/CMAKE_CUDA_STANDARD 17/' ./CMakeLists.txt

mkdir build
cd build

CMAKE_ARGS="-DSPFFT_OMP=ON \
            -DSPFFT_MPI=ON \
            -DSPFFT_INSTALL=ON \
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
            -DCMAKE_FIND_ROOT_PATH='${prefix}/lib/mpich;${prefix}' \
            -DCMAKE_INSTALL_PREFIX=$prefix \
            -DCMAKE_BUILD_TYPE=Release \
            -DBUILD_SHARED_LIBS=ON \
            -DMPI_C_COMPILER=$bindir/mpicc \
            -DMPI_CXX_COMPILER=$bindir/mpicxx"

#TODO: only works with MPItrampoline
if [[ "${target}" == *-apple-* ]]; then
  CMAKE_ARGS="${CMAKE_ARGS} \
               -DMPI_C_LIB_NAMES='mpi;pmpi;hwloc' \
               -DMPI_CXX_LIB_NAMES='mpicxx;mpi;pmpi;hwloc' \
               -DMPI_mpicxx_LIBRARY=${prefix}/lib/mpich/lib/libmpicxx.a \
               -DMPI_mpi_LIBRARY=${prefix}/lib/mpich/lib/libmpi.a \
               -DMPI_pmpi_LIBRARY=${prefix}/lib/mpich/lib/libpmpi.a \
               -DMPI_hwloc_LIBRARY=${prefix}/lib/libhwloc.dylib"
fi

cmake .. ${CMAKE_ARGS}

make -j${nproc} install

"""

augment_platform_block = """
    using Base.BinaryPlatforms
    $(MPI.augment)
    augment_platform!(platform::Platform) = augment_mpi!(platform)
"""
#platforms = supported_platforms()                                                       
platforms = [Platform("x86_64", "linux")]
filter!(!Sys.iswindows, platforms)
filter!(!Sys.isapple, platforms) #TODO: fix apple installation
filter!(!Sys.isfreebsd, platforms)
filter!(p -> !(libc(p) == "musl"), platforms)
platforms = expand_cxxstring_abis(platforms)

platforms = expand_gfortran_versions(platforms)
filter!(p -> !(libgfortran_version(p) < v"5"), platforms)

products = [
   LibraryProduct("libspfft", :libspfft)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("FFTW_jll"), 
    Dependency("CompilerSupportLibraries_jll"), 
    Dependency("LLVMOpenMP_jll", platforms=filter(Sys.isapple, platforms)),
]

platforms, platform_dependencies = MPI.augment_platforms(platforms; MPItrampoline_compat="5.3.0", OpenMPI_compat="4.1.6, 5")
# Avoid platforms where the MPI implementation isn't supported
# OpenMPI
platforms = filter(p -> !(p["mpi"] == "openmpi" && arch(p) == "armv6l" && libc(p) == "glibc"), platforms)
# MPItrampoline
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && (Sys.iswindows(p) || libc(p) == "musl")), platforms)
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && Sys.isfreebsd(p)), platforms)

append!(dependencies, platform_dependencies)

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.6", preferred_gcc_version = v"8")
