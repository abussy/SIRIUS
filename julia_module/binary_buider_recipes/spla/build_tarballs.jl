# Note that this script can accept some limited command-line arguments, run                          
# `julia build_tarballs.jl --help` to see a usage message.                                           
using BinaryBuilder, Pkg                                                                             
using Base.BinaryPlatforms                                                                           
const YGGDRASIL_DIR = "../.."                                                                        
include(joinpath(YGGDRASIL_DIR, "platforms", "mpi.jl"))

name = "spla"
version = v"1.5.5"

sources = [
   GitSource("https://github.com/eth-cscs/spla/", "6fe85e49069ae287e5ec3cfe5487720d85ebe97a")
]

script = raw"""
cd $WORKSPACE/srcdir

mkdir build
cd build

export LDFLAGS="-lblastrampoline"


#if [[ "$nbits" == "64" ]]; then
#    OPENBLAS_LIB="$libdir/libopenblas64_.$dlext"
#
#    # Fix suffixes for 64-bit OpenBLAS
#    SYMB_DEFS=()
#    SYMBOLS=(cblas_sgemm cblas_dgemm cblas_cgemm cblas_zgemm openblas_get_num_threads
#             openblas_set_num_threads openblas_get_parallel)
#    for sym in ${SYMBOLS[@]}; do
#        SYMB_DEFS+=("-D${sym}=${sym}64_")
#    done
#    export CXXFLAGS="${SYMB_DEFS[@]}"
#else
#    OPENBLAS_LIB="$libdir/libopenblas.$dlext"
#fi

CMAKE_ARGS="-DSPLA_OMP=ON \
            -DSPLA_HOST_BLAS=AUTO \
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
            -DCMAKE_FIND_ROOT_PATH='${prefix}/lib/mpich;${prefix}' \
            -DCMAKE_INSTALL_PREFIX=$prefix \
            -DCMAKE_BUILD_TYPE=Release \
            -DBUILD_SHARED_LIBS=ON \
            -DMPI_C_COMPILER=$bindir/mpicc \
            -DMPI_CXX_COMPILER=$bindir/mpicxx"

#TODO: the following only works with MPItrampoline, and therefore not with HDF5
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
   LibraryProduct("libspla", :libspla)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    #Dependency("OpenBLAS32_jll"),
    Dependency("libblastrampoline_jll"),
    Dependency("CompilerSupportLibraries_jll", platforms=filter(!Sys.isapple, platforms)),
    Dependency("LLVMOpenMP_jll", platforms=filter(Sys.isapple, platforms))
]

platforms, platform_dependencies = MPI.augment_platforms(platforms)
# Avoid platforms where the MPI implementation isn't supported
# OpenMPI
platforms = filter(p -> !(p["mpi"] == "openmpi" && arch(p) == "armv6l" && libc(p) == "glibc"), platforms)
# MPItrampoline
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && (Sys.iswindows(p) || libc(p) == "musl")), platforms)
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && Sys.isfreebsd(p)), platforms)

append!(dependencies, platform_dependencies)

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.9", preferred_gcc_version = v"7")
