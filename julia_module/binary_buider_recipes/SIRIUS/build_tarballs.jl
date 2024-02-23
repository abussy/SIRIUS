# Note that this script can accept some limited command-line arguments, run                          
# `julia build_tarballs.jl --help` to see a usage message.                                           
using BinaryBuilder, Pkg                                                                             
using Base.BinaryPlatforms                                                                           
const YGGDRASIL_DIR = "../.."                                                                        
include(joinpath(YGGDRASIL_DIR, "platforms", "mpi.jl"))

name = "SIRIUS"
version = v"7.5.2"

sources = [
   #GitSource("https://github.com/abussy/SIRIUS/", "c45f03811e439c1a49858bcacd849a26932b2cf8")
   GitSource("https://github.com/abussy/SIRIUS/", "6e6b337585e6f37fc09fe7e73dcd319af45e3b3d")
]

script = raw"""
cd $WORKSPACE/srcdir

#Hack, because supremely hard to change CMake version in here
sed -i 's/VERSION 3.23/VERSION 3.21.7/' ./CMakeLists.txt
sed -i 's/sirius PUBLIC cuda_std_14/ sirius PUBLIC cxx_std_14/' ./src/CMakeLists.txt

mkdir build
cd build

#For GSL to be linked correctly to cblas
export LDFLAGS="-lgsl -lgslcblas -lblastrampoline"

#if [[ "$nbits" == "64" ]]; then
#    OPENBLAS_LIB="$libdir/libopenblas64_.$dlext"
#    
#    # Fix suffixes for 64-bit OpenBLAS
#    SYMB_DEFS=()
#    SYMBOLS=(ssyevr cheevd dsyevr zheevd ssyevd dgesv dgetrf dsyevd ilaenv dsytrf dsygvx ssygvx 
#             zhemm chegvx dlamch zhegvx dlartg dgetri zheevx zgetrf cheevx dsytrs zgetri zgemm
#             zpotrf ztrtri dscal dgemm ztrmm dtrmm dtrtri dpotrf)
#    for sym in ${SYMBOLS[@]}; do
#        SYMB_DEFS+=("-D${sym}_=${sym}_64_")
#    done
#    export CXXFLAGS="${SYMB_DEFS[@]}"
#
#else
#    OPENBLAS_LIB="$libdir/libopenblas.$dlext"
#fi

CMAKE_ARGS="-DSIRIUS_CREATE_FORTRAN_BINDINGS=ON \
            -DSIRIUS_USE_OPENMP=ON \
            -DSIRIUS_USE_MEMORY_POOL=OFF \
            -DSIRIUS_BUILD_DOCS=OFF \
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

#somehow need to run cmake twice for MPI to work
cmake .. ${CMAKE_ARGS} || cmake .. ${CMAKE_ARGS}

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
   LibraryProduct("libsirius", :libsirius)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("GSL_jll"), 
    #Dependency("OpenBLAS32_jll"), 
    Dependency("libblastrampoline_jll"), 
    Dependency("Libxc_jll"), 
    Dependency("HDF5_jll"),
    Dependency("Doxygen_jll"), 
    Dependency("spglib_jll"), 
    Dependency(PackageSpec(; name = "spla_jll",  uuid = "c0b6e1fa-1634-5519-836b-ea1958e971dd", 
                           path = "/home/bussya/.julia/dev/spla_jll")),
    Dependency(PackageSpec(; name = "SpFFT_jll",  uuid = "8be71ea8-ca6b-526f-84e8-fb2862ef466b", 
                           path = "/home/bussya/.julia/dev/SpFFT_jll")),
    Dependency(PackageSpec(; name = "COSTA_jll",  uuid = "5f6e5a8a-773c-5339-a5d0-b2b35bd1df39",
                           path = "/home/bussya/.julia/dev/COSTA_jll")),
    Dependency("CompilerSupportLibraries_jll"), 
    Dependency("LLVMOpenMP_jll", platforms=filter(Sys.isapple, platforms))
]

platforms, platform_dependencies = MPI.augment_platforms(platforms; MPItrampoline_compat="5.3.0", OpenMPI_compat="4.1.6, 5")
# Avoid platforms where the MPI implementation isn't supported
# OpenMPI
platforms = filter(p -> !(p["mpi"] == "openmpi" && arch(p) == "armv6l" && libc(p) == "glibc"), platforms)
# MPItrampoline
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && (Sys.iswindows(p) || libc(p) == "musl")), platforms)
platforms = filter(p -> !(p["mpi"] == "mpitrampoline" && Sys.isfreebsd(p)), platforms)

###TODO: temporary, only mpich for testing
#platforms = filter(p -> !(p["mpi"] == "openmpi"), platforms)
#platforms = filter(p -> !(p["mpi"] == "mpitrampoline"), platforms)

append!(dependencies, platform_dependencies)

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.9", preferred_gcc_version = v"10")
