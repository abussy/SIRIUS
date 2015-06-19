#include <sirius.h>
#include <thread>

using namespace sirius;

void test_rho_sum(double alat, double pw_cutoff, double wf_cutoff, int num_bands, std::vector<int>& mpi_grid)
{
    Communicator comm(MPI_COMM_WORLD);
    BLACS_grid blacs_grid(comm, mpi_grid[0], mpi_grid[1]);

    double a1[] = {alat, 0, 0};
    double a2[] = {0, alat, 0};
    double a3[] = {0, 0, alat};

    Simulation_parameters p;
    Unit_cell uc(p, comm);
    uc.set_lattice_vectors(a1, a2, a3);

    auto& rlv = uc.reciprocal_lattice_vectors();
    auto dims = Utils::find_translation_limits(pw_cutoff, rlv);

    MPI_FFT3D fft(dims, Platform::max_num_threads(), blacs_grid.comm_row());

    auto gv_rho = fft.init_gvec(vector3d<double>(0, 0, 0), pw_cutoff, rlv);
    auto gv_wf = fft.init_gvec(vector3d<double>(0, 0, 0), wf_cutoff, rlv);

    if (comm.rank() == 0)
    {
        printf("MPI grid: %i %i\n", mpi_grid[0], mpi_grid[1]);
        printf("FFT dimensions: %i %i %i\n", fft.size(0), fft.size(1), fft.size(2));
        printf("num_gvec_rho: %i\n", gv_rho.num_gvec_);
        printf("num_gvec_wf: %i\n", gv_wf.num_gvec_);
    }

    mdarray<double, 1> rho(fft.local_size());
    rho.zero();

    splindex<block> spl_gv_wf(gv_wf.num_gvec_, blacs_grid.num_ranks_row(), blacs_grid.rank_row());
    DUMP("spl_gv_wf.local_size: %li, block_size: %li", spl_gv_wf.local_size(), spl_gv_wf.block_size());
    DUMP("num_gvec_loc: %i", gv_wf.num_gvec_loc_);
    DUMP("gvec_offset: %i", gv_wf.gvec_offset_);

    dmatrix<double_complex> psi(gv_wf.num_gvec_, num_bands, blacs_grid, (int)spl_gv_wf.block_size(), 1);
    psi.zero();

    for (int i = 0; i < num_bands; i++) psi.set(i, i, double_complex(1, 0));


    if (psi.num_rows_local() != (int)spl_gv_wf.local_size())
    {
        TERMINATE("wrong index splitting");
    }

    //mdarray<double_complex, 1> buf(gv_wf.num_gvec_loc_);
    mdarray<double_complex, 1> buf(gv_wf.gvec_counts_[0]);

    //== auto scounts = spl_gv_wf.counts();
    //== auto soffsets = spl_gv_wf.offsets();

    //== auto rcounts = &gv_wf.counts_(0);
    //== auto roffsets = &gv_wf.offsets_(0);

    //== if (comm.rank() == 0)
    //== {
    //==     for (int i = 0; i < blacs_grid.comm_row().size(); i++)
    //==     {
    //==         printf("rank: %i, sc, so: %i %i, rc, ro: %i %i\n", i, scounts[i], soffsets[i], rcounts[i], roffsets[i]);

    //==     }
    //== }

    auto a2a_desc = blacs_grid.comm_row().map_alltoall(spl_gv_wf.counts(), gv_wf.gvec_counts_); 


    //std::vector<int> sendcounts(blacs_grid.comm_row().size(), 0);
    //std::vector<int> sdispls(blacs_grid.comm_row().size(), -1);
    //std::vector<int> recvcounts(blacs_grid.comm_row().size(), 0);
    //std::vector<int> rdispls(blacs_grid.comm_row().size(), -1);
    //
    ///* loop over sending ranks */
    //for (int sr = 0; sr < blacs_grid.comm_row().size(); sr++)
    //{
    //    if (!scounts[sr]) continue;

    //    /* beginning of index */
    //    int i0 = soffsets[sr];
    //    /* end of index */
    //    int i1 = soffsets[sr] + scounts[sr] - 1;

    //    PRINT("--------------------");
    //    PRINT("sending rank: %i", sr);
    //    PRINT("i0, i1: %i %i", i0, i1);
    //    PRINT("--------------------");

    //    /* loop over receiving ranks */
    //    for (int rr = 0; rr < blacs_grid.comm_row().size(); rr++)
    //    {
    //        if (!rcounts[rr]) continue;

    //        int j0 = roffsets[rr];
    //        int j1 = roffsets[rr] + rcounts[rr] - 1;
    //        PRINT("j0, j1: %i %i", j0, j1);

    //        /* rank rr recieves nothing from rank sr*/
    //        if (j1 < i0 || i1 < j0) continue;

    //        int s_ofs = std::max(j0 - i0, 0);
    //        int r_ofs = std::max(i0 - j0, 0);
    //        int sz = std::min(i1, j1) - std::max(i0, j0) + 1;
    //        
    //        PRINT("rank: %i receives %i elements to position %i", rr, sz, r_ofs);
    //        if (blacs_grid.comm_row().rank() == sr)
    //        {
    //            sendcounts[rr] = sz;
    //            sdispls[rr] = s_ofs;
    //        }
    //        if (blacs_grid.comm_row().rank() == rr)
    //        {
    //            recvcounts[sr] = sz;
    //            rdispls[sr] = r_ofs;
    //        }
    //    }
    //}

    //int n1(0), n2(0);
    //for (int i = 0; i < blacs_grid.comm_row().size(); i++)
    //{
    //    n1 += sendcounts[i];
    //    n2 += recvcounts[i];
    //}
    //if (n1 != scounts[blacs_grid.rank_row()]) TERMINATE("wrong number of send counts");
    //if (n2 != rcounts[blacs_grid.rank_row()]) TERMINATE("wrong number of receive counts");

               




    Timer t2("sum_rho");
    for (int i = 0; i < psi.num_cols_local(); i++)
    {
        double_complex* inp = (gv_wf.num_gvec_loc_ != 0) ? &buf(0) : NULL;

        //blacs_grid.comm_row().allgather(&psi(0, i), &psi_slice(0), (int)spl_gv_wf.global_offset(), (int)spl_gv_wf.local_size());
        blacs_grid.comm_row().alltoall(&psi(0, i), &a2a_desc.sendcounts[0], &a2a_desc.sdispls[0], inp,
                                       &a2a_desc.recvcounts[0], &a2a_desc.rdispls[0]);
        
        /* this is to make assert statements of mdarray happy */
        int* map = (gv_wf.num_gvec_loc_ != 0) ? &gv_wf.index_map_local_to_local_(0) : NULL;
        fft.input_pw(gv_wf.num_gvec_loc_, map, inp);
        fft.transform(1);
        for (int j = 0; j < (int)fft.local_size(); j++) rho(j) += (std::pow(std::real(fft.buffer(j)), 2) +
                                                                   std::pow(std::imag(fft.buffer(j)), 2));
    }
    t2.stop();

    double nel = 0;
    for (int j = 0; j < (int)fft.local_size(); j++) nel += rho(j);
    comm.allreduce(&nel, 1);
    nel /= fft.size();

    PRINT("num_bands: %i, num_electrons: %f", num_bands, nel);

    




    //#ifdef __PRINT_MEMORY_USAGE
    //MEMORY_USAGE_INFO();
    //#endif

    //std::vector<double_complex> pw_coefs(gv.num_gvec_);
    //
    //if (comm.rank() == 0) printf("num_gvec: %i\n", gv.num_gvec_);
    //DUMP("num_gvec_loc: %i", gv.num_gvec_loc_);

    //int n = (fft.size() < 100000) ? gv.num_gvec_ : std::min(50, gv.num_gvec_);

    //for (int ig = 0; ig < n; ig++)
    //{
    //    memset(&pw_coefs[0], 0, pw_coefs.size() * sizeof(double_complex));
    //    auto gvec = fft.gvec_by_index(gv.gvec_index_(ig));

    //    if (comm.rank() == 0) printf("G: %i %i %i\n", gvec[0], gvec[1], gvec[2]);

    //    pw_coefs[ig] = 1.0;
    //    fft.input_pw(gv.num_gvec_loc_, &gv.index_map_local_to_local_(0), &pw_coefs[gv.gvec_offset_]);
    //    fft.transform(1);

    //    mdarray<double_complex, 3> tmp(&fft.buffer(0), fft.size(0), fft.size(1), fft.local_size_z());

    //    double dmax = 0;
    //    // loop over 3D array (real space)
    //    for (int j0 = 0; j0 < fft.size(0); j0++)
    //    {
    //        for (int j1 = 0; j1 < fft.size(1); j1++)
    //        {
    //            for (int j2 = 0; j2 < fft.local_size_z(); j2++)
    //            {
    //                // get real space fractional coordinate
    //                vector3d<double> fv(double(j0) / fft.size(0), 
    //                                    double(j1) / fft.size(1), 
    //                                    double(j2 + fft.offset_z()) / fft.size(2));
    //                double_complex z = std::exp(twopi * double_complex(0.0, (fv * gvec)));
    //                dmax = std::max(dmax, std::abs(tmp(j0, j1, j2) - z));
    //            }
    //        }
    //    }
    //    comm.allreduce<double, op_max>(&dmax, 1);
    //    if (dmax > 1e-12)
    //    {
    //        printf("maximum difference: %18.12f\n", dmax);
    //        exit(-1);
    //    }
    //}
    //printf("OK\n");
}

int main(int argn, char **argv)
{
    cmd_args args;
    args.register_key("--help", "print this help and exit");
    args.register_key("--alat=", "{double} lattice constant");
    args.register_key("--pw_cutoff=", "{double} plane-wave cutoff [a.u.^-1]");
    args.register_key("--wf_cutoff=", "{double} wave-function cutoff [a.u.^-1]");
    args.register_key("--num_bands=", "{int} number of bands");
    args.register_key("--mpi_grid=", "{vector<int>} MPI grid");

    args.parse_args(argn, argv);
    if (args.exist("help"))
    {
        printf("Usage: %s [options]\n", argv[0]);
        args.print_help();
        exit(0);
    }

    double alat = 7.0;
    double pw_cutoff = 16.0;
    double wf_cutoff = 5.0;
    int num_bands = 10;
    std::vector<int> mpi_grid;

    alat = args.value<double>("alat", alat);
    pw_cutoff = args.value<double>("pw_cutoff", pw_cutoff);
    wf_cutoff = args.value<double>("wf_cutoff", wf_cutoff);
    num_bands = args.value<int>("num_bands", num_bands);
    mpi_grid = args.value< std::vector<int> >("mpi_grid", std::vector<int>({1, 1}));


    Platform::initialize(1);

    test_rho_sum(alat, pw_cutoff, wf_cutoff, num_bands, mpi_grid);

    //#ifdef __PRINT_MEMORY_USAGE
    //MEMORY_USAGE_INFO();
    //#endif
    Timer::print();

    Platform::finalize();
}
