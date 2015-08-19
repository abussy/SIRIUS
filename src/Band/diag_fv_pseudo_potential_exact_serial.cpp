#include "band.h"

namespace sirius {

void Band::diag_fv_pseudo_potential_exact_serial(K_point* kp__,
                                                 std::vector<double>& veff_it_coarse__)
{
    STOP();

//    /* cache kinetic energy */
//    std::vector<double> pw_ekin = kp__->get_pw_ekin();
//
//    /* short notation for target wave-functions */
//    mdarray<double_complex, 2>& psi = kp__->fv_states_slab();
//
//    /* short notation for number of target wave-functions */
//    int num_bands = parameters_.num_fv_states();     
//
//    int ngk = kp__->num_gkvec();
//
//    mdarray<double_complex, 2> phi(ngk, ngk);
//    mdarray<double_complex, 2> hphi(ngk, ngk);
//    mdarray<double_complex, 2> ophi(ngk, ngk);
//    mdarray<double_complex, 1> kappa(ngk * ngk);
//    
//    std::vector<double> eval(ngk);
//
//    phi.zero();
//    for (int i = 0; i < ngk; i++) phi(i, i) = complex_one;
//
//    /* offset in the packed array of on-site matrices */
//    mdarray<int, 1> packed_mtrx_offset(unit_cell_.num_atoms());
//    int packed_mtrx_size = 0;
//    for (int ia = 0; ia < unit_cell_.num_atoms(); ia++)
//    {   
//        int nbf = unit_cell_.atom(ia)->mt_basis_size();
//        packed_mtrx_offset(ia) = packed_mtrx_size;
//        packed_mtrx_size += nbf * nbf;
//    }
//    
//    /* pack Q and D matrices */
//    mdarray<double_complex, 1> d_mtrx_packed(packed_mtrx_size);
//    mdarray<double_complex, 1> q_mtrx_packed(packed_mtrx_size);
//
//    for (int ia = 0; ia < unit_cell_.num_atoms(); ia++)
//    {
//        int nbf = unit_cell_.atom(ia)->mt_basis_size();
//        for (int xi2 = 0; xi2 < nbf; xi2++)
//        {
//            for (int xi1 = 0; xi1 < nbf; xi1++)
//            {
//                d_mtrx_packed(packed_mtrx_offset(ia) + xi2 * nbf + xi1) = unit_cell_.atom(ia)->d_mtrx(xi1, xi2);
//                q_mtrx_packed(packed_mtrx_offset(ia) + xi2 * nbf + xi1) = unit_cell_.atom(ia)->type()->uspp().q_mtrx(xi1, xi2);
//            }
//        }
//    }
//    
//    apply_h_o_serial(kp__, veff_it_coarse__, pw_ekin, 0, ngk, phi, hphi, ophi, kappa, packed_mtrx_offset,
//                         d_mtrx_packed, q_mtrx_packed);
//        
//    gen_evp_solver()->solve(ngk, num_bands, num_bands, num_bands, hphi.at<CPU>(), hphi.ld(), ophi.at<CPU>(), ophi.ld(), 
//                            &eval[0], psi.at<CPU>(), psi.ld());
//
//    kp__->set_fv_eigen_values(&eval[0]);
}

};
