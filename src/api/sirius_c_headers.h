#include <stdbool.h>
#include <complex.h>

void
sirius_initialize(bool const* call_mpi_init__, int* error_code__);

void
sirius_finalize(bool const* call_mpi_fin__, bool const* call_device_reset__, bool const* call_fftw_fin__,
                int* error_code__);

void
sirius_start_timer(char const* name__, int* error_code__);

void
sirius_stop_timer(char const* name__, int* error_code__);

void
sirius_print_timers(bool* flatten__, int* error_code__);

void
sirius_serialize_timers(char const* fname__, int* error_code__);

void
sirius_context_initialized(void* const* handler__, bool* status__, int* error_code__);

void
sirius_create_context(int fcomm__, void** handler__, int* fcomm_k__, int* fcomm_band__, int* error_code__);

void
sirius_import_parameters(void* const* handler__, char const* str__, int* error_code__);

void
sirius_set_parameters(void* const* handler__, int const* lmax_apw__, int const* lmax_rho__, int const* lmax_pot__,
                      int const* num_fv_states__, int const* num_bands__, int const* num_mag_dims__,
                      double const* pw_cutoff__, double const* gk_cutoff__, int const* fft_grid_size__,
                      int const* auto_rmt__, bool const* gamma_point__, bool const* use_symmetry__,
                      bool const* so_correction__, char const* valence_rel__, char const* core_rel__,
                      double const* iter_solver_tol_empty__, char const* iter_solver_type__, int const* verbosity__,
                      bool const* hubbard_correction__, int const* hubbard_correction_kind__,
                      bool const* hubbard_full_orthogonalization__, bool const* hubbard_constrained_calculation__,
                      char const* hubbard_orbitals__, int const* sht_coverage__, double const* min_occupancy__,
                      char const* smearing__, double const* smearing_width__, double const* spglib_tol__,
                      char const* electronic_structure_method__, int* error_code__);

void
sirius_get_parameters(void* const* handler__, int* lmax_apw__, int* lmax_rho__, int* lmax_pot__, int* num_fv_states__,
                      int* num_bands__, int* num_spins__, int* num_mag_dims__, double* pw_cutoff__, double* gk_cutoff__,
                      int* fft_grid_size__, int* auto_rmt__, bool* gamma_point__, bool* use_symmetry__,
                      bool* so_correction__, double* iter_solver_tol__, double* iter_solver_tol_empty__,
                      int* verbosity__, bool* hubbard_correction__, double* evp_work_count__, int* num_loc_op_applied__,
                      int* num_sym_op__, char* electronic_structure_method__, int* error_code__);

void
sirius_add_xc_functional(void* const* handler__, char const* name__, int* error_code__);

void
sirius_set_mpi_grid_dims(void* const* handler__, int const* ndims__, int const* dims__, int* error_code__);

void
sirius_set_lattice_vectors(void* const* handler__, double const* a1__, double const* a2__, double const* a3__,
                           int* error_code__);

void
sirius_initialize_context(void* const* handler__, int* error_code__);

void
sirius_update_context(void* const* handler__, int* error_code__);

void
sirius_print_info(void* const* handler__, int* error_code__);

void
sirius_free_object_handler(void** handler__, int* error_code__);

void
sirius_set_periodic_function_ptr(void* const* handler__, char const* label__, double* f_mt__, int const* lmmax__,
                                 int const* nrmtmax__, int const* num_atoms__, double* f_rg__, int const* size_x__,
                                 int const* size_y__, int const* size_z__, int const* offset_z__, int* error_code__);

void
sirius_set_periodic_function(void* const* gs_handler__, char const* label__, double* f_mt__, int const* lmmax__,
                             int const* nrmtmax__, int const* num_atoms__, double* f_rg__, int const* size_x__,
                             int const* size_y__, int const* size_z__, int const* offset_z__, int* error_code__);

void
sirius_get_periodic_function(void* const* gs_handler__, char const* label__, double* f_mt__, int const* lmmax__,
                             int const* nrmtmax__, int const* num_atoms__, double* f_rg__, int const* size_x__,
                             int const* size_y__, int const* size_z__, int const* offset_z__, int* error_code__);

void
sirius_create_kset(void* const* handler__, int const* num_kpoints__, double* kpoints__, double const* kpoint_weights__,
                   bool const* init_kset__, void** kset_handler__, int* error_code__);

void
sirius_create_kset_from_grid(void* const* handler__, int const* k_grid__, int const* k_shift__,
                             bool const* use_symmetry, void** kset_handler__, int* error_code__);

void
sirius_create_ground_state(void* const* ks_handler__, void** gs_handler__, int* error_code__);

void
sirius_initialize_kset(void* const* ks_handler__, int* count__, int* error_code__);

void
sirius_find_ground_state(void* const* gs_handler__, double const* density_tol__, double const* energy_tol__,
                         double const* iter_solver_tol__, bool const* initial_guess__, int const* max_niter__,
                         bool const* save_state__, bool* converged__, int* niter__, double* rho_min__,
                         int* error_code__);

void
sirius_check_scf_density(void* const* gs_handler__, int* error_code__);

void
sirius_update_ground_state(void** gs_handler__, int* error_code__);

void
sirius_add_atom_type(void* const* handler__, char const* label__, char const* fname__, int const* zn__,
                     char const* symbol__, double const* mass__, bool const* spin_orbit__, int* error_code__);

void
sirius_set_atom_type_radial_grid(void* const* handler__, char const* label__, int const* num_radial_points__,
                                 double const* radial_points__, int* error_code__);

void
sirius_set_atom_type_radial_grid_inf(void* const* handler__, char const* label__, int const* num_radial_points__,
                                     double const* radial_points__, int* error_code__);

void
sirius_add_atom_type_radial_function(void* const* handler__, char const* atom_type__, char const* label__,
                                     double const* rf__, int const* num_points__, int const* n__, int const* l__,
                                     int const* idxrf1__, int const* idxrf2__, double const* occ__, int* error_code__);

void
sirius_set_atom_type_hubbard(void* const* handler__, char const* label__, int const* l__, int const* n__,
                             double const* occ__, double const* U__, double const* J__, double const* alpha__,
                             double const* beta__, double const* J0__, int* error_code__);

void
sirius_set_atom_type_dion(void* const* handler__, char const* label__, int const* num_beta__, double* dion__,
                          int* error_code__);

void
sirius_set_atom_type_paw(void* const* handler__, char const* label__, double const* core_energy__,
                         double const* occupations__, int const* num_occ__, int* error_code__);

void
sirius_add_atom(void* const* handler__, char const* label__, double const* position__, double const* vector_field__,
                int* error_code__);

void
sirius_set_atom_position(void* const* handler__, int const* ia__, double const* position__, int* error_code__);

void
sirius_set_pw_coeffs(void* const* gs_handler__, char const* label__, double complex const* pw_coeffs__,
                     bool const* transform_to_rg__, int const* ngv__, int* gvl__, int const* comm__, int* error_code__);

void
sirius_get_pw_coeffs(void* const* gs_handler__, char const* label__, double complex* pw_coeffs__, int const* ngv__,
                     int* gvl__, int const* comm__, int* error_code__);

void
sirius_initialize_subspace(void* const* gs_handler__, void* const* ks_handler__, int* error_code__);

void
sirius_find_eigen_states(void* const* gs_handler__, void* const* ks_handler__, bool const* precompute_pw__,
                         bool const* precompute_rf__, bool const* precompute_ri__, double const* iter_solver_tol__,
                         int const* iter_solver_steps__, int* error_code__);

void
sirius_generate_initial_density(void* const* gs_handler__, int* error_code__);

void
sirius_generate_effective_potential(void* const* gs_handler__, int* error_code__);

void
sirius_generate_density(void* const* gs_handler__, bool const* add_core__, bool const* transform_to_rg__,
                        bool const* paw_only__, int* error_code__);

void
sirius_set_band_occupancies(void* const* ks_handler__, int const* ik__, int const* ispn__,
                            double const* band_occupancies__, int* error_code__);

void
sirius_get_band_occupancies(void* const* ks_handler__, int const* ik__, int const* ispn__, double* band_occupancies__,
                            int* error_code__);

void
sirius_get_band_energies(void* const* ks_handler__, int const* ik__, int const* ispn__, double* band_energies__,
                         int* error_code__);

void
sirius_get_energy(void* const* gs_handler__, char const* label__, double* energy__, int* error_code__);

void
sirius_get_forces(void* const* gs_handler__, char const* label__, double* forces__, int* error_code__);

void
sirius_get_stress_tensor(void* const* gs_handler__, char const* label__, double* stress_tensor__, int* error_code__);

void
sirius_get_num_beta_projectors(void* const* handler__, char const* label__, int* num_bp__, int* error_code__);

void
sirius_get_wave_functions(void* const* ks_handler__, double const* vkl__, int const* spin__, int const* num_gvec_loc__,
                          int const* gvec_loc__, double complex* evec__, int const* ld__,
                          int const* num_spin_comp__, int* error_code__);

void
sirius_add_atom_type_aw_descriptor(void* const* handler__, char const* label__, int const* n__, int const* l__,
                                   double const* enu__, int const* dme__, bool const* auto_enu__, int* error_code__);

void
sirius_add_atom_type_lo_descriptor(void* const* handler__, char const* label__, int const* ilo__, int const* n__,
                                   int const* l__, double const* enu__, int const* dme__, bool const* auto_enu__,
                                   int* error_code__);

void
sirius_set_atom_type_configuration(void* const* handler__, char const* label__, int const* n__, int const* l__,
                                   int const* k__, double const* occupancy__, bool const* core__, int* error_code__);

void
sirius_generate_coulomb_potential(void* const* gs_handler__, double* vh_el__, int* error_code__);

void
sirius_generate_xc_potential(void* const* gs_handler__, int* error_code__);

void
sirius_get_kpoint_inter_comm(void* const* handler__, int* fcomm__, int* error_code__);

void
sirius_get_kpoint_inner_comm(void* const* handler__, int* fcomm__, int* error_code__);

void
sirius_get_fft_comm(void* const* handler__, int* fcomm__, int* error_code__);

void
sirius_get_num_gvec(void* const* handler__, int* num_gvec__, int* error_code__);

void
sirius_get_gvec_arrays(void* const* handler__, int* gvec__, double* gvec_cart__, double* gvec_len__,
                       int* index_by_gvec__, int* error_code__);

void
sirius_get_num_fft_grid_points(void* const* handler__, int* num_fft_grid_points__, int* error_code__);

void
sirius_get_fft_index(void* const* handler__, int* fft_index__, int* error_code__);

void
sirius_get_max_num_gkvec(void* const* ks_handler__, int* max_num_gkvec__, int* error_code__);

void
sirius_get_gkvec_arrays(void* const* ks_handler__, int* ik__, int* num_gkvec__, int* gvec_index__, double* gkvec__,
                        double* gkvec_cart__, double* gkvec_len, double* gkvec_tp__, int* error_code__);

void
sirius_get_step_function(void* const* handler__, double complex* cfunig__, double* cfunrg__, int* num_rg_points__,
                         int* error_code__);

void
sirius_set_h_radial_integrals(void* const* handler__, int* ia__, int* lmmax__, double* val__, int* l1__, int* o1__,
                              int* ilo1__, int* l2__, int* o2__, int* ilo2__, int* error_code__);

void
sirius_set_o_radial_integral(void* const* handler__, int* ia__, double* val__, int* l__, int* o1__, int* ilo1__,
                             int* o2__, int* ilo2__, int* error_code__);

void
sirius_set_o1_radial_integral(void* const* handler__, int* ia__, double* val__, int* l1__, int* o1__, int* ilo1__,
                              int* l2__, int* o2__, int* ilo2__, int* error_code__);

void
sirius_set_radial_function(void* const* handler__, int const* ia__, int const* deriv_order__, double const* f__,
                           int const* l__, int const* o__, int const* ilo__, int* error_code__);

void
sirius_set_equivalent_atoms(void* const* handler__, int* equivalent_atoms__, int* error_code__);

void
sirius_update_atomic_potential(void* const* gs_handler__, int* error_code__);

void
sirius_option_get_number_of_sections(int* length__, int* error_code__);

void
sirius_option_get_section_name(int elem__, char* section_name__, int section_name_length__, int* error_code__);

void
sirius_option_get_section_length(char const* section__, int* length__, int* error_code__);

void
sirius_option_get_info(char const* section__, int elem__, char* key_name__, int key_name_len__, int* type__,
                       int* length__, int* enum_size__, char* title__, int title_len__, char* description__,
                       int description_len__, int* error_code__);

void
sirius_option_get(char const* section__, char const* name__, int const* type__, void* data_ptr__,
                  int const* max_length__, int const* enum_idx__, int* error_code__);

void
sirius_option_set(void* const* handler__, char const* section__, char const* name__, int const* type__,
                  void const* data_ptr__, int const* max_length__, bool const* append__, int* error_code__);

void
sirius_dump_runtime_setup(void* const* handler__, char* filename__, int* error_code__);

void
sirius_get_fv_eigen_vectors(void* const* ks_handler__, int const* ik__, double complex* fv_evec__, int const* ld__,
                            int const* num_fv_states__, int* error_code__);

void
sirius_get_fv_eigen_values(void* const* ks_handler__, int const* ik__, double* fv_eval__, int const* num_fv_states__,
                           int* error_code__);

void
sirius_get_sv_eigen_vectors(void* const* ks_handler__, int const* ik__, double complex* sv_evec__,
                            int const* num_bands__, int* error_code__);

void
sirius_set_rg_values(void* const* gs_handler__, char const* label__, int const* grid_dims__, int const* local_box_origin__,
                     int const* local_box_size__, int const* fcomm__, double const* values__,
                     bool const* transform_to_pw__, int* error_code__);

void
sirius_get_rg_values(void* const* gs_handler__, char const* label__, int const* grid_dims__, int const* local_box_origin__,
                     int const* local_box_size__, int const* fcomm__, double* values__, bool const* transform_to_rg__,
                     int* error_code__);

void
sirius_get_total_magnetization(void* const* gs_handler__, double* mag__, int* error_code__);

void
sirius_get_num_kpoints(void* const* ks_handler__, int* num_kpoints__, int* error_code__);

void
sirius_get_kpoint_properties(void* const* ks_handler__, int const* ik__, double* weight__, double* coordinates__,
                             int* error_code__);

void
sirius_set_callback_function(void* const* handler__, char const* label__, void (*fptr__)(), int* error_code__);

void
sirius_nlcg(void* const* gs_handler__, void* const* ks_handler__, int* error_code__);

void
sirius_nlcg_params(void* const* gs_handler__, void* const* ks_handler__, double const* temp__, char const* smearing__,
                   double const* kappa__, double const* tau__, double const* tol__, int const* maxiter__,
                   int const* restart__, char const* processing_unit__, bool* converged__, int* error_code__);

void
sirius_add_hubbard_atom_pair(void* const* handler__, int* const atom_pair__, int* const translation__, int* const n__,
                             int* const l__, const double* const coupling__, int* error_code__);

void
sirius_set_hubbard_contrained_parameters(void* const* handler__, double const* hubbard_conv_thr__,
                                         double const* hubbard_mixing_beta__, double const* hubbard_strength__,
                                         int const* hubbard_maxstep__, char const* hubbard_constraint_type__,
                                         int* const error_code__);

void
sirius_add_hubbard_atom_constraint(void* const* handler__, int* const atom_id__, int* const n__, int* const l__,
                                   int* const lmax_at__, const double* const occ__, int* const orbital_order__,
                                   int* const error_code__);

void
sirius_create_H0(void* const* gs_handler__, int* error_code__);

void
sirius_linear_solver(void* const* gs_handler__, double const* vkq__, int const* num_gvec_kq_loc__,
                     int const* gvec_kq_loc__, double complex* dpsi__, double complex* psi__,
                     double* eigvals__, double complex* dvpsi__, int const* ld__, int const* num_spin_comp__,
                     double const* alpha_pv__, int const* spin__, int const* nbnd_occ__, double const* tol__,
                     int* niter__, int* error_code__);

void
sirius_generate_rhoaug_q(void* const* gs_handler__, int const* iat__, int const* num_atoms__, int const* num_gvec_loc__,
                         int const* num_spin_comp__, double complex const* qpw__, int const* ldq__,
                         double complex const* phase_factors_q__, int const* mill__,
                         double complex const* dens_mtrx__, int const* ldd__, double complex* rho_aug__,
                         int* error_code__);

void
sirius_generate_d_operator_matrix(void* const* gs_handler__, int* error_code__);

void
sirius_save_state(void** gs_handler__, const char* file_name__, int* error_code__);

void
sirius_load_state(void** gs_handler__, const char* file_name__, int* error_code__);

void
sirius_set_density_matrix(void** gs_handler__, int const* ia__, double complex* dm__, int const* ld__,
                          int* error_code__);

void
sirius_get_major_version(int* version);

void
sirius_get_minor_version(int* version);

void
sirius_get_revision(int* version);

void
sirius_is_initialized(bool* status__, int* error_code__);

void
sirius_create_context_from_json(int fcomm__, void** handler__, char const* fname__, int* error_code__);

void
sirius_create_context_from_json_commworld(void** handler__, char const* fname__, int* error_code__);

void
sirius_get_num_atoms(void* const* gs_handler__, int* num_atoms__, int* error_code__);

void
sirius_get_kp_params_from_ctx(void* const* handler__, int* k_grid__, int* k_shift__, bool* use_symmetry__,
                              int* error_code__);

void
sirius_get_scf_params_from_ctx(void* const* handler__, double* density_tol__, double* energy_tol__,
                               double* iter_solver_tol__, int* max_niter, int* error_code__);

void
sirius_get_nlcg_params_from_ctx(void* const* handler__, double* temp__, char* smearing__, double* kappa__,
                                double* tau__, double* tol__, int* maxiter__, int* restart__, 
                                char* processing_unit__, int* error_code__);

void
sirius_get_fft_local_z_offset(void* const* handler__, int* local_z_offset__, int* error_code__);

void
sirius_create_hamiltonian(void* const* gs_handler__, void** H0_handler__, int* error_code__);

void 
sirius_diagonalize_hamiltonian(void* const* gs_handler__, void* const* H0_handler__, 
                                    double* const iter_solver_tol__, int* const max_steps__, 
                                    bool* converged__, int* niter__, int* error_code__);

void
sirius_find_band_occupancies(void* const* ks_handler__, int* error_code__);

void
sirius_set_num_bands(void* const* handler__, int* const num_bands__, int* error_code__);

void
sirius_fft_transform(void* const* gs_handler__, char const* label__, int* direction__, int* error_code__);

void
sirius_get_psi(void* const* ks_handler__, int* ik__, int* ispin__, double complex* psi__, 
               int* error_code__);

void
sirius_get_gkvec(void* const* ks_handler__, int* ik__, double* gvec__, int* error_code__);

