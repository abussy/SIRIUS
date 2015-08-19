// Copyright (c) 2013-2014 Anton Kozhevnikov, Thomas Schulthess
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
// the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the 
//    following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
//    and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR 
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/** \file periodic_function.h
 *   
 *  \brief Contains declaration and partial implementation of sirius::Periodic_function class.
 */

#ifndef __PERIODIC_FUNCTION_H__
#define __PERIODIC_FUNCTION_H__

#include "simulation_context.h"
#include "mdarray.h"
#include "spheric_function.h"
#include "mixer.h"

// TODO: this implementation is better, however the distinction between local and global periodic functions is
//       still not very clear
namespace sirius
{

/// Representation of the periodical function on the muffin-tin geometry.
/** Inside each muffin-tin the spherical expansion is used:
 *   \f[
 *       f({\bf r}) = \sum_{\ell m} f_{\ell m}(r) Y_{\ell m}(\hat {\bf r})
 *   \f]
 *   or
 *   \f[
 *       f({\bf r}) = \sum_{\ell m} f_{\ell m}(r) R_{\ell m}(\hat {\bf r})
 *   \f]
 *   In the interstitial region function is stored on the real-space grid or as a Fourier series:
 *   \f[
 *       f({\bf r}) = \sum_{{\bf G}} f({\bf G}) e^{i{\bf G}{\bf r}}
 *   \f]
 *
 *   The following terminology is used to describe the distribution of the function:
 *       - global function: the whole copy of the function is stored on each MPI rank. Ranks should take care about the
 *         syncronization of the data.
 *       - local function: the function is distributed across the MPI ranks. 
 *
 *   \note In order to check if the function is defined as global or as distributed, check the f_mt_ and f_it_ pointers.
 *         If the function is global, the pointers should not be null.
 */
template<typename T> 
class Periodic_function
{ 
    protected:

        /* forbid copy constructor */
        Periodic_function(const Periodic_function<T>& src) = delete;
        
        /* forbid assigment operator */
        Periodic_function<T>& operator=(const Periodic_function<T>& src) = delete;

    private:
       
        /// Complex counterpart for a given type T.
        typedef typename type_wrapper<T>::complex_t complex_t; 

        Simulation_parameters const& parameters_;
        
        Unit_cell const& unit_cell_;

        Step_function const* step_function_;

        Communicator const& comm_;

        /// Alias for FFT driver.
        FFT3D<CPU>* fft_;

        Gvec const& gvec_;

        /// Local part of muffin-tin functions.
        mdarray<Spheric_function<spectral, T>, 1> f_mt_local_;
        
        /// Global muffin-tin array 
        mdarray<T, 3> f_mt_;

        /// Interstitial part of the function.
        mdarray<T, 1> f_it_;

        /// plane-wave expansion coefficients
        mdarray<complex_t, 1> f_pw_;

        int angular_domain_size_;
        
        /// number of plane-wave expansion coefficients
        int num_gvec_;

        //splindex<block> spl_fft_size_;

        /// Set pointer to local part of muffin-tin functions
        void set_local_mt_ptr()
        {
            for (int ialoc = 0; ialoc < (int)unit_cell_.spl_num_atoms().local_size(); ialoc++)
            {
                int ia = unit_cell_.spl_num_atoms(ialoc);
                f_mt_local_(ialoc) = Spheric_function<spectral, T>(&f_mt_(0, 0, ia), angular_domain_size_, unit_cell_.atom(ia)->radial_grid());
            }
        }
        
    public:

        /// Constructor
        Periodic_function(Simulation_context& ctx__,
                          int angular_domain_size__,
                          bool alloc_pw__ = true);
        
        /// Destructor
        ~Periodic_function();
        
        /// Allocate memory
        void allocate(bool allocate_global_mt);

        /// Zero the function.
        void zero();
        
        /// Syncronize global function.
        void sync(bool sync_mt, bool sync_it);

        /// Copy from source
        //void copy(Periodic_function<T>* src);
        inline void copy_to_global_ptr(T* f_mt__, T* f_it__);

        /// Add the function
        void add(Periodic_function<T>* g);

        T integrate(std::vector<T>& mt_val, T& it_val);

        template <index_domain_t index_domain>
        inline T& f_mt(int idx0, int idx1, int ia);
        
        /** \todo write and read distributed functions */
        void hdf5_write(HDF5_tree h5f);

        void hdf5_read(HDF5_tree h5f);

        size_t size();

        size_t pack(size_t offset, Mixer<double>* mixer);
        
        size_t unpack(T const* array);
       
        /// Set the global pointer to the muffin-tin part
        void set_mt_ptr(T* mt_ptr__)
        {
            f_mt_ = mdarray<T, 3>(mt_ptr__, angular_domain_size_, unit_cell_.max_num_mt_points(), unit_cell_.num_atoms());
            set_local_mt_ptr();
        }

        /// Set the global pointer to the interstitial part
        void set_it_ptr(T* it_ptr__)
        {
            f_it_ = mdarray<T, 1>(it_ptr__, fft_->local_size());
            //set_local_it_ptr();
        }

        inline Spheric_function<spectral, T>& f_mt(int ialoc)
        {
            return f_mt_local_(ialoc);
        }

        inline Spheric_function<spectral, T> const& f_mt(int ialoc) const
        {
            return f_mt_local_(ialoc);
        }

        //template <index_domain_t index_domain>
        inline T& f_it(int ir)
        {
            return f_it_(ir);
        }

        //template <index_domain_t index_domain>
        inline T const& f_it(int ir) const
        {
            return f_it_(ir);
        }
        
        inline complex_t& f_pw(int ig)
        {
            return f_pw_(ig);
        }

        double value(vector3d<double>& vc)
        {
            int ja, jr;
            double dr, tp[2];
        
            if (unit_cell_.is_point_in_mt(vc, ja, jr, dr, tp)) 
            {
                int lmax = Utils::lmax_by_lmmax(angular_domain_size_);
                std::vector<double> rlm(angular_domain_size_);
                SHT::spherical_harmonics(lmax, tp[0], tp[1], &rlm[0]);
                double p = 0.0;
                for (int lm = 0; lm < angular_domain_size_; lm++)
                {
                    double d = (f_mt_(lm, jr + 1, ja) - f_mt_(lm, jr, ja)) / 
                               (unit_cell_.atom(ja)->type()->radial_grid(jr + 1) - unit_cell_.atom(ja)->type()->radial_grid(jr));
        
                    p += rlm[lm] * (f_mt_(lm, jr, ja) + d * dr);
                }
                return p;
            }
            else
            {
                double p = 0.0;
                for (int ig = 0; ig < num_gvec_; ig++)
                {
                    vector3d<double> vgc = gvec_.cart(ig);
                    p += std::real(f_pw_(ig) * std::exp(double_complex(0.0, vc * vgc)));
                }
                return p;
            }
        }

        int64_t hash()
        {
            STOP();

            int64_t h = Utils::hash(&f_it_(0), fft_->local_size() * sizeof(T));
            h += Utils::hash(&f_pw_(0), num_gvec_ * sizeof(double_complex), h);
            return h;
        }

        void fft_transform(int direction__)
        {
            switch (direction__)
            {
                case 1:
                {
                    fft_->input(gvec_.num_gvec_loc(), gvec_.index_map(), &f_pw(gvec_.gvec_offset()));
                    fft_->transform(1);
                    fft_->output(&f_it(0));
                    break;
                }
                case -1:
                {
                    fft_->input(&f_it(0));
                    fft_->transform(-1);
                    fft_->output(gvec_.num_gvec_loc(), gvec_.index_map(), &f_pw(gvec_.gvec_offset()));
                    fft_->comm().allgather(&f_pw(0), gvec_.gvec_offset(), gvec_.num_gvec_loc());
                    break;
                }
                default:
                {
                    TERMINATE("wrong fft direction");
                }
            }
        }
        
        mdarray<T, 3>& f_mt()
        {
            return f_mt_;
        }

        mdarray<T, 1>& f_it()
        {
            return f_it_;
        }

        mdarray<complex_t, 1>& f_pw()
        {
            return f_pw_;
        }

        static T inner(Periodic_function<T> const* f__, Periodic_function<T> const* g__)
        {
            assert(f__->fft_ == g__->fft_);
            assert(f__->step_function_ == g__->step_function_);
            assert(&f__->unit_cell_ == &g__->unit_cell_);
            assert(&f__->comm_ == &g__->comm_);
            
            //splindex<block> spl_fft_size(f__->fft_->size(), f__->comm_.size(), f__->comm_.rank());
        
            T result = 0.0;
            T ri = 0.0;
        
            if (f__->step_function_ == nullptr)
            {
                //for (int irloc = 0; irloc < (int)spl_fft_size.local_size(); irloc++)
                for (int irloc = 0; irloc < f__->fft_->local_size(); irloc++)
                    ri += type_wrapper<T>::conjugate(f__->f_it(irloc)) * g__->f_it(irloc);
            }
            else
            {
                //for (int irloc = 0; irloc < (int)spl_fft_size.local_size(); irloc++)
                for (int irloc = 0; irloc < f__->fft_->local_size(); irloc++)
                {
                    //int ir = (int)spl_fft_size[irloc];
                    ri += type_wrapper<T>::conjugate(f__->f_it(irloc)) * g__->f_it(irloc) * 
                          f__->step_function_->theta_r(irloc);
                }
            }
                    
            ri *= (f__->unit_cell_.omega() / f__->fft_->size());
            f__->fft_->comm().allreduce(&ri, 1);
            
            if (f__->parameters_.full_potential())
            {
                for (int ialoc = 0; ialoc < (int)f__->unit_cell_.spl_num_atoms().local_size(); ialoc++)
                    result += sirius::inner(f__->f_mt(ialoc), g__->f_mt(ialoc));
            }
        
            f__->comm_.allreduce(&result, 1);
        
            return result + ri;
        }
};

#include "periodic_function.hpp"

};

#endif // __PERIODIC_FUNCTION_H__

