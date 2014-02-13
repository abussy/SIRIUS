// This file is part of SIRIUS
//
// Copyright (c) 2013 Anton Kozhevnikov, Thomas Schulthess
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

#ifndef __DESCRIPTORS_H__
#define __DESCRIPTORS_H__

/** \file descriptors.h

    \brief Descriptors for various data structures
*/

/// describes single atomic level
struct atomic_level_descriptor
{
    /// principal quantum number
    int n;

    /// angular momentum quantum number
    int l;
    
    /// quantum number k
    int k;
    
    /// level occupancy
    double occupancy;

    /// true if this is a core level
    bool core;
};

/// describes radial solution
struct radial_solution_descriptor
{
    /// principal quantum number
    int n;
    
    /// angular momentum quantum number
    int l;
    
    /// order of energy derivative
    int dme;
    
    /// energy of the solution
    double enu;
    
    /// automatically determine energy
    int auto_enu;
};

/// set of radial solution descriptors, used to construct augmented waves or local orbitals
typedef std::vector<radial_solution_descriptor> radial_solution_descriptor_set;

/// descriptor of a local orbital radial function
struct local_orbital_descriptor
{
    local_orbital_t type;

    int l;

    radial_solution_descriptor_set rsd_set;

    int p1;
    int p2;
};


class uspp_descriptor
{
    public:
        std::vector<double> r;
        std::vector<double> vloc;

        /// maximum angular momentum for |beta> projectors
        int lmax;

        /// number of radial functions for |beta> projectors
        int num_beta_radial_functions;
        
        /// number of Q coefficients
        int num_q_coefs; 
        
        /// Q coefficients
        mdarray<double, 4> q_coefs;

        std::vector<double> q_functions_inner_radii;

        mdarray<double, 2> q_radial_functions;

        std::vector<int> beta_l;

        std::vector<int> num_beta_radial_points;

        mdarray<double, 2> beta_radial_functions;

        std::vector<double> core_charge_density;

        std::vector<double> total_charge_density;

        mdarray<double_complex, 2> d_mtrx_ion;

        mdarray<double_complex, 2> d_mtrx;

        mdarray<double_complex, 2> q_mtrx;

        mdarray<double_complex, 2> q_pw;
};

struct nearest_neighbour_descriptor
{
    /// id of neighbour atom
    int atom_id;

    /// translation along each lattice vector
    vector3d<int> translation;

    /// distance from the central atom
    double distance;
};

#endif
