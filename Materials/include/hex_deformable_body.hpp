// Tao Du
// taodu@csail.mit.edu
// Oct 12, 2016
#pragma once
#include "deformable_body.hpp"
#include "hexahedral_mesh.hpp"
#include "typedefs.hpp"
#include <stdio.h>

namespace materials {
    template<typename T>
    class HexDeformableBody : public DeformableBody<3, T> {
    public:
        HexDeformableBody(const Material<3, T>& material,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(material, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {}






        HexDeformableBody(const std::vector<std::reference_wrapper<const Material<3, T>>>& materials,
                          const std::vector<int>& material_id,
                          const Matrix3X<T>& initial_vertex_position,
                          const T density, const HexahedralMesh<T>& undeformed_hex_mesh)
                          : DeformableBody<3, T>(materials, material_id, initial_vertex_position, undeformed_hex_mesh),
                            hex_size_((undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(0)) - undeformed_hex_mesh.vertex(undeformed_hex_mesh.element(0)(1))).norm()) {}

        HexDeformableBody(const HexDeformableBody& hex_body) : DeformableBody<3, T>(hex_body),
                                                               hex_size_(hex_body.hex_size_) {}

        ~HexDeformableBody() {}


        //TODO: Students should fill this out
	    //vertices is a matrix of the current vertex positions (3 x n)        
        const Eigen::SparseMatrix<T> ComputeStiffnessMatrix(
                const Matrix3X<T>& vertices) const {

            std::vector<Eigen::Triplet<T>> triplet_list;
            const int vertex_num = static_cast<int>(this->vertex_position_.size() / 3);
           
            Eigen::SparseMatrix<T> K(vertex_num * 3, vertex_num * 3);

            auto &elements = this->undeformed_mesh_.element();
            int num_elements = elements.cols();
            triplet_list.reserve(24 * 24 * num_elements);   // reserve space

            T pos = 0.5 * (1 + 1 / std::sqrt(3));
            T neg = 0.5 * (1 - 1 / std::sqrt(3));
            Eigen::Matrix<T,3,8> quad;
            quad << neg,neg,neg,neg,pos,pos,pos,pos,
                     neg,neg,pos,pos,neg,neg,pos,pos,
                     neg,pos,neg,pos,neg,pos,neg,pos;

            Eigen::Matrix<T,3,3> phi;
            const Material<3,T> &material = this->materials_[0]; // let's just take the first material
            auto stress_diff = material.StressDifferential(phi);

            // local stiffness
            Eigen::Matrix<T, 24, 24> local_K;
            local_K.setZero();
            for (int i=0; i<8; i++) {
                auto dFdX = DeformationGradientPartialx(quad.col(i), 1.0 / hex_size_);
                auto tmp = dFdX.transpose() * stress_diff * dFdX;
                local_K += (tmp + tmp.transpose()) / 2.0;
            }
            local_K *= hex_size_ * hex_size_ * hex_size_ / 8.0;

            // global stiffness by summing local stiffness
            for (int e_i=0; e_i<num_elements; e_i++) {
                for (int i=0; i<8; i++) {
                    for (int j=0; j<8; j++) {
                        for (int x=0; x<3; x++) {
                            for (int y=0; y<3; y++) {
                                Eigen::Triplet<T> entry(
                                    3 * elements(j,e_i) + y, 
                                    3 * elements(i,e_i) + x,
                                    local_K(3*j + y, 3*i + x)
                                ); 
                                triplet_list.push_back(entry);
                            }
                        }
                    }
                }
            }
           
            K.setFromTriplets(triplet_list.begin(), triplet_list.end());
            // Make sure K is symmetric.
            K = (K + Eigen::SparseMatrix<T>(K.transpose())) / 2.0;
            return K;
        }


        //return dphi (the deformation gradient) for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 3, 3> DeformationGradient(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){
            // Rename variables.
            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            const Eigen::Matrix<T, 3, 8>& x0 = deformed_cube;
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = {1 - rx, rx};
            const T y_factor[2] = {1 - ry, ry};
            const T z_factor[2] = {1 - rz, rz};
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        F.col(0) += x0.col(4 * i + 2 * j + k)
                                    * (i == 0 ? -inv_dx : inv_dx) * y_factor[j] * z_factor[k];
                        F.col(1) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        F.col(2) += x0.col(4 * i + 2 * j + k)
                                    * x_factor[i] * y_factor[j] * (k == 0 ? -inv_dx : inv_dx);
                    }
            return F;
        }

        static const Eigen::Matrix<T, 9, 24> DeformationGradientPartialx(
                const Vector3<T> &r,
                const T inv_dx) {
            
            const T x_factor[2] = { 1 - r(0), r(0) };
            const T y_factor[2] = { 1 - r(1), r(1) };
            const T z_factor[2] = { 1 - r(2), r(2) };

            Eigen::Matrix<T,3,8> coeff;
            for (int i=0; i<=1; i++) {
                for (int j=0; j<=1; j++) {
                    for (int k=0; k<=1; k++) {
                        int ind = 4*i + 2*j + k;
                        coeff(0, ind) = (i == 0 ? -inv_dx : inv_dx) * y_factor[j] * z_factor[k];
                        coeff(1, ind) = x_factor[i] * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        coeff(2, ind) = x_factor[i] * y_factor[j] * (k == 0 ? -inv_dx : inv_dx);
                    }
                }
            }

            Eigen::Matrix<T,9,24> Jacobian;
            Jacobian.setZero();
            for (int i=0; i<8; i++) {
                for (int j=0; j<3; j++) {
                    for (int k=0; k<3; k++) {
                        Jacobian(j*3 + k, i*3 + k) = coeff(j,i);
                    }
                }
            }
            
            return Jacobian;
        }
        

        //return dphi/dx for a given voxel:
        //undeformed_vertex is a point in material space.
        //undeformed_cube are the vertices of an undeformed voxel
        //deformed_cube are the vertices of a the deformed voxel.
        static const Eigen::Matrix<T, 9, 24> DeformationGradientPartialx(
                const Vector3<T>& undeformed_vertex,
                const Eigen::Matrix<T, 3, 8>& undeformed_cube,
                const Eigen::Matrix<T, 3, 8>& deformed_cube){

            const Vector3<T>& X = undeformed_vertex;
            const Eigen::Matrix<T, 3, 8>& X0 = undeformed_cube;
            Eigen::Matrix<T, 9, 24> Jacobian = MatrixX<T>::Zero(9, 24);
            const T dx = X0(0, 4) - X0(0, 0);
            const T inv_dx = 1.0 / dx;
            const Vector3<T> offset = X - X0.col(0);
            const T rx = offset.x() / dx;
            const T ry = offset.y() / dx;
            const T rz = offset.z() / dx;
            Eigen::Matrix<T, 3, 3> F = Eigen::Matrix<T, 3, 3>::Zero();
            const T x_factor[2] = { 1 - rx, rx };
            const T y_factor[2] = { 1 - ry, ry };
            const T z_factor[2] = { 1 - rz, rz };
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < 2; ++k) {
                        const int index = 4 * i + 2 * j + k;
                        const T scale_first_column = (i == 0 ? -inv_dx : inv_dx)
                                                          * y_factor[j] * z_factor[k];
                        Jacobian(0, 3 * index) += scale_first_column;
                        Jacobian(1, 3 * index + 1) += scale_first_column;
                        Jacobian(2, 3 * index + 2) += scale_first_column;
                        const T scale_second_column = x_factor[i]
                                                           * (j == 0 ? -inv_dx : inv_dx) * z_factor[k];
                        Jacobian(3, 3 * index) += scale_second_column;
                        Jacobian(4, 3 * index + 1) += scale_second_column;
                        Jacobian(5, 3 * index + 2) += scale_second_column;
                        const T scale_third_column = x_factor[i] * y_factor[j]
                                                          * (k == 0 ? -inv_dx : inv_dx);
                        Jacobian(6, 3 * index) += scale_third_column;
                        Jacobian(7, 3 * index + 1) += scale_third_column;
                        Jacobian(8, 3 * index + 2) += scale_third_column;
                    }
            return Jacobian;



        }

        /*
            Permutation matrix to deal with arbitrary fixed vertices
        */
        Eigen::PermutationMatrix<Eigen::Dynamic> GetPermutation(const std::vector<int> &fixed) {
            int n_vert = this->undeformed_mesh_.NumOfVertex();
            std::vector<bool> ff(n_vert, false);
            for (int vi : fixed) {
                ff[vi] = true;
            }
            

            Eigen::VectorXi P_indices(3 * n_vert);
            int P_ii = 0;
            for (int i=0; i<n_vert; i++) {
                if (ff[i]) {
                    Eigen::Vector3i ind(3*i, 3*i+1, 3*i+2);
                    P_indices.segment<3>(P_ii) = ind;
                    P_ii += 3;
                }
            }
            for (int i=0; i<n_vert; i++) {
                if ( ! ff[i] ) {
                    Eigen::Vector3i ind(3*i, 3*i+1, 3*i+2);
                    P_indices.segment<3>(P_ii) = ind;
                    P_ii += 3;
                }
            }

            Eigen::PermutationMatrix<Eigen::Dynamic> P(P_indices);
            return P;
        }

        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;

        /*
            Get U given indices of fixed vertices + external forces
        */
        VectorXT GetDeformation(
                const std::vector<int> &fixed, 
                const VectorXT &F) {

            int n_vert = this->undeformed_mesh_.NumOfVertex();
            int n_fixed = fixed.size();
            int n_free = n_vert - n_fixed;
           
            Eigen::PermutationMatrix<Eigen::Dynamic> P = GetPermutation(fixed);

            const Eigen::SparseMatrix<T> K = ComputeStiffnessMatrix(
                this->undeformed_mesh_.vertex()
            );

            Eigen::SparseMatrix<T> K_P(3*n_vert, 3*n_vert);
            K_P = K.twistedBy(P.transpose()); // P^-1 K P

            Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Upper|Eigen::Lower> cg;
            cg.compute(
                K_P.bottomRightCorner(
                    3*n_free,
                    3*n_free
                )
            );

            VectorXT U(3*n_vert);
            U.setZero();
            U.tail(3*n_free) = cg.solve(
                (P.transpose() * F).tail(
                    3*n_free
                )
            );

            return P * U;
        }


        virtual void getInitialNodes(Matrix3X<T>& initial_nodes){
            initial_nodes = this->undeformed_mesh_.vertex();
        }

    private:
        HexDeformableBody& operator=(const HexDeformableBody&);


        //TODO: Studnets fill this function out
        static const Eigen::Matrix<T, 8, 8> GaussIntegrationFactor() {
            // \int_{-1}^{1} f(x) dx \approx f(-1/sqrt(3)) + f(1/sqrt(3)).
            Eigen::Matrix<T, 8, 8> X0_coeff = Eigen::MatrixXd::Zero(8, 8);
           
            // See Equation 11.2 from supplemental reading
            const T pos = 1 + 1 / std::sqrt(3);
            const T neg = 1 - 1 / std::sqrt(3);

            // N_{4(i-1)+2(j)+k} = (1+s*(-1)^i)(1+s*(-1)^j)(1+s*(-1)^k)
            const T N[8] = {
                neg*neg*neg/8,
                neg*neg*pos/8,
                neg*pos*neg/8,
                neg*pos*pos/8,
                pos*neg*neg/8,
                pos*neg*pos/8,
                pos*pos*neg/8,
                pos*pos*pos/8
            };

            // follow rules for 2-point quadrature
            // bitwise operations can be used for simplicity
            for (int i=0; i<8; i++) {
                for (int j=0; j<8; j++) {
                    X0_coeff(i,j) = N[i^j];
                }
            }
            return X0_coeff;
        }



        const T hex_size_;
    };

}
