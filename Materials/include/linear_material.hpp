// Tao Du
// taodu@csail.mit.edu
// Sept 22, 2016
#pragma once
#include "material.hpp"

namespace materials {

    template<int dim, typename T>
    class LinearElasticityMaterial : public Material<dim, T> {
    public:
        LinearElasticityMaterial(const T young_modulus,
                                 const T poisson_ratio) : Material<dim, T>(young_modulus, poisson_ratio) {}

        LinearElasticityMaterial(const LinearElasticityMaterial<dim, T>& material)  : Material<dim, T>(material) {}

        ~LinearElasticityMaterial() {}


        const T EnergyDensity(const typename Material<dim, T>::MatrixDimT& F) const {
            const typename Material<dim, T>::MatrixDimT small_strain_tensor =
                    0.5 * (F + F.transpose()) - Material<dim, T>::MatrixDimT::Identity();
            const T trace = small_strain_tensor.trace();
            return Material<dim, T>::mu() * small_strain_tensor.array().square().sum()
                   + Material<dim, T>::lambda() / 2.0 * trace * trace;
        }



        const typename Material<dim, T>::MatrixDimT StressTensor(
                const typename Material<dim, T>::MatrixDimT& F) const {
            const typename Material<dim, T>::MatrixDimT I = Material<dim, T>::MatrixDimT::Identity();
            return Material<dim, T>::mu() * (F + F.transpose() - 2 * I) +
                   Material<dim, T>::lambda() * (F.trace() - dim) * I;
        }


        const typename Material<dim, T>::MatrixDimT StressDifferential(
                const typename Material<dim, T>::MatrixDimT& F,
                const typename Material<dim, T>::MatrixDimT& dF) const {

	    /*Note: This is wrong.  normally, this would return the differential dP/dF * dF.  
	      We have deleted the code (since it would give away the answer to the other StressDifferential).  
              You don't need this function to complete the assignment.
            */
            const typename Material<dim, T>::MatrixDimT I = Material<dim, T>::MatrixDimT::Identity();
            return I;

        };

        //TODO: Students fill this out (and change the return value)
        // F is the deformation gradient.
        /*
            - Number all elements of F (3x3 matrix) from 1-9 in the following way:
                [1 2 3]
                [4 5 6]
                [7 8 9]
            - Rewrite it in vector form [1 2 3 4 5 6 7 8 9]^T
            - Take derivative of P tensor accordingly. 
        */
        const typename Material<dim, T>::MatrixDim2T StressDifferential(
                const typename Material<dim, T>::MatrixDimT& F) const {

            const T mu = Material<dim, T>::mu();
            const T lambda = Material<dim, T>::lambda();
            const typename Material<dim, T>::MatrixDim2T I2 = Material<dim, T>::MatrixDim2T::Identity();
            
            // d/dF(mu*F)
            typename Material<dim, T>::MatrixDim2T stressDiff = mu * I2;

            // d/dF(mu*F^T)
            int pos = 0;
            for (int c=0; c<dim; c++) {
                for (int r=0; r<dim; r++) {
                    stressDiff(r*dim + c, pos) += mu;
                    pos++;
                }
            }
            
            // d/dF(lambda*Tr(F)*I)
            for (int r=0; r<dim; r++) {
                for (int c=0; c<dim; c++) {
                    stressDiff(r * (dim+1), c * (dim+1)) += lambda;
                }
            }

            return stressDiff;
        }

    private:
        LinearElasticityMaterial<dim, T>& operator=(
                const LinearElasticityMaterial<dim, T>&);



    };

}
