#pragma once

#include "read_stl.hpp"
#include "BasicGeometry.hpp"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <unordered_set>
#include <fstream>
#include <ctime>

namespace mesh {

    template<typename T>
    class Voxelizer {
    public:
        Voxelizer(const std::string& stl_file_name, const T dx)
            : _dx(dx) {
            // Randomness.
            srand(static_cast<unsigned>(time(0)));
            // Load triangles from the stl file.
            std::vector<Vector3<T>> normals;
            if (!ReadSTL(stl_file_name, _triangles, normals)) {
                std::cout << "ERROR: cannot read " << stl_file_name << std::endl;
                return;
            }
            // Compute the bounding box of _triangle and save the results into _pmin.
            _pmin = _triangles[0][0];
            Vector3<T> pmax = _triangles[0][0];
            for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                    _pmin = _pmin.cwiseMin(v);
                    pmax = pmax.cwiseMax(v);
                }
            for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
            }
            // Compute the number of voxels along each direction.
            for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;
            // Initialize the voxel array.
            _voxels = std::vector<std::vector<std::vector<bool>>>(_nvoxel.x(),
                std::vector<std::vector<bool>>(_nvoxel.y(),
                    std::vector<bool>(_nvoxel.z(), false)));
        }

        const Vector3<T> pmin() const { return _pmin; }
        const T dx() const { return _dx; }
        const Vector3<T> pmax() const { return _pmin + Vector3<T>(_nvoxel.x(), _nvoxel.y(), _nvoxel.z()) * _dx; }
        const Vector3<int> voxel_num() const { return _nvoxel; }

        void BasicVoxelization() {
            /* Assignment 2, Part 2.1. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            
            std::vector<Vector3<T>> intersections;
            for (int i = 0; i < nx; ++i) {
                printf("Doing x voxel #%i out of %i\n", i, nx);
                for (int j = 0; j < ny; ++j) {
                    printf("Doing y voxel #%i out of %i\n", j, nx);
                    for (int k = 0; k < nz; ++k) {
                        printf("Doing z voxel #%i out of %i\n", k, nx);
                        intersections.clear();

                        // Determine voxel center
                        Vector3<T> center(
                            _pmin[0] + _dx * (i + 0.5),
                            _pmin[1] + _dx * (j + 0.5),
                            _pmin[2] + _dx * (k + 0.5)
                        );

                        // Specify direction of ray
                        Vector3<T> dir(1,1,1);

                        // Send an arbitrary ray outward, intersecting w/ all triangles in mesh
                        for (auto iter=_triangles.begin(); iter!=_triangles.end(); ++iter) {
                            std::vector<Vector3<T>> vertices = *iter;
                            geometry::Triangle<T> triangle(
                                vertices[0],
                                vertices[1],
                                vertices[2]
                            );

                            // gather intersections of ray with triangles in mesh
                            triangle.IntersectRay(center, dir, intersections);
                        }
                        
                        // Remove duplicate intersections to account for strange cases
                        // where intersections lie inbetween multiple triangles.
                        geometry::RemoveDuplicates<T>(intersections);                        

                        // voxel is filled iff number of intersections is odd
                        _voxels[i][j][k] = (intersections.size() % 2 == 1);
                    }
                }
            }
        }

        void AdvancedVoxelization() {
            /* Assignment 2, Part 2.2. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = false;
        }

        void AdvancedVoxelizationWithApproximation() {
            /* Assignment 2, Part 2.3. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        _voxels[i][j][k] = false;
        }

        void WriteVoxelToMesh(const std::string& stl_file_name) const {
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            std::vector<std::vector<Vector3<int>>> faces;
            std::vector<Vector3<int>> corners({
                Vector3<int>(0, 0, 0),
                Vector3<int>(0, 0, 1),
                Vector3<int>(0, 1, 0),
                Vector3<int>(0, 1, 1),
                Vector3<int>(1, 0, 0),
                Vector3<int>(1, 0, 1),
                Vector3<int>(1, 1, 0),
                Vector3<int>(1, 1, 1)
            });
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k) {
                        if (!_voxels[i][j][k]) continue;
                        // Check -x direction.
                        Vector3<int> cmin(i, j, k);
                        if (i == 0 || !_voxels[i - 1][j][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[1], cmin + corners[3] });
                            faces.push_back({ cmin + corners[0], cmin + corners[3], cmin + corners[2] });
                        }
                        if (i == nx - 1 || !_voxels[i + 1][j][k]) {
                            faces.push_back({ cmin + corners[4], cmin + corners[6], cmin + corners[7] });
                            faces.push_back({ cmin + corners[4], cmin + corners[7], cmin + corners[5] });
                        }
                        if (j == 0 || !_voxels[i][j - 1][k]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[4], cmin + corners[5] });
                            faces.push_back({ cmin + corners[0], cmin + corners[5], cmin + corners[1] });
                        }
                        if (j == ny - 1 || !_voxels[i][j + 1][k]) {
                            faces.push_back({ cmin + corners[2], cmin + corners[3], cmin + corners[7] });
                            faces.push_back({ cmin + corners[2], cmin + corners[7], cmin + corners[6] });
                        }
                        if (k == 0 || !_voxels[i][j][k - 1]) {
                            faces.push_back({ cmin + corners[0], cmin + corners[2], cmin + corners[6] });
                            faces.push_back({ cmin + corners[0], cmin + corners[6], cmin + corners[4] });
                        }
                        if (k == nz - 1 || !_voxels[i][j][k + 1]) {
                            faces.push_back({ cmin + corners[5], cmin + corners[7], cmin + corners[3] });
                            faces.push_back({ cmin + corners[5], cmin + corners[3], cmin + corners[1] });
                        }
                    }
            std::ofstream fout(stl_file_name);
            fout << "solid vcg" << std::endl;
            for (const auto& f : faces) {
                std::vector<Vector3<T>> p;
                for (const auto& fi : f) {
                    Vector3<T> v = _pmin + fi.cast<T>() * _dx;
                    p.push_back(v);
                }
                const Vector3<T> n = (p[1] - p[0]).cross(p[2] - p[1]).normalized();
                fout << "  facet normal " << n.x() << " " << n.y() << " " << n.z() << std::endl;
                fout << "    outer loop" << std::endl;
                for (const auto& v : p) {
                    fout << "      vertex " << v.x() << " " << v.y() << " " << v.z() << std::endl;
                }
                fout << "    endloop" << std::endl;
                fout << "  endfacet" << std::endl;
            }
            fout << "endsolid vcg" << std::endl;
        }

    private:
        std::vector<std::vector<Vector3<T>>> _triangles;
        T _dx;  // The size of each voxel.
        Vector3<T> _pmin;    // The min and max corner of the bounding box.
        Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.
        std::vector<std::vector<std::vector<bool>>> _voxels;   // True <-> voxel is occupied.
    };

}
