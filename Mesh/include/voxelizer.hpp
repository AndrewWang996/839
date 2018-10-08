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
 
            for (int i = 0; i < nx; ++i) {
                printf("Doing x voxel #%i out of %i\n", i, nx);
                for (int j = 0; j < ny; ++j) {
                    int k = 0;

                    // Determine voxel center
                    Vector3<T> center(
                        _pmin[0] + _dx * (i + 0.5),
                        _pmin[1] + _dx * (j + 0.5),
                        _pmin[2] + _dx * 0.5                  // it SHOULD be okay if it's slightly lower
                    );

                    // Specify direction of ray
                    Vector3<T> dir(0,0,1);

                    // Send an arbitrary ray outward, intersecting w/ all triangles in mesh
                    std::vector<Vector3<T>> intersections;
                    intersections.clear();
                    for (int ti=0; ti<_triangles.size(); ti++) {
                        std::vector<Vector3<T>> verts = _triangles[ti];

                        geometry::Triangle<T> triangle(
                            verts[0],
                            verts[1],
                            verts[2]
                        );

                        // gather intersections of ray with triangles in mesh
                        triangle.IntersectRay(center, dir, intersections);
                    }
                   
                    // set the corresponding voxels for these Z-axis aligned
                    // intersection points 
                    SetZVoxels(i, j, intersections);
                }
            }
        }

        void SetZVoxels(int i, int j, std::vector<Vector3<T>>& intersections) {
             // Remove duplicate intersections to account for strange cases
             // where intersections lie inbetween multiple triangles.
             geometry::RemoveDuplicates<T>(intersections);
 
             // Place z-coordinates into a new vector heights
             // Then sort them.
             std::vector<T> heights;
             geometry::GetComponents<T>(intersections, 2, heights);
             std::sort(heights.begin(), heights.end()); 

             if (intersections.size() % 2 == 1) {
                 printf("Num intersections is %lu. This is ODD!\n", intersections.size());
             }                      

             // Loop over list of intersections, only processing regions
             // where voxels should be set.
             for (int ind = 1; ind < heights.size(); ind += 2) {
                 int lo = ceil( (heights[ind-1] - _pmin[2]) / _dx );
                 int hi = floor( (heights[ind] - _pmin[2]) / _dx );

                 for (int zInd=lo; zInd<=hi; zInd++) {
                     _voxels[i][j][zInd] = true;
                 } 
             }
        }


        void GetTriangleBounds(
            geometry::Triangle<T>& triangle, 
            Vector3<T>& tmin, 
            Vector3<T>& tmax) {

            for (int i=0; i<3; i++) {
                Vector3<T> v = triangle.vertices(i);
                tmin = tmin.cwiseMin(v);
                tmax = tmax.cwiseMax(v);
            } 
        }

        void AdvancedVoxelization() {
            /* Assignment 2, Part 2.2. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
            
            // Store a list of intersection points for each grid location
            // on the xy-plane
            std::vector<std::vector<std::vector<Vector3<T>>>> intersections;
            for (int i = 0; i < nx; ++i) {
                std::vector<std::vector<Vector3<T>>> emptyX;
                for (int j = 0; j < ny; ++j) {
                    std::vector<Vector3<T>> emptyY;
                    emptyX.push_back(emptyY);
                }
                intersections.push_back(emptyX);
            }
            

            float inc = 0.1;

            // loop over every triangle
            for (int ti=0; ti<_triangles.size(); ti++) {
                // print out progress
                if ((float)ti / _triangles.size() > inc) {
                    printf("%i/%lu\n", ti, _triangles.size());
                    inc += 0.1;
                }

                // Instantiate triangle
                std::vector<Vector3<T>> verts = _triangles[ti];
                geometry::Triangle<T> triangle(verts[0], verts[1], verts[2]);

                // Get bounds of triangle to easily project it onto any
                // axis plane (eg. xy-plane)
                T inf = (1 << 30);
                Vector3<T> tmin(inf, inf, inf);
                Vector3<T> tmax(-inf, -inf, -inf);
                GetTriangleBounds(triangle, tmin, tmax);
                int xmin = (int)floor( (tmin[0] - _pmin[0]) / _dx - 0.5);
                int xmax = (int)ceil( (tmax[0] - _pmin[0]) / _dx + 0.5);
                int ymin = (int)floor( (tmin[1] - _pmin[1]) / _dx - 0.5);
                int ymax = (int)ceil( (tmax[1] - _pmin[1]) / _dx + 0.5);

                // loop over projected area
                for (int i = xmin; i <= xmax; ++i) {
                    for (int j = ymin; j <= ymax; ++j) {
                        // Determine voxel center
                        Vector3<T> center(
                            _pmin[0] + _dx * (i + 0.5),
                            _pmin[1] + _dx * (j + 0.5),
                            _pmin[2] + _dx * 0.5                  // it SHOULD be okay if it's slightly lower
                        );

                        // Specify direction of ray
                        Vector3<T> dir(0,0,1);

                        // intersect triangle with ray
                        triangle.IntersectRay(center, dir, intersections[i][j]);
                    }
                }
            }

            // set the voxels for each xy grid location
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    SetZVoxels(i, j, intersections[i][j]);    
                }
            }
        }

        void AdvancedVoxelizationWithApproximation() {
            /* Assignment 2, Part 2.3. */
            /* Implement your code here. */
            // Fill the _voxels array with the correct flag.
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];
        

            // 11 is suggested number of iterations by the pset
            for (int k=0; k<11; k++) {

                // generate a random vector direction
                Vector3<T> vec = Vector3<T>::Random().normalized();


                // loop through all triangles
                for (auto iter=_triangles.begin(); iter!=_triangles.end(); ++iter) {
                    std::vector<Vector3<T>> verts = *iter;
                    geometry::Triangle<T> triangle(verts[0], verts[1], verts[2]);

                    
                }


            }        



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
