#pragma once

#include "read_stl.hpp"
#include "BasicGeometry.hpp"
#include "hexahedral_mesh.hpp"
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
    template <typename T>
    using Matrix3 = Eigen::Matrix<T, 3, 3>;

    /**
        This class is used to store counts in each voxel, which
        is helpful for voxelizing non-watertight meshes. After storing
        counts, you may call FilterThreshold() to convert to the
        standard boolean voxelization.
    */
    template<typename T>
    class VoxelPainter {
    public:
        /**
            Increment values of voxels. To convert it into a standard
            boolean voxel grid, we fill in any voxel with a count / value
            greater than the input threshold
        */
        VoxelPainter(T dx, Vector3<T> pmin, Eigen::Vector3i nvoxel) {
            _dx = dx;
            _pmin = pmin;
            _nvoxel = nvoxel;
            _voxels = std::vector<std::vector<std::vector<int>>>(_nvoxel.x(),
                std::vector<std::vector<int>>(_nvoxel.y(),
                    std::vector<int>(_nvoxel.z(), 0)));
        }

        Eigen::Vector3i ContainingVoxel(Vector3<T> point) {
            Vector3<T> relativePoint = (point - _pmin);
            Eigen::Vector3i gridPoint;
            for (int i=0; i<3; i++) {
                gridPoint[i] = static_cast<int>(relativePoint[i] / _dx);
            }
            return gridPoint;
        }

        void IncrementVoxel(Eigen::Vector3i vox) {
            for (int i=0; i<3; i++) {
                if (vox(i) >= _nvoxel(i) || vox(i) < 0) {
                    printf("this vox is out of bounds bro\n");
                }
            }
            _voxels[vox(0)][vox(1)][vox(2)] ++;
        }

        void FilterThreshold(
            T threshold, 
            std::vector<std::vector<std::vector<bool>>>& boolVoxels) {   // True <-> voxel is occupied.
            
            for (int i=0; i<_nvoxel[0]; ++i) {
                for (int j=0; j<_nvoxel[1]; ++j) {
                    for (int k=0; k<_nvoxel[2]; ++k) {
                        boolVoxels[i][j][k] = (_voxels[i][j][k] >= threshold);
                    }
                }
            }
        }

    private:
        T _dx;
        Vector3<T> _pmin;
        Vector3<int> _nvoxel; 
        std::vector<std::vector<std::vector<int>>> _voxels;   
    };


    /**
            This class is used for voxelization with non-axis aligned
        rays. That is why the constructor takes in a vector z.
        
            Note: When we refer to x'y'z' space, we are using the 
        non-axis aligned rays. When we refer to xyz space we are 
        referring to standard, not transformed space.
    */
    template<typename T>
    class RotatedProjector {
    public:
        // list of intersections for each x'y' grid location
        std::vector<std::vector<std::vector<Vector3<T>>>> Intersections;
        
        RotatedProjector(
            std::vector<std::vector<Vector3<T>>>& triangles, 
            T dx,
            Vector3<T>& z) {

            _triangles = triangles; 
            _dx = dx;
            _z = z.normalized();

            // Generate two other vectors _x, _y for our orthonormal axes
            Vector3<T> random = Vector3<T>::Random();
            _x = _z.cross(random).normalized();
            _y = _z.cross(_x).normalized();

            // Calculate change of basis matrix
            // TODO: Test this...
            for (int i=0; i<3; i++) {
                _COB_INV(i,0) = _x(i);
                _COB_INV(i,1) = _y(i);
                _COB_INV(i,2) = _z(i);
            }
            _COB = _COB_INV.inverse();

            // Compute the bounding box of _triangle and save the results into _pmin.
            T inf = (1 << 30);
            _pmin = Vector3<T>(inf, inf, inf);
            Vector3<T> pmax(-inf, -inf, -inf);
            for (const auto& triangle : _triangles)
                for (const auto& v : triangle) {
                    // change basis to x'y'z'
                    Vector3<T> projectedV = ProjectOntoAxes(v);
                    
                    // set min, max
                    _pmin = _pmin.cwiseMin(projectedV);
                    pmax = pmax.cwiseMax(projectedV);
                }
            for (int i = 0; i < 3; ++i) {
                _pmin[i] -= _dx;
                pmax[i] += _dx;
            }

            // Compute the number of voxels along each direction.
            for (int i = 0; i < 3; ++i)
                _nvoxel[i] = static_cast<int>((pmax[i] - _pmin[i]) / _dx) + 1;


            // Instantiate empty list of intersections
            for (int i = 0; i < _nvoxel[0]; ++i) {
                std::vector<std::vector<Vector3<T>>> emptyX;
                for (int j = 0; j < _nvoxel[1]; ++j) {
                    std::vector<Vector3<T>> emptyY;
                    emptyX.push_back(emptyY);
                }
                Intersections.push_back(emptyX);
            }
        }

        const Vector3<int> voxel_num() const { return _nvoxel; }

        /**
            Project a point in xyz-space into x'y'z' space
        */ 
        Vector3<T> ProjectOntoAxes(const Vector3<T>& point) {
            return _COB * point;
        }

        /**
            Project a point in x'y'z' space into xyz-space.
        */
        Vector3<T> ProjectOntoRealSpace(const Vector3<T>& point) {
            return _COB_INV * point;
        }

        /**
            Find the voxel grid coordinates (discretized) of a point
            given in x'y'z' space. It is IMPORTANT that you do not
            pass in a point in xyz space (aka point directly taken 
            from _triangles.)
        */
        Eigen::Vector3i GridConversion(Vector3<T>& point) {
            Vector3<T> relativePoint = (point - _pmin);
            Eigen::Vector3i gridPoint;
            for (int i=0; i<3; i++) {
                gridPoint[i] = static_cast<int>(relativePoint[i] / _dx);
            }
            return gridPoint;
        }

        /**
            Convert voxel grid coordinates (discretized) into a point
            given in x'y'z' space. To convert back to xyz space, use
            the ProjectOntoRealSpace() method.
        */
        Vector3<T> PointConversion(Eigen::Vector3i& grid) {
            Vector3<T> point;
            for (int i=0; i<3; i++) {
                point[i] = _pmin[i] + (grid[i] + 0.5) * _dx;
            }
            return point;
        }

        void RayWalk(int i, int j, VoxelPainter<T>& voxelpainter) {
            T step_size = (_dx / 4);    // TODO: Figure out better place to put this
            std::vector<Vector3<T>> intersections = Intersections[i][j];
            
            // Remove duplicate intersections (eg. when ray hits boundary between triangles)
            geometry::RemoveDuplicates<T>(intersections);

            // if the parity is odd, then we skip.
            if (intersections.size() % 2 == 1) {
                return;
            }
                   
            // Place z-coordinates into a new vector heights
            // Then sort them.
            std::vector<T> heights;
            geometry::GetComponents<T>(intersections, 2, heights);
            std::sort(heights.begin(), heights.end()); 


            // Determine corresponding x'y' coordinates
            T x = _pmin(0) + (i + 0.5) * _dx;
            T y = _pmin(1) + (j + 0.5) * _dx;
            
            // Loop over list of intersections, processing "in" regions
            // eg. between every other pair of intersections  
            for (int ind=1; ind<heights.size(); ind+=2) {
                T zmin = heights[ind-1];
                T zmax = heights[ind];
                
                // find the start and end in xyz-space
                Vector3<T> start = ProjectOntoRealSpace(Vector3<T>(x,y,zmin));
                Vector3<T> end = ProjectOntoRealSpace(Vector3<T>(x,y,zmax));
                Vector3<T> inc = step_size * (end - start).normalized();
                
                // Increment a point along this segment in real space
                // until point passes the end. If we hit the same voxel
                // as before, then SKIP!
                Vector3<T> point = start;
                Eigen::Vector3i lastVoxel(-1,-1,-1);
                while ( (end - point).dot(end - start) > 0 ) {
                    Eigen::Vector3i currVoxel = voxelpainter.ContainingVoxel(point);

                    // only process the current voxel if it is 
                    // different than the previous voxel
                    if (currVoxel != lastVoxel) {
                        voxelpainter.IncrementVoxel(currVoxel);
                        lastVoxel = currVoxel;
                    }
                    point += inc;
                }
                
            } 
        }


    private:
        Vector3<T> _x, _y, _z;
        Matrix3<T> _COB; // change of basis from standard xyz to x'y'z'
        Matrix3<T> _COB_INV;
        std::vector<std::vector<Vector3<T>>> _triangles;
        T _dx;  // The size of each voxel.
        Vector3<T> _pmin;    // The min and max corner of the bounding box.
        Eigen::Vector3i _nvoxel;   // The number of voxels along each direction.

    };


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

        const std::vector<std::vector<Vector3<T>>> triangles() const { return _triangles; }
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
            

            // loop over every triangle
            for (int ti=0; ti<_triangles.size(); ti++) {
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
                            _pmin[2] + _dx * 0.5                 
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
        
            // Store counts in each voxel to be processed later
            VoxelPainter<T> voxelpainter(_dx, _pmin, _nvoxel);
            
            // 11 is suggested number of iterations by the pset
            const int NUM_RAYS = 11;
            for (int k=0; k < NUM_RAYS; k++) {

                // generate a random vector direction
                Vector3<T> ray = Vector3<T>::Random().normalized();
                RotatedProjector<T> RP(_triangles, _dx, ray);

                // loop through all triangles
                for (auto iter=_triangles.begin(); iter!=_triangles.end(); ++iter) {
                    std::vector<Vector3<T>> verts = *iter;
    
                    // convert to x'y'z' space
                    geometry::Triangle<T> triangle(
                        RP.ProjectOntoAxes(verts[0]), 
                        RP.ProjectOntoAxes(verts[1]), 
                        RP.ProjectOntoAxes(verts[2]) 
                    );

                    // Find bounds in x'y'z' space               
                    T inf = (1 << 30);
                    Vector3<T> tmin(inf, inf, inf);
                    Vector3<T> tmax(-inf, -inf, -inf);
                    GetTriangleBounds(triangle, tmin, tmax);
                    Eigen::Vector3i gridmin = RP.GridConversion(tmin);
                    Eigen::Vector3i gridmax = RP.GridConversion(tmax);

                    // loop over projected area
                    for (int i = gridmin(0); i <= gridmax(0); ++i) {
                        for (int j = gridmin(1); j <= gridmax(1); ++j) {
                            // Determine voxel center in x'y'z' space
                            Eigen::Vector3i centerGrid(i, j, 0);
                            Vector3<T> center = RP.PointConversion(centerGrid);                
        
                            // Specify direction of ray in x'y'z' space
                            Vector3<T> dir(0,0,1);

                            // intersect triangle with ray
                            triangle.IntersectRay(center, dir, RP.Intersections[i][j]);
                        }
                    }
                }
                // Loop over all grid locations in the x'y' plane and
                // process rays in the z' direction
                Eigen::Vector3i nvoxelRotated = RP.voxel_num();
                for (int i=0; i < nvoxelRotated(0); i++) {
                    for (int j=0; j < nvoxelRotated(1); j++) {
                        RP.RayWalk(i, j, voxelpainter);
                    }
                }
            }
      
            // fill in all voxels with count above the threshold of
            // NUM_RAYS / 2 
            voxelpainter.FilterThreshold(NUM_RAYS/2.0, _voxels);
        }


        /*
            Convert mesh into a hex mesh
        */
        materials::HexahedralMesh<T> ConvertToHexMesh() {
            const int nx = _nvoxel[0], ny = _nvoxel[1], nz = _nvoxel[2];

            // locate elements
            int e[nx][ny][nz] = {0};

            int num_elements = 0;
            for (int i = 0; i < nx; ++i)
                for (int j = 0; j < ny; ++j)
                    for (int k = 0; k < nz; ++k)
                        if (_voxels[i][j][k])
                            e[i][j][k] = ++num_elements;
            Eigen::MatrixXi element(8, num_elements);

            // vertex array
            std::vector<Eigen::Vector3d> vertices;
            vertices.reserve(num_elements * 2);

            const int dx[8] = {0, 0, 0, 0, 1, 1, 1, 1};
            const int dy[8] = {0, 0, 1, 1, 0, 0, 1, 1};
            const int dz[8] = {0, 1, 0, 1, 0, 1, 0, 1};

            // fill elements matrix
            int vi = 0;
            for (int i = 0; i <= nx; ++i)
                for (int j = 0; j <= ny; ++j)
                    for (int k = 0; k <= nz; ++k) {
                        bool flag = false;
                        for (int d = 0; d < 8; d++) {
                            int i1 = i - dx[d], j1 = j - dy[d], k1 = k - dz[d];
                            if (i1 >= 0 && i1 < nx && j1 >= 0 && j1 < ny && k1 >= 0 && k1 < nz && _voxels[i1][j1][k1]) {
                                element(d, e[i1][j1][k1] - 1) = vi;
                                flag = true;
                            }
                        }
                        if (flag) {
                            vertices.push_back(_pmin + Eigen::Vector3d(i * _dx, j * _dx, k * _dx));
                            vi++;
                        }
                    }

            // create vertex matrix
            const int num_vertices = vertices.size();
            Eigen::MatrixXd vertex(3, num_vertices);
            for (int i = 0; i < num_vertices; ++i)
                vertex.col(i) = vertices[i];

            return std::move(materials::HexahedralMesh<T>(vertex, element));
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
