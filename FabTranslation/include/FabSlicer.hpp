#pragma once
#include "tri_mesh.hpp"
#include "BasicGeometry.hpp"
#include "IntervalTree.hpp"
#include "cinolib/meshes/meshes.h"
#include <ctime>
#include <math.h>
#include <algorithm>

namespace fab_translation {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    template <typename T>
    class FabSlicer {
        
    public:
        FabSlicer(mesh::TriMesh<T> tri_mesh, T bottom, T top, T dx,
            T infill_dx)
            : _tri_mesh(tri_mesh), _bottom(bottom), _top(top), _dx(dx), _infill_dx(infill_dx) {

            /* Implement your code here */
            /* 1. Initialize your variables
               2. Build interval tree */
        }

        /* Main entrance for FabSlicer
            return contour and infill_edges, which can be directed sent to corresponding 
                visualization function to visualize in MeshLab
            1. each contour contains a list of point in loop order, each layer can have 
                multiple contours, so the contour is vector<vector<vector<Point>>>.
            2. infill edges is a edge soup for each layer, thus it's vector<vector<Edge>>
        */
        void RunTranslation(std::vector<std::vector<std::vector<Vector3<T>>>>& contour,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {

            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> intersection_edges;

            // Slicing_bruteforce(_tri_mesh, intersection_edges);
            Slicing_accelerated(_tri_mesh, intersection_edges);

            CreateContour(_tri_mesh, intersection_edges, contour);

            Infill(contour, infill_edges);
        }

        /* Slicing algorithms
            goal: slice the triangle mesh by a set of parallel planes, 
                  output an intersection edge soup for each layer */
        void Slicing_bruteforce(mesh::TriMesh<T>& tri_mesh, 
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges) {

            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();

            for (T h = _bottom; h <= _top; h += _dx) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> intersections_one_plane;
                intersections_one_plane.clear();

                geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));
                for (int i = 0;i < elements.size();++i) {
                    geometry::Triangle<T> triangle(vertices[elements[i](0)], vertices[elements[i](1)], vertices[elements[i](2)]);
                    std::vector<Vector3<T>> intersections = triangle.IntersectPlane(plane);

                    /* Implement your code here */
                    if (intersections.size() == 2) {
                        std::pair<Vector3<T>, Vector3<T>> seg(intersections[0], intersections[1]);
                        intersections_one_plane.push_back(seg);
                    }
                }

                intersection_edges.push_back(intersections_one_plane);
            }
        }


        /* converts a z value into a bucket by dividing by the z increment size */
        int To_Z_Bucket(T z) {
            return (int)round(z / _dx);
        }


        /*
            Preprocessing method to get the z buffer in the form of 
            map[ bucket ][ list of triangles ]
        */
        std::map<int, std::vector<geometry::Triangle<T>>> Get_Z_Buffer(mesh::TriMesh<T>& tri_mesh) {
            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();

            std::map<int, std::vector<geometry::Triangle<T>>> zBuff;

            // initialize the relevant buckets as empty vectors
            for (T h = _bottom; h <= _top; h += _dx) {
                int bucket = To_Z_Bucket(h);
                std::vector<geometry::Triangle<T>> empty;
                zBuff[bucket] = empty;
            }
 
            for (auto iter=elements.begin(); iter!=elements.end(); ++iter) {
                Eigen::Vector3i elem = *iter;
                Vector3<T> A = vertices[elem(0)];
                Vector3<T> B = vertices[elem(1)];
                Vector3<T> C = vertices[elem(2)];

                geometry::Triangle<T> triangle(A,B,C);

                // determine unique set of buckets to insert triangle into
                int bucketA = To_Z_Bucket(A[2]);
                int bucketB = To_Z_Bucket(B[2]);
                int bucketC = To_Z_Bucket(C[2]);

                int minBucket = std::min(std::min(bucketA, bucketB), bucketC);
                int maxBucket = std::max(std::max(bucketA, bucketB), bucketC);

                // insert into all buckets that contain part of the triangle
                for (int bucket = minBucket; bucket <= maxBucket; bucket++) {
                    zBuff[bucket].push_back(triangle);
                }
            
            }

            return zBuff;
        }

        void Slicing_accelerated(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges) {
            
            std::vector<Eigen::Vector3i>& elements = tri_mesh.elements();
            std::vector<Vector3<T>>& vertices = tri_mesh.vertices();
            std::vector<Eigen::Vector3i>& edges = tri_mesh.edges();

            intersection_edges.clear();

            // Preprocess by storing the triangles in elements into Z-buckets
            std::map<int, std::vector<geometry::Triangle<T>>> zBuff = Get_Z_Buffer(tri_mesh);

            for (T h = _bottom; h <= _top; h += _dx) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> intersections_one_plane;
                intersections_one_plane.clear();

                geometry::Plane<T> plane(Vector3<T>(0, 0, h), Vector3<T>(0, 0, 1));

                // convert the current height h to a z bucket.
                int bucket = To_Z_Bucket(h); 

                // grab the triangles that fall into this bucket.
                std::vector<geometry::Triangle<T>> candidate_triangles = zBuff[bucket]; 

                // loop over all triangles, checking to see which ones intersect the plane                
                for (auto iter=candidate_triangles.begin(); iter!=candidate_triangles.end(); ++iter) {
                    geometry::Triangle<T> triangle = *iter;
                    std::vector<Vector3<T>> intersections = triangle.IntersectPlane(plane);
                    
                    // add any triangles with 2 intersections with the plane
                    if (intersections.size() == 2) {
                        std::pair<Vector3<T>, Vector3<T>> seg(intersections[0], intersections[1]);
                        intersections_one_plane.push_back(seg);
                    }
                }

                intersection_edges.push_back(intersections_one_plane);
            }
        }

        /* Find contours
            Goal: Given an intersetion edge soup for each layer, link those edges one by one to form contour loops.
                  Each layer probably has several disjoint contour loops. */
        void CreateContour(mesh::TriMesh<T>& tri_mesh,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> &intersection_edges,
            std::vector<std::vector<std::vector<Vector3<T>>>>& contours) {
            
            /* Implement your code here */
            /* Input is a edge soup, your task is to generate the loop by linking those edges one by one.
               Thinking about how to find two edge sharing a same end point. set a threshold? or a more clever way? */
            
            // Iterate over every contour's list of intersectione edges
            for (auto contour_inter_iter = intersection_edges.begin(); 
                contour_inter_iter != intersection_edges.end(); 
                contour_inter_iter++) {
                
                // a single z-plane / contour's intersection edges
                std::vector<std::pair<Vector3<T>, Vector3<T>>> contour_inter_edges = *contour_inter_iter;

                // a single z-plane's [possibly multiple] loop(s)
                std::vector<std::vector<Vector3<T>>> contour_loops;

                // add the loops to the above vector
                CreateLoops(contour_inter_edges, contour_loops);

                // add the z-plane's loop(s) to contours
                contours.push_back(contour_loops);
            }


        }


        /**
            Given a single z-plane's list of intersection edges, there might be multiple 
            contour loops (imagine the sequence of pillars in Killian court, each of which
            would produce its own separate contour loop.) Add these contour loops to the
            contour_loops passed in as a parameter.
        */
        void CreateLoops(
            std::vector<std::pair<Vector3<T>, Vector3<T>>>& intersection_edges,
            std::vector<std::vector<Vector3<T>>>& contour_loops) {

            // mark edges visited when we get to them...            
            std::vector<bool> visited(intersection_edges.size());

            // loop through all unvisited edges
            for (int i=0; i<intersection_edges.size(); i++) {
                if (visited[i]) {       
                    continue;
                }
                
                // create a single loop
                std::vector<Vector3<T>> single_loop;
                CreateSingleLoop(intersection_edges, visited, i, single_loop);

                // add single loop to all contour_loops
                contour_loops.push_back(single_loop);
            }
            
        }


        void CreateSingleLoop(
            std::vector<std::pair<Vector3<T>, Vector3<T>>>& intersection_edges,
            std::vector<bool>& visited,
            int idx,
            std::vector<Vector3<T>>& loop) {

            // indicator variable to say whether our current vertex is the 
            // first or second element of our current edge (a C++ pair)
            bool first = true;

            while ( ! visited[idx] ) {
                visited[idx] = true;        // mark edge as visited
            
                // set current and next vertices by looking at current edge
                std::pair<Vector3<T>, Vector3<T>> edge = intersection_edges[idx];
                Vector3<T> v_curr = (first ? edge.first : edge.second);
                Vector3<T> v_next = (first ? edge.second : edge.first);

                // add current vertex to loop
                loop.push_back(v_curr);

                // look for unvisited edge that has "v_next" as a endpoint
                idx = GetUnvisitedEdgeWithEndpoint(intersection_edges, idx, !first);
                
                // check to see if our current index is now the first or second
                // element of our NEW current edge
                if ( (intersection_edges[idx].first - v_next).norm() < EPS ) { first = true; }
                else { first = false; }
            }

            // complete the loop
            if (loop.size() > 0) {
                // get last edge
                std::pair<Vector3<T>, Vector3<T>> lastEdge = intersection_edges[idx];

                // add the endpoint that was not the last vertex
                Vector3<T> lastPoint = (first ? lastEdge.first : lastEdge.second);
                loop.push_back(lastPoint);
            }    
        }

        /**
            Searches through all unvisited edges to check if any have the
            given endpoint. The check is done by seeing if the L2 norm is
            less than some epsilon value
        */
        int GetUnvisitedEdgeWithEndpoint(
            std::vector<std::pair<Vector3<T>, Vector3<T>>>& intersection_edges,
            int idx,        // index of the current edge 
            bool first) {   // which endpoint of the current edge
               
            std::pair<Vector3<T>, Vector3<T>> edge = intersection_edges[idx];
            Vector3<T> endpoint = (first ? edge.first : edge.second);
 
            for (int i=0; i<intersection_edges.size(); i++) {
                if (i == idx) {
                    continue;
                }
                Vector3<T> vA = intersection_edges[i].first;
                Vector3<T> vB = intersection_edges[i].second;
               
                if ((vA - endpoint).norm() < EPS || (vB - endpoint).norm() < EPS) {
                    return i;
                }
            }
           

            printf("WHAT THE FUCKKKK!!!\n");
            // no unvisited edges with the given endpoint found...
            // actually this function should never reach here.
            return -1;
        }
 

        /* Generate infill pattern
           Goal: Given the contours at each layer, this function aims to infill the internal part by a pattern which
                 can be procedurally built. (e.g. grid, honey comb, Fermat spiral) 
           The code is for grid pattern, you can rewrite the whole function based on your need. */
        void Infill(std::vector<std::vector<std::vector<Vector3<T>>>>& contours,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            
            infill_edges.clear();
            

            for (int i = 0;i < contours.size();++i) {
                std::vector<std::pair<Vector3<T>, Vector3<T>>> infill_edges_one_layer;
                infill_edges_one_layer.clear();

                /* Implement your code here */
                /* 1. find all intersections between contours and your infill pattern 
                2. infill internal space with desired pattern */

                infill_edges.push_back(infill_edges_one_layer);
            }
        }  

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            intersection_edges: edge soup for each layer's slicing */
        void VisualizeSlicing(std::string file_name, 
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>> intersection_edges) {
            
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < intersection_edges.size();++i)
                for (int j = 0;j < intersection_edges[i].size();++j) {
                    Vector3<T> s_pos = intersection_edges[i][j].first;
                    Vector3<T> t_pos = intersection_edges[i][j].second;
                    int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                    for (int step = 0;step <= num_steps;++step) {
                        Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                        points.push_back(pos);
                    }
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's contour
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            contour: each layer's contour list, each contour is a list a point in loop order */
        void VisualizeContour(std::string file_name,
            T point_density, 
            std::vector<std::vector<std::vector<Vector3<T>>>>& contour) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < contour.size();++i)
                for (int j = 0;j < contour[i].size();++j) 
                    for (int k = 0;k < contour[i][j].size();++k) {
                        Vector3<T> s_pos = contour[i][j][k];
                        Vector3<T> t_pos = contour[i][j][(k + 1) % contour[i][j].size()];
                        int num_steps = (int)((t_pos - s_pos).norm() / point_density) + 1;
                        for (int step = 0;step <= num_steps;++step) {
                            Vector3<T> pos = s_pos * ((T)step / num_steps) + t_pos * ((T)1.0 - (T)step / num_steps);
                            points.push_back(pos);
                        }
                    }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

        /* Point cloud visualization for each layer's slicing
            file_name: output .ply filename (should be with .ply extension) 
            point_density: the smaller, the denser 
            infill_edges: edge soup for each layer's slicing */
        void VisualizeInfill(std::string file_name,
            T point_density,
            std::vector<std::vector<std::pair<Vector3<T>, Vector3<T>>>>& infill_edges) {
            // generate point cloud for ply
            std::vector<Vector3<T>> points;
            points.clear();
            for (int i = 0;i < infill_edges.size();++i)
                for (int j = 0;j < infill_edges[i].size();++j) {
                    int num_steps = (int)((infill_edges[i][j].first - infill_edges[i][j].second).norm() / point_density) + 1;
                    for (int k = 0;k <= num_steps;++k)
                        points.push_back(infill_edges[i][j].first + (infill_edges[i][j].second - infill_edges[i][j].first) * (T)k / (T)num_steps);
                }

            // output to ply
            FILE* fp = fopen(file_name.c_str(), "w");
            fprintf(fp, "ply\nformat ascii 1.0\n");
            fprintf(fp, "element vertex %d\n", (int)points.size());
            fprintf(fp, "property float32 x\nproperty float32 y\nproperty float32 z\n");
            fprintf(fp, "end_header\n");
            for (int i = 0;i < points.size();++i)
                if (std::is_same<T, float>::value)
                    fprintf(fp, "%.6f %.6f %.6f\n", points[i](0), points[i](1), points[i](2));
                else
                    fprintf(fp, "%.6lf %.6lf %.6lf\n", points[i](0), points[i](1), points[i](2));
            fclose(fp);
        }

    private:

        float EPS = 1e-6;

        mesh::TriMesh<T> _tri_mesh;

        /* Variables for slicing */
        T _bottom, _top, _dx;

        /* Variables for infill algorithm */
        T _infill_dx;                                   // infill pattern will be equal-length grid
        T _infill_x_lower_bound, _infill_x_upper_bound;
        T _infill_y_lower_bound, _infill_y_upper_bound;

        /* accelerated data structure */
        data_structure::IntervalTree<T> _interval_tree;
    };
}
