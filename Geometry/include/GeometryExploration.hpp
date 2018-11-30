#pragma once
#include <Eigen/Dense>
#include <vector>
#include <set>
#include <iostream>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;
    template <typename T>
    using Vector2 = Eigen::Matrix<T, 2, 1>;
    template <typename T>
    using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    /*
    Returns v where
     - v > 0 if ccw
     - v < 0 if clockwise
     - v == 0 if collinear
    */
    template <typename T>
    T CCW(const Vector2<T>& p1, const Vector2<T>& p2, const Vector2<T>& p3) {
        return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
    }


    template <typename T>
    int lowestPointIndex(const std::vector<Vector2<T>>& points) {
        int i_best = 0;
        for (int i=0; i<points.size(); i++) {
            if (points[i][1] < points[i_best][1]) {
                i_best = i;
            }
        }
        return i_best;
    }



    // Q1: implement Graham Scan for 2D convex hull
    // The input is a vector of 2D points
    // The output should be the 2D points on the convex hull
    // Remove duplicates and sort your return vector
    // for comparing with the reference result
    template <typename T>
    std::vector<Vector2<T>> ConvexHull2D(const std::vector<Vector2<T>> &points) {
        int N = points.size();
        
        // make copy of points
        std::vector<Vector2<T>> ps = points;

        // Determine point with lowest y coordinate
        int i_best = lowestPointIndex(ps);
        
        // Swap with first element
        std::swap(ps[0], ps[i_best]);

        // Sort points by polar angle
        std::sort(
            ps.begin() + 1, 
            ps.end(), 
            [pts=ps](Vector2<T> p1, Vector2<T> p2) {
                return CCW(pts[0], p1, p2) > 0;
            }
        );        

        // Instantiate stack
        std::vector<Vector2<T>> stack { ps[0], ps[1], ps[2] };
       
        // Get convex hull 
        for (int i=3; i<ps.size(); i++) {
            while (CCW(
                    stack[stack.size()-2], 
                    stack[stack.size()-1], 
                    ps[i]) <= 0) {
                stack.pop_back();
            }
            stack.push_back(ps[i]);
        }

        // return convex hull (elements in stack)
        return stack;
    }

    // Q2: implement brute-force method for Nd Pareto front
    // The input is a vector of Nd points
    // The output should be the Nd points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<VectorX<T>> ParetoFrontNdNaive(const std::vector<VectorX<T>> &points) {
        return points;
    }

    // Q3: implement nlog(n) method for 2d Pareto front
    // The input is a vector of 2d points
    // The output should be the 2d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector2<T>> ParetoFront2D(const std::vector<Vector2<T>> &points) {
        return points;
    }

    // bonus question: implement nlog(n) method for 3d Pareto front
    // The input is a vector of 3d points
    // The output should be the 3d points on the Pareto front
    // sort your return vector for comparing with 
    // the reference result
    template <typename T>
    std::vector<Vector3<T>> ParetoFront3D(const std::vector<Vector3<T>> &points) {
        return points;
    }



}
