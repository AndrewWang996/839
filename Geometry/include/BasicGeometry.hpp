#pragma once
#include <Eigen/Dense>
#include <vector>

namespace geometry {
    template <typename T>
    using Vector3 = Eigen::Matrix<T, 3, 1>;

    float EPS = 1e-6;
    
    // namespace functions
    template <typename T>
    void RemoveDuplicates(std::vector<Vector3<T>>& points) {
        for (int i=points.size() - 1; i>=0; i--) {
            for (int j=0; j<i; j++) {
                if ( (points[i] - points[j]).norm() < EPS ) {
                    points.erase(points.begin() + i);
                    break;
                }
            }
        }
    }

    // the plane is represented by (x - _p) /dot _normal = 0
    template <typename T>
    class Plane {
    public:
        Plane(Vector3<T> p, Vector3<T> normal) {
            _p = p;
            _normal = normal;
            _normal.normalize();
        }

        Vector3<T>& p() { return _p; }
        Vector3<T>& normal() { return _normal; }
        
        // return if the point is on plane
        // also fill parameter dist field as the signed distance from point to plane
        bool onPlane(Vector3<T> point, T& dist) {
            dist = (point - _p).dot(_normal);
            if (std::fabs(dist) < 1e-6) {
                return true;
            } else {
                return false;
            }
        }

    private:
        Vector3<T> _p;
        Vector3<T> _normal;
    };




    template <typename T>
    class LineSegment {
    public:
        LineSegment(Vector3<T> v0, Vector3<T> v1) {
            _vertices[0] = v0;
            _vertices[1] = v1;
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }


        /**
            We're literally following word-for-word what the Wikipedia page for 
            line plane intersections says.

            https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
        */
        void IntersectPlane(Plane<T> p, std::vector<Vector3<T>>& intersections) {
            Vector3<T> n = p.normal();
            Vector3<T> p0 = p.p();
            Vector3<T> l0 = _vertices[0];
            Vector3<T> l = _vertices[1] - _vertices[0];
            
            // if line and plane are parallel, then no intersections
            if ( abs(l.dot(n)) < EPS) {
                return;
            }

            // if intersection does not land on segment, skip
            float d = (p0 - l0).dot(n) / l.dot(n);
            if (d < -EPS || d > 1+EPS) {
                return;
            }

            // add poi (point of intersection) to list
            Vector3<T> poi = d * l + l0;
            intersections.push_back(poi);
        }

    private:
        float EPS = 1e-6;

        Vector3<T> _vertices[2];
    }; 




    template <typename T>
    class Triangle {
    public:
        Triangle(Vector3<T> v0, Vector3<T> v1, Vector3<T> v2) {
            _vertices[0] = v0;
            _vertices[1] = v1;
            _vertices[2] = v2;
            _normal = (v2 - v0).cross(v1 - v0).normalized(); 
            _area = (v2 - v0).cross(v1 - v0).norm() / 2.0;
        }

        Vector3<T>* vertices() { return _vertices; }
        Vector3<T>& vertices(int idx) { return _vertices[idx]; }


        // check sums of areas of 3 components == area of triangle
        bool ContainsPoint(Vector3<T>& point) const {
            Vector3<T> subarea1 = (point - _vertices[0]).cross(point - _vertices[2]);
            Vector3<T> subarea2 = (point - _vertices[2]).cross(point - _vertices[1]);
            Vector3<T> subarea3 = (point - _vertices[1]).cross(point - _vertices[0]);
            Vector3<T> area = (_vertices[2] - _vertices[0]).cross(_vertices[1] - _vertices[0]);
            
            T diffNorm = ((subarea1 + subarea2 + subarea3) - area).norm();
            return diffNorm < EPS;
        }

        // Assignment 2: implement ray-triangle intersection.
        // The ray is defined as r(t) = origin + t * dir.
        // You should return a scalar t such that r(t) is the intersection point. Which value
        // to return for the case of no intersections is up to you. You can also change the
        // signature of this function as you see fit.
        // const T IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir) const {
        void IntersectRay(const Vector3<T>& origin, const Vector3<T>& dir, std::vector<Vector3<T>>& intersections) const {
            /* Assignment 2. */
            Plane<T> plane(_vertices[0], _normal);

            Vector3<T> n = _normal;
            Vector3<T> p0 = _vertices[0];
            Vector3<T> l0 = origin;
            Vector3<T> l = dir;
            
            // if line and plane are parallel, then no intersections
            if ( abs(l.dot(n)) < EPS) {
                return;
            }

            // get poi (point of intersection)
            float d = (p0 - l0).dot(n) / l.dot(n);
            Vector3<T> poi = d * l + l0;

            // test to see if poi is on ray
            if (d < 0) {
                return;
            }

            // test to see if intersection is contained in the triangle
            if (ContainsPoint(poi)) {
                intersections.push_back(poi);
            }
        }
        
        /* Implement triangle plane intersection.
            Input is a plane, output is a list of intersection points
            Hints:
                1. Take some considerations: handle duplicate intersections (when plane intersect with 
                    triangle at corner but you find them by edge-plane intersection).
                2. You can modify the function input/return type to any format you want.
        */
        std::vector<Vector3<T>> IntersectPlane(Plane<T> p) {
            std::vector<Vector3<T>> intersections;            

            LineSegment<T> seg1(_vertices[0], _vertices[1]);
            LineSegment<T> seg2(_vertices[0], _vertices[2]);
            LineSegment<T> seg3(_vertices[1], _vertices[2]);
           
            seg1.IntersectPlane(p, intersections);
            seg2.IntersectPlane(p, intersections);
            seg3.IntersectPlane(p, intersections);

            RemoveDuplicates<T>(intersections);
            return intersections;
        }


    private:
        float EPS = 1e-6;
        Vector3<T> _vertices[3];
        Vector3<T> _normal;
        T _area;
    };



}
