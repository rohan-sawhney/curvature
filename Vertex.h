#ifndef VERTEX_H
#define VERTEX_H

#include "Types.h"

class Vertex {
public:
    // outgoing halfedge
    HalfEdgeIter he;
    
    // location in 3d
    Eigen::Vector3d position;
    
    // id between 0 and |V|-1
    int index;
       
    // gaussian curvature, normalized between -1 and 1
    double gaussCurvature;
    
    // mean curvature, normalized between -1 and 1
    double meanCurvature;
    
    // checks if vertex is contained in any edge or face
    bool isIsolated() const;
    
    // returns area of barycentric dual cell associated with the vertex
    double dualArea() const;
    
    // 2pi - ∑ø
    double angleDefect() const;
};

#endif
