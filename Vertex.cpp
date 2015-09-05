#include "Vertex.h"
#include "HalfEdge.h"
#include "Face.h"

std::vector<HalfEdge> isolated;

bool Vertex::isIsolated() const
{
    return he == isolated.begin();
}

double Vertex::dualArea() const
{
    double area = 0.0;
    
    HalfEdgeCIter h = he;
    do {
        area += h->face->area();
        h = h->flip->next;
        
    } while (h != he);
    
    return area / 3.0;
}

double Vertex::angleDefect() const
{
    double defect = 2 * M_PI;
    
    HalfEdgeCIter h = he;
    do {
        Eigen::Vector3d u = h->next->vertex->position - h->vertex->position;
        Eigen::Vector3d v = h->next->next->vertex->position - h->vertex->position;
        
        double theta = acos(u.dot(v) / (u.norm() * v.norm()));
        defect -= theta;
        
        h = h->flip->next;
        
    } while (h != he);
    
    return defect;
}