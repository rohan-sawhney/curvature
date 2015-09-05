#ifndef MESH_H
#define MESH_H

#include "Types.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"
#include <Eigen/SparseCore>

class Mesh {
public:
    // default constructor
    Mesh();
    
    // copy constructor
    Mesh(const Mesh& mesh);
        
    // read mesh from file
    bool read(const std::string& fileName);
    
    // write mesh to file
    bool write(const std::string& fileName) const;
    
    // computes principal, gaussian and mean curvatures
    void computeCurvatures();
    
    // member variables
    std::vector<HalfEdge> halfEdges;
    std::vector<Vertex> vertices;
    std::vector<Eigen::Vector3d> uvs;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<HalfEdgeIter> boundaries;

private:
    // computes gaussian curvature per vertex
    double computeGaussCurvature(Eigen::VectorXd& K);
    
    // builds Laplace Beltrami operator
    void buildLaplacian(Eigen::SparseMatrix<double>& L) const;
    
    // computes mean curvature per vertex
    double computeMeanCurvature(Eigen::VectorXd& H);
    
    // center mesh about origin and rescale to unit radius
    void normalize();
};

#endif