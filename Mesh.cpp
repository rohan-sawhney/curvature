#include "Mesh.h"
#include "MeshIO.h"

Mesh::Mesh()
{
    
}

Mesh::Mesh(const Mesh& mesh)
{
    *this = mesh;
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
    }
    
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    return false;
}

void Mesh::computeGaussCurvature(Eigen::VectorXd& K)
{
    for (size_t i = 0; i < vertices.size(); i++) {
        K(i) = vertices[i].angleDefect() / vertices[i].dualArea();
    }
}

void Mesh::buildLaplacian(Eigen::SparseMatrix<double>& L)
{
    L.resize((int)vertices.size(), (int)vertices.size());
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double dualArea = v->dualArea();
        do {
            // (cotA + cotB) / 2A
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan()) / dualArea;
            
            L.insert(v->index, he->flip->vertex->index) = coefficient;
            L.coeffRef(v->index, v->index) -= coefficient;

            he = he->flip->next;
        } while (he != v->he);
    }
    
    L.makeCompressed();
}

void Mesh::computeMeanCurvature(Eigen::VectorXd& H)
{
    Eigen::SparseMatrix<double> L;
    buildLaplacian(L);
    
    Eigen::MatrixXd x;
    x.resize((int)vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); i++) {
        x.row(i) = vertices[i].position;
    }
    x = L * x;

    // set absolute mean curvature
    for (size_t i = 0; i < vertices.size(); i++) {
        H(i) = 0.5 * x.row(i).norm();
    }
}

void Mesh::computeCurvatures()
{
    int v = (int)vertices.size();
    
    Eigen::VectorXd K(v);
    computeGaussCurvature(K);
    
    Eigen::VectorXd H(v);
    computeMeanCurvature(H);
    
    Eigen::VectorXd k1(v);
    Eigen::VectorXd k2(v);
    
    // compute principal curvatures
    double maxGauss = -INFINITY;
    double maxMean = -INFINITY;
    for (int i = 0; i < v; i++) {
        double dis = sqrt(H(i)*H(i) - K(i));
        
        k1(i) = H(i) + dis;
        k2(i) = H(i) - dis;
        
        vertices[i].gaussCurvature = k1(i) * k2(i);
        if (maxGauss < fabs(vertices[i].gaussCurvature)) maxGauss = fabs(vertices[i].gaussCurvature);
        
        vertices[i].meanCurvature = (k1(i) + k2(i)) / 2.0;
        if (maxMean < fabs(vertices[i].meanCurvature)) maxMean = fabs(vertices[i].meanCurvature);
    }
    
    // normalize curvatures
    for (int i = 0; i < v; i++) {
        vertices[i].gaussCurvature /= maxGauss;
        vertices[i].meanCurvature /= maxMean;
    }
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
    }
    
    // determine radius
    double rMax = 0;
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
