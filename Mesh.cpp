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

double Mesh::computeGaussCurvature(Eigen::VectorXd& K)
{
    double maxGauss = -INFINITY;
    for (size_t i = 0; i < vertices.size(); i++) {
        K(i) = vertices[i].angleDefect() / vertices[i].dualArea();
        
        if (maxGauss < fabs(K(i))) maxGauss = fabs(K(i));
    }
    
    return maxGauss;
}

void Mesh::buildLaplacian(Eigen::SparseMatrix<double>& L) const
{
    L.resize((int)vertices.size(), (int)vertices.size());
    
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        
        HalfEdgeCIter he = v->he;
        double dualArea = v->dualArea();
        double sumCoefficients = 0.0;
        do {
            // (cotA + cotB) / 2A
            double coefficient = 0.5 * (he->cotan() + he->flip->cotan()) / dualArea;
            sumCoefficients += coefficient;
            
            L.insert(v->index, he->flip->vertex->index) = coefficient;
            
            he = he->flip->next;
        } while (he != v->he);
        
        L.insert(v->index, v->index) = -sumCoefficients;
    }
    
    L.makeCompressed();
}

double Mesh::computeMeanCurvature(Eigen::VectorXd& H)
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
    double maxMean = -INFINITY;
    for (size_t i = 0; i < vertices.size(); i++) {
        H(i) = 0.5 * x.row(i).norm();
        
        if (maxMean < H(i)) maxMean = H(i);
    }
    
    return maxMean;
}

void Mesh::computeCurvatures()
{
    int v = (int)vertices.size();
    
    Eigen::VectorXd K(v);
    double maxGauss = computeGaussCurvature(K);
    
    Eigen::VectorXd H(v);
    double maxMean = computeMeanCurvature(H);
    
    // compute principal curvatures and normalize gauss, mean curvature
    for (int i = 0; i < v; i++) {
        double dis = sqrt(H(i)*H(i) - K(i));
        
        vertices[i].k1 = H(i) + dis;
        vertices[i].k2 = H(i) - dis;
        
        vertices[i].gaussCurvature = K(i) / maxGauss;
        vertices[i].meanCurvature = H(i) / maxMean;
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
