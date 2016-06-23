#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Mesh.h"

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.;
double x = 0;
double y = 0;
double z = -2.5;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/halfedge-mesh/bunny.obj";

Mesh mesh;
bool success = true;
int curvature = 0; // 0: Gaussian, 1: Mean, 2: Normal
double theta = 0;
std::vector<double> normalCurvatures;

void normalizeNormalCurvatures()
{
    int v = (int)mesh.vertices.size();
    normalCurvatures.reserve(v);
    
    double maxNormal = -INFINITY;
    for (int i = 0; i < v; i++) {
        normalCurvatures[mesh.vertices[i].index] = mesh.vertices[i].normalCurvature(theta);
        
        if (maxNormal < fabs(normalCurvatures[mesh.vertices[i].index])) {
            maxNormal = fabs(normalCurvatures[mesh.vertices[i].index]);
        }
    }
    
    for (int i = 0; i < v; i++) {
        normalCurvatures[i] /= maxNormal;
    }
}

void printInstructions()
{
    std::cerr << "→/←: toggle between gaussian, mean and normal curvature\n"
              << "o/p: decrement/increment theta value for normal curvature between 0 to 180\n"
              << "↑/↓: move in/out\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "escape: exit program\n"
              << std::endl;
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
}

void draw()
{
    glBegin(GL_LINES);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        
        int s = 3;
        Eigen::Vector3d v = e->he->vertex->position;
        Eigen::Vector3d u = e->he->flip->vertex->position - v;
        u.normalize();
        double dl = e->length() / (double)s;
        
        double c = 0.0;
        double c2 = 0.0;
        if (curvature == 0) {
            c = e->he->vertex->gaussCurvature;
            c2 = e->he->flip->vertex->gaussCurvature;
            
        } else if (curvature == 1) {
            c = e->he->vertex->meanCurvature;
            c2 = e->he->flip->vertex->meanCurvature;
        
        } else {
            c = normalCurvatures[e->he->vertex->index];
            c2 = normalCurvatures[e->he->flip->vertex->index];
        }
        double dc = (c2 - c) / (double)s;
        
        for (int i = 0; i < s; i++) {
            if (c < 0) glColor4f(0.0, 0.0, fabs(c), 0.6);
            else glColor4f(c, 0.0, 0.0, 0.6);
            glVertex3d(v.x(), v.y(), v.z());
            
            c += dc;
            if (c < 0) glColor4f(0.0, 0.0, fabs(c), 0.6);
            else glColor4f(c, 0.0, 0.0, 0.6);
            v += u * dl;
            glVertex3d(v.x(), v.y(), v.z());
        }
    }
    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 2.0, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        if (curvature == 2) {
            normalizeNormalCurvatures();
        }
        
        draw();
    }

    glutSwapBuffers();
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case 'o':
            if (curvature == 2) {
                theta -= 10.0;
                if (theta < 0.0) theta = 180.0;
                std::string title = "Normal Curvature, ø = " + std::to_string(theta);
                glutSetWindowTitle(title.c_str());
            }
            break;
        case 'p':
            if (curvature == 2) {
                theta += 10.0;
                if (theta > 180.0) theta = 0.0;
                std::string title = "Normal Curvature, ø = " + std::to_string(theta);
                glutSetWindowTitle(title.c_str());
            }
            break;
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z += 0.03;
            break;
        case GLUT_KEY_DOWN:
            z -= 0.03;
            break;
        case GLUT_KEY_LEFT:
            curvature --;
            if (curvature < 0) curvature = 2;
            break;
        case GLUT_KEY_RIGHT:
            curvature ++;
            if (curvature == 3) curvature = 0;
            break;
    }
    
    if (curvature == 0) {
        glutSetWindowTitle("Gaussian Curvature");
        
    } else if (curvature == 1) {
        glutSetWindowTitle("Mean Curvature");
    
    } else {
        std::string title = "Normal Curvature, ø = " + std::to_string(theta);
        glutSetWindowTitle(title.c_str());
    }
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    mesh.computeCurvatures();
    
    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    glutCreateWindow("Gaussian Curvature");
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}
