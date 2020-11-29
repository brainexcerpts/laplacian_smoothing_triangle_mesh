#include <GL/glut.h>

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#include "mesh.hpp"
#include "topology/vertex_to_face.hpp"
#include "topology/vertex_to_1st_ring_vertices.hpp"
#include "smooth_alg.hpp"


// compatibility with original GLUT
#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

// =============================================================================

float _table_angle = 0.0f;
bool _3d_view = false;
static int _win_number;
Mesh* _mesh;
Vertex_to_1st_ring_vertices _first_ring;
std::vector<Vec3> _rest_pose_vertices;

// =============================================================================

void compute_normals(Mesh& mesh)
{
    unsigned nb_vertices = 0;
    mesh._normals.assign( nb_vertices, Vec3(0.0f));
    for(unsigned i = 0; i < mesh._triangles.size(); i++ ) {
        const Tri_face& tri = mesh._triangles[i];
        // Compute normals:
        Vec3 v1 = mesh._vertices[tri.b] - mesh._vertices[tri.a];
        Vec3 v2 = mesh._vertices[tri.c] - mesh._vertices[tri.a];

        Vec3 n = v1.cross( v2 );
        n.normalize();

        mesh._normals[ tri.a ] += n;
        mesh._normals[ tri.b ] += n;
        mesh._normals[ tri.c ] += n;
    }
}

// -----------------------------------------------------------------------------

/// @brief Load mesh from a file
static Mesh* build_mesh(const char* file_name)
{
    Mesh* ptr = new Mesh();
    Mesh& mesh = *ptr;

    std::ifstream file( file_name );

    if( !file.is_open() ){
        std::cerr << "Can't open file: " << file_name << std::endl;
        return nullptr;
    }

    std::string format;
    unsigned nb_vertices = 0;
    unsigned nb_faces = 0;
    int nil;

    file >> format;
    file >> nb_vertices >> nb_faces >> nil;

    mesh._vertices.resize( nb_vertices );
    for (unsigned i=0; i < nb_vertices; i++ ){
        Vec3 p;
        file >> p.x >> p.y >> p.z;
        mesh._vertices[i] = p;
    }

    mesh._triangles.resize( nb_faces );
    mesh._colors.assign( nb_vertices, Vec3(0.0f));

    mesh._normals.assign( nb_vertices, Vec3(0.0f));
    for(unsigned i = 0; i < nb_faces; i++ )
    {
        int nb_verts_face;
        file >> nb_verts_face;

        if(nb_verts_face != 3) {
            std::cerr << "We only handle triangle faces" << std::endl;
            continue;
        }

        Tri_face& tri = mesh._triangles[i];
        file >> tri.a >> tri.b >> tri.c;
    }
    compute_normals(mesh);
    file.close();
    return ptr;
}

// ----------

void draw_mesh(const Mesh& mesh) {
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, mesh._vertices.data());
    glNormalPointer(GL_FLOAT, 0, mesh._normals.data());
    glColorPointer(3, GL_FLOAT, 0, mesh._colors.data());
    glDrawElements(GL_TRIANGLES, mesh._triangles.size()*3, GL_UNSIGNED_INT, mesh._triangles.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
}

// -----------------------------------------------------------------------------

void key_stroke (unsigned char c, int mouseX, int mouseY) {
    static bool wires  = false;

    switch (c) {
    case 27 :
        glFinish();
        glutDestroyWindow(_win_number);
        exit (0);
        break;
    case 'w' :
        if(!wires) {
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_LIGHTING);
            glutPostRedisplay();
            wires = true;
        }else {
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_LIGHTING);
            glutPostRedisplay();
            wires = false;
        }
        break;
    }
}

// -----------------------------------------------------------------------------

void mouse_keys (int button, int state, int x, int y)
{
    if(button == GLUT_WHEEL_UP)
        _table_angle += 1.0f;
    if(button == GLUT_WHEEL_DOWN)
        _table_angle -= 1.0f;
}

// -----------------------------------------------------------------------------

void smooth_mesh(Mesh& mesh, std::vector<Vec3> in_vertices)
{
    static float iter = 2;
    //mesh._vertices = smooth_iterative(in_vertices, _first_ring._rings_per_vertex, int(iter));
    //mesh._vertices = explicit_laplacian_smoothing(in_vertices, _first_ring._rings_per_vertex, int(iter), 1.0); // Unstable
    //mesh._vertices = explicit_laplacian_smoothing(in_vertices, _first_ring._rings_per_vertex, int(iter), 0.5); // more stable
    mesh._vertices = implicit_laplacian_smoothing(in_vertices, _first_ring._rings_per_vertex, int(iter), 200.0f); // Unconditionally stable

    if(_3d_view){
        compute_normals(*_mesh);
        for(unsigned v = 0; v < mesh.nb_vertices(); ++v)
            mesh._colors[v] = Vec3(0.6f, 0.2f, 0.f);//(mesh._normals[v]+1.0f)*0.5f;
    }

}

// -----------------------------------------------------------------------------

void display(void)
{
    glClearColor ( 0.93,  0.93,  0.93,  0.93);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glViewport(0,0,900,900);

    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60., 1., 0.5, 100.);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    if(_3d_view){
        GLfloat lightPos0[] = {0.0f, 0.0f, 0.0f, 1.0f};
        glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

        glTranslatef(0.0, 0.0, -1.5);
        glRotatef(45.0f, -1.0f, 0.0f, 0.0f);
        glRotatef(_table_angle, -1.0f, 0.0f, 0.0f);

        static float angle = 0.0f;
        angle = fmodf(angle+0.3f, 360.f);
        //glRotatef(angle, 1.0f, 0.0f, 0.0f);
        glRotatef(angle, 0.0f, 0.0f, 1.0f);
        glutPostRedisplay();
    }
    else
        glTranslatef(0.0, 0.0, -1.0);
    float s = 0.5f;
    glScalef(s, s, s);

    draw_mesh(*_mesh);

    glutSwapBuffers();
    glFlush ();
}

// -----------------------------------------------------------------------------

void load_mesh()
{
    _mesh = build_mesh("samples/buddha.off");
    //g_mesh = build_mesh("samples/donut.off");
    Mesh& mesh = *_mesh;
    _rest_pose_vertices = mesh._vertices;
    // Compute first ring
    Vertex_to_face v_to_face;
    v_to_face.compute( mesh );
    _first_ring.compute(mesh, v_to_face );
}

// -----------------------------------------------------------------------------

int main (int argc, char** argv)
{

#ifndef NDEBUG
   std::cout << "debug" << std::endl;
#else
    std::cout << "release" << std::endl;
#endif
    load_mesh();
    smooth_mesh(*_mesh, _rest_pose_vertices);

    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize (900, 900);
    glutInitWindowPosition (240, 212);
    _win_number = glutCreateWindow (argv[0]);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLineWidth(2.2f);
    //glEnable(GL_LINE_SMOOTH);
    //glLightfv(GL_LIGHT0, GL_AMBIENT,  {});
    //glLightfv(GL_LIGHT0, GL_DIFFUSE,  Global::_light0_diffuse);
    //GLfloat lightp[] = {1.0f,1.0f,1.0f,1.0f};
    //glLightfv(GL_LIGHT0, GL_SPECULAR, lightp);

    glutKeyboardFunc(key_stroke);
    glutMouseFunc(mouse_keys);
    glutDisplayFunc(display);
    glutMainLoop();

    return (0);
}

