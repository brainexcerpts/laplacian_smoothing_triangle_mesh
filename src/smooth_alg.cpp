#include "smooth_alg.hpp"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Eigen/QR>
#include <Eigen/LU>
#include <Eigen/SparseLU>

// -----------------------------------------------------------------------------

typedef Eigen::Triplet<double, int> Triplet;
/// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<double> Sparse_mat;

// -----------------------------------------------------------------------------

/// @return A sparse representation of the normalized Laplacian
/// list[ith_row][list of columns] = Triplet(ith_row, jth_column, matrix value)
static
std::vector<std::vector<Triplet>>
get_normalized_laplacian(const std::vector< Vec3 >& vertices,
                         const std::vector< std::vector<int> >& edges )
{
    unsigned nv = unsigned(vertices.size());
    std::vector<std::vector<Triplet>> mat_elemts(nv);
    for(int i = 0; i < nv; ++i)
        mat_elemts[i].reserve(10);

    for(int i = 0; i < nv; ++i)
    {
        const Vec3 c_pos = vertices[i];

        //get laplacian
        double sum = 0.;
        int nb_edges = edges[i].size();
        for(int e = 0; e < nb_edges; ++e)
        {
            int next_edge = (e + 1           ) % nb_edges;
            int prev_edge = (e + nb_edges - 1) % nb_edges;

            Vec3 v1 = c_pos                 - vertices[edges[i][prev_edge]];
            Vec3 v2 = vertices[edges[i][e]] - vertices[edges[i][prev_edge]];
            Vec3 v3 = c_pos                 - vertices[edges[i][next_edge]];
            Vec3 v4 = vertices[edges[i][e]] - vertices[edges[i][next_edge]];

            double cotan1 = (v1.dot(v2)) / (1e-6 + (v1.cross(v2)).norm() );
            double cotan2 = (v3.dot(v4)) / (1e-6 + (v3.cross(v4)).norm() );

            double w = (cotan1 + cotan2)*0.5;
            sum += w;
            mat_elemts[i].push_back( Triplet(i, edges[i][e], w) );
        }

        for( Triplet& t : mat_elemts[i] )
            t = Triplet( t.row(), t.col(), t.value() / sum);

        mat_elemts[i].push_back( Triplet(i, i, -1.0) );
    }
    return mat_elemts;
}

// -----------------------------------------------------------------------------

static
std::vector< std::vector<float> >
get_cotan_weights(const std::vector< Vec3 >& vertices,
                  const std::vector< std::vector<int> >& edges)
{
    unsigned nv = unsigned(vertices.size());
    std::vector< std::vector<float> > weights(nv);

    for(int i = 0; i < nv; ++i)
    {
        const Vec3 c_pos = vertices[i];
        double sum = 0.;
        int nb_edges = edges[i].size();
        weights[i].resize( nb_edges );
        for(int e = 0; e < nb_edges; ++e)
        {
            int next_edge = (e + 1           ) % nb_edges;
            int prev_edge = (e + nb_edges - 1) % nb_edges;

            Vec3 v1 = c_pos                 - vertices[edges[i][prev_edge]];
            Vec3 v2 = vertices[edges[i][e]] - vertices[edges[i][prev_edge]];
            Vec3 v3 = c_pos                 - vertices[edges[i][next_edge]];
            Vec3 v4 = vertices[edges[i][e]] - vertices[edges[i][next_edge]];

            double cotan1 = (v1.dot(v2)) / (1e-6 + (v1.cross(v2)).norm() );
            double cotan2 = (v3.dot(v4)) / (1e-6 + (v3.cross(v4)).norm() );

            double w = (cotan1 + cotan2)*0.5;
            weights[i][e] = w;
        }
    }
    return weights;
}

//------------------------------------------------------------------------------

std::vector< Vec3 >
smooth_iterative(const std::vector< Vec3 >& in_vertices,
                 const std::vector< std::vector<int> >& edges,
                 int nb_iter,
                 float alpha)
{
    unsigned nb_vertices = unsigned(in_vertices.size());
    std::vector< std::vector<float> > cotan_weights = get_cotan_weights(in_vertices, edges);

    std::vector< Vec3 > buffer_vertices( nb_vertices );
    std::vector< Vec3 > source = in_vertices;

    if(nb_iter == 0){
        return in_vertices;
    }

    Vec3* src_vertices = source.data();
    Vec3* dst_vertices = buffer_vertices.data();
    for(int k = 0; k < nb_iter; k++)
    {
        for( int i = 0; i < nb_vertices; i++)
        {
#if 0
            if( _topo->is_vert_on_side( i ) ){
                dst_vertices[i] = src_vertices[i];
                continue;
            }
#endif

            Vec3 cog(0.f);
            float sum = 0.f;
            size_t nb_neighs = edges[i].size();
            for(size_t n = 0; n < nb_neighs; n++)
            {
                Vert_idx neigh = edges[i][n];
                float w = cotan_weights[i][n];
                cog += src_vertices[neigh] * w;
                sum += w;
            }
            float t = alpha;
            dst_vertices[i] = (cog / sum) * t + src_vertices[i] * (1.f - t);

        }

        std::swap(dst_vertices, src_vertices);
    }

    return (nb_iter%2 == 1) ? buffer_vertices : source;
}

//------------------------------------------------------------------------------

//explicit
std::vector< Vec3 >
explicit_laplacian_smoothing(const std::vector< Vec3 >& in_vertices,
                             const std::vector< std::vector<int> >& edges,
                             int nb_iter,
                             float alpha)
{
    unsigned nb_vertices = unsigned(in_vertices.size());

    std::vector<Eigen::VectorXd> xyz;
    std::vector<Eigen::VectorXd> rhs;

    xyz.resize(3, Eigen::VectorXd::Zero(nb_vertices));
    rhs.resize(3, Eigen::VectorXd::Zero(nb_vertices));

    for(int i = 0; i < nb_vertices; ++i)
    {
        Vec3 pos = in_vertices[i];
        xyz[0][i] = pos.x;
        xyz[1][i] = pos.y;
        xyz[2][i] = pos.z;
    }

    // Build laplacian
    std::vector<std::vector<Triplet>> mat_elemts = get_normalized_laplacian(in_vertices, edges);
    Eigen::SparseMatrix<double> L(nb_vertices, nb_vertices);
    std::vector<Triplet> triplets;
    triplets.reserve(nb_vertices * 10);
    for( const std::vector<Triplet>& row : mat_elemts)
        for( const Triplet& elt : row )
            triplets.push_back( elt );

    L.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> I = Eigen::MatrixXd::Identity(nb_vertices, nb_vertices).sparseView();
    L = I + L*alpha;

    //L = L*L*L*L*L*L*L*L*L*L;
    rhs = xyz;
    for(int n = 0; n < nb_iter; n++){
        for(int k = 0; k < 3; k++){
            xyz[k] = (L * xyz[k]);
        }
    }

    std::vector< Vec3 > out_verts(nb_vertices);
    for(int i = 0; i < nb_vertices; ++i){
        Vec3 v;
        v.x = xyz[0][i];
        v.y = xyz[1][i];
        v.z = xyz[2][i];
        out_verts[i] = v;
    }

    return out_verts;
}

//------------------------------------------------------------------------------

//implicit
std::vector< Vec3 >
implicit_laplacian_smoothing(const std::vector< Vec3 >& in_vertices,
                             const std::vector< std::vector<int> >& edges,
                             int nb_iter,
                             float alpha)
{
    unsigned nb_vertices = unsigned(in_vertices.size());

    std::vector<Eigen::VectorXd> xyz;
    std::vector<Eigen::VectorXd> rhs;

    xyz.resize(3, Eigen::VectorXd::Zero(nb_vertices));
    rhs.resize(3, Eigen::VectorXd::Zero(nb_vertices));

    for(int i = 0; i < nb_vertices; ++i)
    {
        Vec3 pos = in_vertices[i];
        rhs[0][i] = pos.x;
        rhs[1][i] = pos.y;
        rhs[2][i] = pos.z;
    }


    // Build laplacian
    std::vector<std::vector<Triplet>> mat_elemts = get_normalized_laplacian(in_vertices, edges);
    Eigen::SparseMatrix<double> L(nb_vertices, nb_vertices);
    std::vector<Triplet> triplets;
    triplets.reserve(nb_vertices * 10);
    for( const std::vector<Triplet>& row : mat_elemts)
        for( const Triplet& elt : row )
            triplets.push_back( elt );

    L.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseMatrix<double> I = Eigen::MatrixXd::Identity(nb_vertices, nb_vertices).sparseView();
    L = I - L*alpha;

    L = L*L*L;

    // Solve for x, y, z
    Eigen::SparseLU<Sparse_mat> solver;
    solver.compute( L );

    for(int k = 0; k < 3; k++){
        xyz[k] = solver.solve(rhs[k]);
    }

    std::vector< Vec3 > out_verts(nb_vertices);
    for(int i = 0; i < nb_vertices; ++i){
        Vec3 v;
        v.x = xyz[0][i];
        v.y = xyz[1][i];
        v.z = xyz[2][i];
        out_verts[i] = v;
    }

    return out_verts;
}


