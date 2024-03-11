#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    std::vector<Vector2D> nextPts; 
    for (size_t i = 0; i < points.size() - 1; ++i) { 
      Vector2D p1 = points[i]; 
      Vector2D p2 = points[i + 1]; 
      Vector2D nextPt = p1 * (1-t) + p2 * t; 
      nextPts.push_back(nextPt); 
    } 
    return nextPts;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    std::vector<Vector3D> nextPts; 
    for (size_t i = 0; i < points.size() - 1; ++i) { 
      Vector3D p1 = points[i]; 
      Vector3D p2 = points[i + 1]; 
      Vector3D nextPt = p1 * (1-t) + p2 * t;  
      nextPts.push_back(nextPt); 
    } 
    return nextPts;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    std::vector<Vector3D> nextPts(points);
    for (size_t i = 0; i < points.size() - 1; ++i) {
      nextPts = evaluateStep(nextPts, t);
    }
    return nextPts[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    std::vector<Vector3D> curvePts;
    for (size_t i = 0; i < controlPoints.size(); ++i) {
      std::vector<Vector3D> rowPts = controlPoints[i];
      Vector3D ptsOnRow = evaluate1D(rowPts, u);
      curvePts.push_back(ptsOnRow);
    }
    return evaluate1D(curvePts, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    double totalArea = 0.0;
    Vector3D normalSum(0, 0, 0);
    HalfedgeCIter halfEdge = halfedge();
    do {
      if (!halfEdge->face()->isBoundary()){
        Vector3D p0 = position;
        Vector3D p1 = halfEdge->next()->vertex()->position;
        Vector3D p2 = halfEdge->next()->next()->vertex()->position;
        Vector3D e1 = p1 - p0;
        Vector3D e2 = p2 - p0;
        double area = 0.5 * cross(e1, e2).norm();
        normalSum += area * halfEdge->face()->normal();
        totalArea += area;
      }
      halfEdge = halfEdge->twin()->next();
    } while (halfEdge != halfedge()); 

    if (totalArea != 0.0) {
      return normalSum / totalArea;
    } 
    return Vector3D(0, 0, 1);
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (!e0->isBoundary()) {
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h5->twin();
      HalfedgeIter h7 = h4->twin();
      HalfedgeIter h8 = h2->twin();
      HalfedgeIter h9 = h1->twin();
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h2->vertex();
      VertexIter v3 = h5->vertex();
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();
      v0->halfedge() = h5;
      v1->halfedge() = h2;
      v2->halfedge() = h4;
      v3->halfedge() = h1;
      e0->halfedge() = h0;
      e1->halfedge() = h2;
      e2->halfedge() = h4;
      e3->halfedge() = h5;
      e4->halfedge() = h1;
      f0->halfedge() = h0;
      f1->halfedge() = h3;
      h0->setNeighbors(h1, h3, v2, e0, f0);
      h1->setNeighbors(h2, h6, v3, e4, f0);
      h2->setNeighbors(h0, h9, v1, e1, f0);
      h3->setNeighbors(h4, h0, v3, e0, f1);
      h4->setNeighbors(h5, h8, v2, e2, f1);
      h5->setNeighbors(h3, h7, v0, e3, f1);
      h6->setNeighbors(h6->next(), h1, v1, e4, h6->face());
      h7->setNeighbors(h7->next(), h5, v3, e3, h7->face());
      h8->setNeighbors(h8->next(), h4, v0, e2, h8->face());
      h9->setNeighbors(h9->next(), h2, v2, e1, h9->face());
    }
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if (e0->halfedge()->isBoundary()) {
      return e0->halfedge()->vertex();
    }
    EdgeIter e1 = newEdge();
    EdgeIter e2 = newEdge();
    EdgeIter e3 = newEdge();
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter n0 = newHalfedge();
    HalfedgeIter n1 = newHalfedge();
    HalfedgeIter n2 = newHalfedge();
    HalfedgeIter n3 = newHalfedge();
    HalfedgeIter n4 = newHalfedge();
    HalfedgeIter n5 = newHalfedge();
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();
    FaceIter f0 = h0->face();
    FaceIter f1 = newFace();
    FaceIter f2 = newFace();
    FaceIter f3 = newFace();
    VertexIter midPt = newVertex();
    midPt->position = (v0->position + v1->position) / 2.0;
    midPt->halfedge() = n0;
    h0->setNeighbors(n2, n3, v0, e0, f0);
    h3->setNeighbors(n4, n0, v1, e1, f3);
    h5->setNeighbors(h3, h5->twin(), v3, h5->edge(), f3);
    h1->setNeighbors(n1, h1->twin(), v1, h1->edge(), f1);
    h2->setNeighbors(h0, h2->twin(), v2, h2->edge(), f0);
    h4->setNeighbors(n5, h4->twin(), v0, h4->edge(), f2);
    n0->setNeighbors(h1, h3, midPt, e1, f1);
    n1->setNeighbors(n0, n2, v2, e2, f1);
    n2->setNeighbors(h2, n1, midPt, e2, f0);
    n3->setNeighbors(h4, h0, midPt, e0, f2);
    n4->setNeighbors(h5, n5, midPt, e3, f3);
    n5->setNeighbors(n3, n4, v3, e3, f2);
    f0->halfedge() = h0;
    f1->halfedge() = n0;
    f2->halfedge() = n3;
    f3->halfedge() = h3;
    e1->halfedge() = n0;
    e2->halfedge() = n1;
    e3->halfedge() = n4;
    e0->isNew = false;
    e1->isNew = false;
    e2->isNew = true;
    e3->isNew = true;
    return midPt;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    for (VertexIter e = mesh.verticesBegin(); e != mesh.verticesEnd(); e++){

      e->isNew = false;
      float deg = e->degree();
      float div = 3 / (8 * deg);

      if (deg == 3){
        div = 3 / 16;
      }

      Vector3D neighborSum(0, 0, 0);
      HalfedgeCIter halfEdge = e->halfedge();

      do {
        halfEdge = halfEdge->twin();
        neighborSum += halfEdge->vertex()->position;
        halfEdge = halfEdge->next();
      } while (halfEdge != e->halfedge());

      e->newPosition = e->position * (1 - deg * div) + div * neighborSum;
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++){

      Vector3D v0 = e->halfedge()->vertex()->position;
      Vector3D v1 = e->halfedge()->twin()->vertex()->position;
      Vector3D v2 = e->halfedge()->next()->next()->vertex()->position;
      Vector3D v3 = e->halfedge()->twin()->next()->next()->vertex()->position;

      e->newPosition = 3 / 8 * (v0 + v1) + 1 / 8 * (v2 + v3);
      e->isNew = 0;

    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    int numEdges = mesh.nEdges();
    int i = 0;

    for (EdgeIter e = mesh.edgesBegin(); i < numEdges; i++) {

      if (e->isNew == 0) {
      VertexIter numSEdge = mesh.splitEdge(e);
      numSEdge->isNew = 1;
      numSEdge->position = e->newPosition;
      }

      e++;
    }
    // 4. Flip any new edge that connects an old and new vertex.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      int ctr = 0;
      if (e->isNew) {
        HalfedgeIter h = e->halfedge();

        if (h->vertex()->isNew) {
          ctr++;
        }

        if (h->twin()->vertex()->isNew) {
          ctr++;
        }

        if (ctr == 1) {
          e = mesh.flipEdge(e);
          e->isNew = 0;
        }

      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {

      if (v->isNew == 0) {
        v->position = v->newPosition;
      } else {
        v->isNew = 0;
      }

    }

  }
}
