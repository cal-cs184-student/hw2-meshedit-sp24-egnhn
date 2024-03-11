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
      HalfedgeIter he0 = e0->halfedge();
      HalfedgeIter he1 = he0->next();
      HalfedgeIter he2 = he1->next();
      HalfedgeIter he3 = he0->twin();
      HalfedgeIter he4 = he3->next();
      HalfedgeIter he5 = he4->next();
      HalfedgeIter he6 = he5->twin();
      HalfedgeIter he7 = he4->twin();
      HalfedgeIter he8 = he2->twin();
      HalfedgeIter he9 = he1->twin();

      VertexIter v0 = he0->vertex();
      VertexIter v1 = he3->vertex();
      VertexIter v2 = he2->vertex();
      VertexIter v3 = he5->vertex();

      EdgeIter e1 = he1->edge();
      EdgeIter e2 = he2->edge();
      EdgeIter e3 = he4->edge();
      EdgeIter e4 = he5->edge();

      FaceIter f0 = he0->face();
      FaceIter f1 = he3->face();

      v0->halfedge() = he5;
      v1->halfedge() = he2;
      v2->halfedge() = he4;
      v3->halfedge() = he1;

      e0->halfedge() = he0;
      e1->halfedge() = he2;
      e2->halfedge() = he4;
      e3->halfedge() = he5;
      e4->halfedge() = he1;

      f0->halfedge() = he0;
      f1->halfedge() = he3;

      he0->setNeighbors(he1, he3, v2, e0, f0);
      he1->setNeighbors(he2, he6, v3, e4, f0);
      he2->setNeighbors(he0, he9, v1, e1, f0);
      he3->setNeighbors(he4, he0, v3, e0, f1);
      he4->setNeighbors(he5, he8, v2, e2, f1);
      he5->setNeighbors(he3, he7, v0, e3, f1);
      he6->setNeighbors(he6->next(), he1, v1, e4, he6->face());
      he7->setNeighbors(he7->next(), he5, v3, e3, he7->face());
      he8->setNeighbors(he8->next(), he4, v0, e2, he8->face());
      he9->setNeighbors(he9->next(), he2, v2, e1, he9->face());
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

    HalfedgeIter n0 = newHalfedge();
    HalfedgeIter n1 = newHalfedge();
    HalfedgeIter n2 = newHalfedge();
    HalfedgeIter n3 = newHalfedge();
    HalfedgeIter n4 = newHalfedge();
    HalfedgeIter n5 = newHalfedge();

    HalfedgeIter he0 = e0->halfedge();
    HalfedgeIter he1 = he0->next();
    HalfedgeIter he2 = he1->next();
    HalfedgeIter he3 = he0->twin();
    HalfedgeIter he4 = he3->next();
    HalfedgeIter he5 = he4->next();

    VertexIter v0 = he0->vertex();
    VertexIter v1 = he3->vertex();
    VertexIter v2 = he2->vertex();
    VertexIter v3 = he5->vertex();

    FaceIter f0 = he0->face();
    FaceIter f1 = newFace();
    FaceIter f2 = newFace();
    FaceIter f3 = newFace();

    VertexIter midPt = newVertex();
    midPt->position = (v0->position + v1->position) / 2.0;
    midPt->halfedge() = n0;

    he0->setNeighbors(n2, n3, v0, e0, f0);
    he3->setNeighbors(n4, n0, v1, e1, f3);
    he5->setNeighbors(he3, he5->twin(), v3, he5->edge(), f3);
    he1->setNeighbors(n1, he1->twin(), v1, he1->edge(), f1);
    he2->setNeighbors(he0, he2->twin(), v2, he2->edge(), f0);
    he4->setNeighbors(n5, he4->twin(), v0, he4->edge(), f2);

    n0->setNeighbors(he1, he3, midPt, e1, f1);
    n1->setNeighbors(n0, n2, v2, e2, f1);
    n2->setNeighbors(he2, n1, midPt, e2, f0);
    n3->setNeighbors(he4, he0, midPt, e0, f2);
    n4->setNeighbors(he5, n5, midPt, e3, f3);
    n5->setNeighbors(n3, n4, v3, e3, f2);

    f0->halfedge() = he0;
    f1->halfedge() = n0;
    f2->halfedge() = n3;
    f3->halfedge() = he3;
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
