#include <Eigen/Dense>
#include <vector>

// Define Face class. A face will have a set of vertices,
// an owner cell, a neighbour cell, an area-scaled centroid,
// and a normal (Sf) pointing out from its owner cell and scaled
// by the face area
class Face {
 public:
  Eigen::ArrayXi vertices;  // Vertices making up face
  int ownerCell;            // Owner cell of face
  int neighbourCell;        // Neighbour cell of face
  Eigen::Vector3d Sf;  // Normal vector to face, pointing out from owner and
                       // scaled by face area
  Eigen::Vector3d centroid;  // areal centroid of face
  double f_x;                // Face interpolation factor
  double delta;  // delta (1/distance between owner and nbr centroids)
  int boundary;  // Index of adjoining boundary patch. If not at boundary it is
                 // set as -1
  bool empty_face;  // Boolean to indicate if face is at empty boundary (i.e. z
                    // boundary in a 2D simulation)
  double flux;      // Flux across the face
};
